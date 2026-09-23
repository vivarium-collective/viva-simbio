"""Generate a self-contained HTML demo report for viva-simbio.

Runs three simbio CRN configurations through process-bigraph composites,
collects time-series snapshots, and renders an interactive report with Plotly
charts, metrics cards, a bigraph architecture diagram, and a collapsible PBG
document tree. Opens the report in the default browser when done.
"""

from __future__ import annotations

import base64
import json
import os
import time

from process_bigraph import Composite, allocate_core, gather_emitter_results

from viva_simbio.composites.crn import brusselator, lotka_volterra, repressilator

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "report.html")
TIMEOUT_S = 120.0

CONFIGS = [
    {
        "id": "brusselator",
        "title": "Brusselator",
        "subtitle": "a chemical clock",
        "description": "The Prigogine oscillator: an autocatalytic loop drives "
                       "sustained limit-cycle oscillations.",
        "accent": "#6366f1",
        "builder": brusselator,
        "params": {"k2": 3.0, "k3": 1.0, "interval": 0.25},
        "total_time": 25.0,
    },
    {
        "id": "lotka",
        "title": "Lotka–Volterra",
        "subtitle": "predator & prey",
        "description": "Predator-prey dynamics expressed as a reaction network; "
                       "prey and predator populations cycle out of phase.",
        "accent": "#0ea5e9",
        "builder": lotka_volterra,
        "params": {"k_birth": 1.1, "k_predation": 0.4, "k_death": 0.4,
                   "k_repro": 0.1, "interval": 0.25},
        "total_time": 40.0,
    },
    {
        "id": "repressilator",
        "title": "Repressilator",
        "subtitle": "a synthetic genetic oscillator (Hill kinetics)",
        "description": "Three genes repressing each other in a ring (Elowitz & "
                       "Leibler 2000). Hill-function rate laws — not mass-action.",
        "accent": "#ef4444",
        "builder": repressilator,
        "params": {"a": 10.0, "n": 3.0, "beta": 1.0, "interval": 0.5},
        "total_time": 50.0,
    },
]


def run_config(cfg):
    core = allocate_core()
    doc = cfg["builder"](**cfg["params"])
    sim = Composite({"state": doc}, core=core)
    t0 = time.perf_counter()
    sim.run(cfg["total_time"])
    elapsed = time.perf_counter() - t0
    if elapsed > TIMEOUT_S:
        raise TimeoutError(f"{cfg['id']} exceeded {TIMEOUT_S}s")
    rows = gather_emitter_results(sim)[("emitter",)]
    times = [r["time"] for r in rows]
    species = sorted(rows[-1]["concentrations"].keys())
    series = {s: [r["concentrations"].get(s, 0.0) for r in rows] for s in species}
    return {"doc": doc, "times": times, "series": series, "elapsed": elapsed,
            "n": len(rows)}


def make_diagram(doc, accent):
    """Render a simplified bigraph PNG; return data-URI or None on failure."""
    try:
        from bigraph_viz import plot_bigraph

        node_colors = {
            ("simbio",): accent,
            ("emitter",): "#8b5cf6",
            ("stores",): "#e0e7ff",
        }
        plot_bigraph(
            state=doc,
            out_dir=HERE,
            filename="_bigraph",
            file_format="png",
            remove_process_place_edges=True,
            rankdir="LR",
            node_fill_colors=node_colors,
            node_label_size="16pt",
            port_labels=False,
            dpi="150",
        )
        png = os.path.join(HERE, "_bigraph.png")
        with open(png, "rb") as f:
            uri = "data:image/png;base64," + base64.b64encode(f.read()).decode()
        os.remove(png)
        return uri
    except Exception as exc:  # pragma: no cover - diagram is best-effort
        print(f"  (diagram skipped: {exc})")
        return None


def plotly_div(cfg, result):
    traces = []
    for name, ys in result["series"].items():
        traces.append({
            "x": result["times"], "y": ys, "name": name, "mode": "lines",
            "line": {"width": 2.5},
        })
    layout = {
        "margin": {"t": 10, "r": 10, "b": 40, "l": 50},
        "xaxis": {"title": "time", "gridcolor": "#eef2f7"},
        "yaxis": {"title": "concentration", "gridcolor": "#eef2f7"},
        "plot_bgcolor": "white", "paper_bgcolor": "white",
        "legend": {"orientation": "h", "y": -0.2},
        "height": 360,
    }
    div_id = f"chart-{cfg['id']}"
    return (
        f'<div id="{div_id}"></div>'
        f'<script>Plotly.newPlot("{div_id}", {json.dumps(traces)}, '
        f'{json.dumps(layout)}, {{displayModeBar:false, responsive:true}});</script>'
    )


def render(report_data):
    sections = []
    nav = []
    for cfg, result, diagram in report_data:
        nav.append(f'<a href="#{cfg["id"]}">{cfg["title"]}</a>')
        final = {k: round(v, 4) for k, v in
                 {s: result["series"][s][-1] for s in result["series"]}.items()}
        metrics = "".join(
            f'<div class="metric"><div class="mv">{v}</div>'
            f'<div class="ml">{k}</div></div>' for k, v in final.items()
        )
        diagram_html = (
            f'<img src="{diagram}" alt="bigraph" style="max-width:100%;border:1px solid #e5e7eb;border-radius:8px"/>'
            if diagram else '<p class="muted">architecture diagram unavailable (graphviz)</p>'
        )
        antimony = result["doc"]["simbio"]["config"].get("antimony", "")
        antimony_html = (
            f'<pre class="antimony">{html_escape(antimony.strip())}</pre>'
            if antimony else '<p class="muted">(reaction-spec model)</p>'
        )
        doc_tree = (
            f'<pre class="doc-tree">{html_escape(json.dumps(cfg_doc_view(result["doc"]), indent=2))}</pre>'
        )
        sections.append(f"""
        <section id="{cfg['id']}" style="--accent:{cfg['accent']}">
          <h2>{cfg['title']} <span class="sub">{cfg['subtitle']}</span></h2>
          <p class="muted">{cfg['description']}</p>
          <div class="cards">
            <div class="metric"><div class="mv">{result['n']}</div><div class="ml">snapshots</div></div>
            <div class="metric"><div class="mv">{result['elapsed']*1000:.0f} ms</div><div class="ml">wall time</div></div>
            <div class="metric"><div class="mv">{cfg['params'].get('interval')}</div><div class="ml">interval</div></div>
            {metrics}
          </div>
          <div class="grid">
            <div class="panel"><h3>Time series</h3>{plotly_div(cfg, result)}</div>
            <div class="panel"><h3>Bigraph architecture</h3>{diagram_html}</div>
          </div>
          <div class="grid">
            <div class="panel"><h3>Antimony model</h3><p class="muted">loaded into a genuine simbio model via libantimony → libSBML → simbio core</p>{antimony_html}</div>
            <div class="panel"><h3>PBG composite document <span class="sub">ports &amp; wires</span></h3>{doc_tree}</div>
          </div>
        </section>""")

    return TEMPLATE.format(nav="".join(nav), sections="".join(sections))


def html_escape(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def cfg_doc_view(doc):
    """JSON-friendly view of the composite document, keeping ports and wires.

    The antimony string is replaced by a short marker so the wiring (inputs /
    outputs port -> store paths) stays the focus of the tree.
    """
    view = {}
    for k, v in doc.items():
        if isinstance(v, dict) and v.get("_type") in ("process", "step"):
            node = {
                "_type": v["_type"],
                "address": v.get("address"),
            }
            if "interval" in v:
                node["interval"] = v["interval"]
            config = dict(v.get("config", {}))
            if "antimony" in config:
                config["antimony"] = "<antimony string — shown alongside>"
            if config:
                node["config"] = config
            # The whole point: keep the port -> store-path wiring visible.
            if "inputs" in v:
                node["inputs"] = v["inputs"]
            if "outputs" in v:
                node["outputs"] = v["outputs"]
            view[k] = node
        else:
            view[k] = v
    return view


TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>viva-simbio demo report</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body {{ font-family:-apple-system,Segoe UI,Roboto,sans-serif; margin:0;
         color:#0f172a; background:#f8fafc; }}
  header {{ padding:32px 40px; background:white; border-bottom:1px solid #e5e7eb; }}
  header h1 {{ margin:0 0 4px; }}
  header p {{ margin:0; color:#64748b; }}
  nav {{ position:sticky; top:0; background:white; padding:12px 40px;
        border-bottom:1px solid #e5e7eb; z-index:10; }}
  nav a {{ margin-right:20px; color:#475569; text-decoration:none; font-weight:600; }}
  nav a:hover {{ color:#6366f1; }}
  main {{ padding:24px 40px; max-width:1100px; }}
  section {{ background:white; border:1px solid #e5e7eb; border-left:4px solid var(--accent);
            border-radius:12px; padding:24px; margin-bottom:28px; }}
  h2 {{ margin:0 0 4px; }} .sub {{ color:#94a3b8; font-weight:400; font-size:.7em; }}
  .muted {{ color:#64748b; }}
  .cards {{ display:flex; flex-wrap:wrap; gap:12px; margin:16px 0; }}
  .metric {{ background:#f1f5f9; border-radius:8px; padding:10px 16px; min-width:80px; }}
  .mv {{ font-size:1.3em; font-weight:700; color:var(--accent); }}
  .ml {{ font-size:.75em; color:#64748b; text-transform:uppercase; letter-spacing:.04em; }}
  .grid {{ display:grid; grid-template-columns:1fr 1fr; gap:20px; }}
  @media (max-width:820px) {{ .grid {{ grid-template-columns:1fr; }} }}
  .panel {{ background:white; border:1px solid #eef2f7; border-radius:10px; padding:16px; }}
  .panel h3 {{ margin:0 0 12px; font-size:.95em; }}
  .doc-tree {{ background:#0f172a; color:#e2e8f0; padding:16px; border-radius:8px;
              overflow:auto; font-size:.8em; }}
  .antimony {{ background:#1e293b; color:#bae6fd; padding:16px; border-radius:8px;
              overflow:auto; font-size:.8em; line-height:1.45; }}
  details summary {{ cursor:pointer; font-weight:600; }}
</style></head>
<body>
<header>
  <h1>viva-simbio</h1>
  <p>process-bigraph wrapper for the <b>simbio</b> Chemical Reaction Network simulator — real bridge to simbio's LSODA solver.</p>
</header>
<nav>{nav}</nav>
<main>{sections}</main>
</body></html>"""


def main():
    print("Running simbio CRN configurations...")
    report_data = []
    for cfg in CONFIGS:
        print(f"  - {cfg['id']}")
        result = run_config(cfg)
        diagram = make_diagram(result["doc"], cfg["accent"])
        report_data.append((cfg, result, diagram))

    html = render(report_data)
    with open(OUTPUT, "w") as f:
        f.write(html)
    print(f"Wrote {OUTPUT}")

    import webbrowser
    webbrowser.open("file://" + OUTPUT)


if __name__ == "__main__":
    main()
