# Next steps for simbio

This workspace IS the model. You build it incrementally through phases.
Most work happens in the **dashboard** — a side-rail of tabs, each with
the buttons you need. Skills (Claude Code) are the alternative for
code-writing tasks. The set of tabs evolves as the dashboard grows;
see what's live at the URL printed by `scripts/serve.sh`.

```
bash scripts/serve.sh   # opens browser at http://localhost:<port>
```

## 0 — One-time setup

- [ ] Create the venv: `uv venv .venv`
- [ ] Activate it: `source .venv/bin/activate`
- [ ] Install workspace deps: `uv pip install -e ".[dev]"`
- [ ] Lint: `python3 scripts/lint-workspace.py` should print `workspace lint: OK`
- [ ] Commit + (eventually) push: `git init && git add -A && git commit -m "feat: workspace bootstrap"`

> **If `uv pip install` fails with "vivarium-workbench was not found in
> the package registry":** vivarium-workbench isn't on PyPI yet. The
> template's init script auto-pins it when a sibling `../vivarium-dashboard/`
> checkout exists, but if you scaffolded elsewhere you'll need to add it
> manually to `pyproject.toml`:
>
> ```toml
> [tool.uv.sources]
> vivarium-workbench = { path = "/path/to/vivarium-dashboard", editable = true }
> ```
>
> Re-run `uv pip install -e ".[dev]"` after the edit.

## The dashboard tabs

> The dashboard's tab set evolves as the platform grows. The list below
> reflects the rail at the time this template was rendered; the source
> of truth is the live UI you opened at `http://localhost:<port>/`.

| Rail label | Hash route | What it's for |
|---|---|---|
| **Workspace inputs** | `#workspace-inputs` | External resources — datasets, references (PDFs auto-extract metadata), expert docs. Simulation modules live in the Registry tab. |
| **Registry** | `#registry` | Browse the curated catalog (`scripts/_catalog/modules.json`), install pbg-* modules, and view live `build_core()` introspection of discovered Processes/Steps/Types. |
| **Composites** | `#simulation-setup` | Browse and explore composites available in the workspace. Each composite packages a runnable simulation state with parameters you can configure. |
| **Investigations** | `#investigations` | Declarative research recipes — pick a composite, declare simulations (single / sweep / seeds), name observables, choose visualizations. Run, save, compare. |
| **Visualizations** | `#visualizations` | Charts rendered from observable trajectories. Configure a registered Visualization class with specific settings, or describe a chart in natural language and let `/pbg-viz` generate it. |
| **GitHub Branches** | `#branches` | All `stage/*` branches in the workspace with one-click merge / PR / diff actions. Branches accumulate as the dashboard creates them on actions (Add observable, Install module, …); use this tab to land them on main. |

(Note: the **Composites** tab still routes under `#simulation-setup` for
backwards compatibility — the rail label was renamed but the URL fragment
wasn't. Use the table above when bookmarking focused panels.)

## 1 — Workspace inputs (any time)

Inputs aren't a sequential stage — load them whenever they're useful.

- **Datasets** — experimental data the model validates against. Drag-drop a file.
- **References** — paper PDFs. Drop the file; pypdf extracts metadata; bibtex auto-generates.
- **Expert docs** — lab notes, curated reviews, working drafts. Drag-drop a PDF.

Simulation modules (Python packages that contribute Processes/Steps) are
installed from the **Registry** tab, not here.

Each `+ Add` lands on a `stage/*` branch you can merge from the GitHub Branches tab.

## 2 — Registry

Two coordinated views on one tab:

- **Available modules** — the curated `pbg-*` catalog in `scripts/_catalog/modules.json`, regenerated from the `vivarium-collective` GitHub org by `scripts/_catalog/sync-catalog.py`. Click **Install** to `git submodule add` into `external/<name>/`, `pip install -e` it into the workspace venv, and append it to `pyproject.toml` deps + `workspace.yaml.imports`.
- **Discovered processes / types** — live introspection of `build_core()`. Every Process/Step/Emitter/Visualization/Type that `allocate_core()` finds via [bigraph-schema's discovery convention][discovery]. Click **Refresh registry** after an Install to see new entries.

Empty Discovered view usually means the venv isn't built (`uv pip install -e ".[dev]"`) or the import isn't installed.

## 3 — Composites

Composite documents — process-bigraph state trees — live here. Pick one, inspect its wiring (process tree, store paths, emitters), and run it from the **Composite Explorer** sub-page to confirm it produces sensible output before wiring it into an investigation. Composites you write directly into `viva_simbio/composites/` show up automatically.

## 4 — Studies and Investigations

A **study** (`studies/<slug>/study.yaml`) is one research question, with a
`phase:` of `Design | Build | Simulate | Evaluate | Decide` and an 8-section
canonical structure (`purpose`, `simulation_set`, `model_change`,
`key_assumptions`, `readouts`, `behavior_tests`, `conclusion_logic`,
`limitations`). Use `/pbg-study <slug>` from the [pbg-superpowers] plugin to
scaffold one. Decide-phase studies can record `followup_proposals[]`; seed a
child study from any proposal with `/pbg-study seed-from-followup
<parent>/<proposal_id>` (the child gets a `seeded_from:` provenance stamp).

An **investigation** (`investigations/<slug>/investigation.yaml`) groups
related studies and orders them as a DAG via `pipeline_gate.prerequisites`.
Use `/pbg-investigation <slug>` to scaffold.

The legacy per-investigation "run a composite + observables + visualizations"
flow still lives in the dashboard's Investigations tab — pick a composite,
declare simulations (single / sweep / seeds), name observables, run, save,
compare.

## 5 — Visualizations

Charts rendered from observable trajectories. Two creation paths:

- **Configure a registered class** — pick from the Visualization classes discovered in the Registry (subclasses of `pbg_superpowers.visualization.Visualization` in any installed pbg-* package), give it settings.
- **Generate from natural language** — describe what you want; the dashboard writes a request file and prompts you to run `/pbg-viz <name>`, which scaffolds a new Visualization function into `viva_simbio/visualizations/` and commits it on a stage branch.

## 6 — GitHub Branches

Every dashboard action (install a module, add an observable, generate a viz, …) creates a `stage/*` branch and commits the change there. This tab is the landing UI:

- **Copy gh** — `gh pr create` command for that branch.
- **Copy git merge** — local merge command, no PR.
- **Show diff** — preview without leaving the tab.

When all your in-flight stage branches are merged, the rail badge clears.

### Container image (GHCR)

Each push triggers `.github/workflows/build-and-push.yml`, which builds your workspace `Dockerfile` and publishes `ghcr.io/<owner>/<repo>:<tag>` (tags: `<branch>`, `pr-N`, `sha-<short>`, `latest` on default). The workflow also **auto-syncs the GHCR package's visibility to the source repo's visibility** on each push — so if your repo is public, `docker pull ghcr.io/<owner>/<repo>:<tag>` works anonymously without anyone running a manual "Change visibility" toggle. If you intentionally want the image private even though the repo is public, either set the repo to private or remove the `Sync GHCR package visibility` step from the workflow.

## Two paths to the same workspace

- **Browser dashboard** is the primary UI. Pure-CLI workflows are supported via `python3 scripts/lint-workspace.py` + direct YAML edits, but the buttons handle branch-creation and PR setup for you.
- **Claude Code + [pbg-superpowers]** adds skills that pair with the dashboard: `/pbg-study` and `/pbg-investigation` (scaffold a study or investigation), `/pbg-viz` (generate a visualization), `/pbg-expert` (wrap a simulator; pass `--lightweight` for in-workspace, or omit for a sibling repo + composite form via `<name> <tools…>`), `/pbg-report` (regenerate per-model reports). The dashboard and the skills share the same `workspace.yaml`.

## Config files to know

- `workspace.yaml` — canonical state (observables, simulations, visualizations, imports, datasets, …)
- `studies/<slug>/study.yaml` — one research question (8 canonical sections, `phase:` lifecycle)
- `investigations/<slug>/investigation.yaml` — DAG of studies via `pipeline_gate.prerequisites`
- `references/{papers.bib, claims.yaml}` — bibliography + claim → paper mapping
- `.pbg/schemas/{workspace,study,investigation}.schema.json` — JSON schemas (lint enforces them)

## Focused dashboard panels

Open just one panel for a targeted interaction:

| URL | What |
|---|---|
| `http://localhost:<port>/?focus=workspace-inputs` | Datasets / references / expert docs |
| `http://localhost:<port>/?focus=registry` | Module catalog + installed modules |
| `http://localhost:<port>/?focus=simulation-setup` | Composites browser (legacy route name) |
| `http://localhost:<port>/?focus=investigations` | Investigation specs + runs |
| `http://localhost:<port>/?focus=visualizations` | Visualization lifecycle |
| `http://localhost:<port>/?focus=branches` | GitHub Branches lander |
| `http://localhost:<port>/?focus=composite-explore&id=<composite>` | A single composite's explorer |

Skills in Claude Code can `open <url>` to surface a focused interaction without dropping users into the full rail.

## When stuck

1. `python3 scripts/lint-workspace.py` — most workspace-shape issues surface here.
2. Look at the **Next step** banner at the top of the dashboard — it points at the right tab.
3. File an issue at <https://github.com/vivarium-collective/pbg-superpowers/issues>.

[discovery]: https://github.com/vivarium-collective/pbg-superpowers/blob/main/docs/conventions/discovery.md
[pbg-superpowers]: https://github.com/vivarium-collective/pbg-superpowers
