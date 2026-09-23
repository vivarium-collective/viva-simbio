# viva-simbio

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph)
wrapper for [**simbio**](https://github.com/dyscolab/simbio), a Python package
for simulating **Chemical Reaction Networks (CRNs)** built on top of
[`poincare`](https://github.com/dyscolab/poincare).

`SimbioProcess` is a **real bridge**: every `update()` assembles a genuine
`simbio.Compartment` and integrates it with simbio's LSODA solver. No kinetics
are reimplemented — simbio builds and solves the ODEs.

## What it does

You describe a CRN one of two ways:

- **Antimony** (recommended) — a human-readable reaction string, e.g. the
  Brusselator, Lotka–Volterra, or the repressilator. It is compiled to SBML by
  **libantimony** and rebuilt as a genuine simbio model (see *Loading from
  Antimony* below).
- **Reaction spec** — a dict of `species` + mass-action `reactions`, built
  directly with simbio's core.

The process exposes the species concentrations as a `map[string,float]` store
that the surrounding bigraph owns, integrates the network forward by each
`interval`, and emits the **change** in concentration per species as an additive
delta — so a sibling process (influx, dilution, transport, another reaction
module, a controller modulating rate constants) composes naturally with it.

### Loading from Antimony

> simbio 1.1.0 ships an SBML importer (`simbio.io.sbml.loads`), but it is
> **broken against the `poincare`/`symbolite` versions it resolves on PyPI**
> (a removed `symbolite.core.as_function`, an `add_species` that builds a bare
> `Variable` instead of a `Species`, a `Species(var, stoich)` call that assumes
> the old `Reactant` signature, and a `MathMLSymbol` that fails poincare's unit
> translation). Rather than vendor a fork of it, `viva-simbio` keeps the bridge
> **fully real** without touching that importer: the model is parsed by
> **libantimony** (real), the network is extracted with **libSBML** (real), and
> the model is rebuilt with simbio's own working core (`Species`, `Parameter`,
> `RateLaw`). simbio still assembles and integrates the ODEs. Each reaction's
> kinetic law is evaluated into a genuine simbio expression, so arbitrary rate
> laws — mass-action **and** Hill / Michaelis–Menten — are supported.

## Installation

```bash
# From PyPI (once published):
pip install viva-simbio
# or with uv:
uv pip install viva-simbio

# For development (editable):
uv venv --python 3.12 .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

simbio requires Python ≥ 3.12.

> Once installed, the process registers automatically via
> `bigraph_schema.package.discover` — no manual `register_link()` calls needed.

## Quick start

```python
from process_bigraph import Composite, allocate_core, gather_emitter_results
from viva_simbio.composites.crn import brusselator

core = allocate_core()
doc = brusselator(k2=3.0, interval=0.25)   # a chemical oscillator
sim = Composite({"state": doc}, core=core)
sim.run(25.0)

results = gather_emitter_results(sim)
print(results[("emitter",)][-1])   # final concentrations + time
```

Or drive the process directly with your own Antimony model:

```python
from viva_simbio import SimbioProcess

proc = SimbioProcess(config={"antimony": """
model my_oscillator
  species X = 1, Y = 1;
  k1 = 1; k2 = 3; k3 = 1; k4 = 1;
  J1: -> X; k1;
  J2: X -> Y; k2 * X;
  J3: 2 X + Y -> 3 X; k3 * X^2 * Y;
  J4: X ->; k4 * X;
end
"""}, core=core)

# read absolute concentrations, get back per-species deltas over the interval
delta = proc.update({"concentrations": {"X": 1.0, "Y": 1.0}}, interval=0.5)

# a sibling could perturb any model parameter through the `parameters` port:
delta = proc.update({"concentrations": {"X": 1.0, "Y": 1.0},
                     "parameters": {"k2": 8.0}}, interval=0.5)
```

The reaction-spec path is also available (`config={"species": ..., "reactions": ...}`,
built via `build_crn_model`) for purely mass-action networks.

## API reference

### `SimbioProcess` (process-bigraph `Process`)

| Config | Type | Default | Meaning |
|---|---|---|---|
| `antimony` | `string` | `""` | Antimony model string (primary path) |
| `species` | `map[string,float]` | `{}` | Species → initial concentration (reaction-spec path) |
| `reactions` | `list` | `[]` | Mass-action reaction specs (reaction-spec path) |
| `volume` | `float` | `1.0` | Compartment volume (reaction-spec path) |

If `antimony` is non-empty it takes precedence; otherwise `species` + `reactions`
are used. A reaction spec is a dict:

```python
{
  "name": "bind",                  # optional, defaults to r0, r1, ...
  "reactants": ["A", "B"],          # name | ["A", 2] | {"species": "A", "stoichiometry": 2}
  "products":  ["AB"],
  "rate": 1.0,                      # mass-action rate constant
}
```

| Port | Direction | Type | Meaning |
|---|---|---|---|
| `concentrations` | input | `map[string,float]` | Current absolute concentrations (integration initial condition) |
| `parameters` | input | `map[string,float]` | Optional overrides for model parameters (rate constants, Hill coefficients, …), keyed by name |
| `concentrations` | output | `map[string,float]` | Per-species concentration **delta** over the interval (composes additively) |

### Composite generators

Discoverable via `pbg_superpowers.composite_generator.discover_generators()`,
each defined as an Antimony model and wiring both input ports (`concentrations`
and `parameters`):

- `simbio_brusselator` — the Brusselator chemical oscillator.
- `simbio_lotka_volterra` — predator–prey oscillations as a reaction network.
- `simbio_repressilator` — three-gene repressilator with **Hill** kinetics.

## Architecture

```
parameters store ──┐ (rate constants, Hill coeffs)
                    ▼
shared store "concentrations" (map[string,float], absolute)
        │  read (initial condition)            ▲  delta (accumulates)
        ▼                                       │
   SimbioProcess.update(interval) ── builds simbio model, solves LSODA
        ▲
   Antimony string ──(libantimony → libSBML → simbio core)
```

simbio concept → PBG mapping:

| simbio | viva-simbio |
|---|---|
| Antimony / SBML model | parsed by libantimony + libSBML, rebuilt with simbio core |
| `Compartment` / `Species` / `Parameter` / `RateLaw` | built from the parsed network |
| `Simulator(Model).solve(values=…, t_span=…)` | called once per `update()` |
| absolute concentration trajectory | emitted as per-step **deltas** so the bigraph composes |

## Demo

```bash
source .venv/bin/activate
python demo/demo_report.py     # writes and opens demo/report.html
```

The report runs three oscillators (Brusselator, Lotka–Volterra, repressilator),
each showing its **Antimony string**, time-series charts, metrics, an
architecture diagram, and the PBG composite document with its **ports and
wires**.

## Limitations and assumptions

- Arbitrary rate laws are supported through the Antimony path (mass-action,
  Hill, Michaelis–Menten, …) by evaluating each kinetic law into a simbio
  expression. SBML *function definitions* and assignment/rate rules are not
  translated.
- The process re-solves over `(0, interval)` each step from the store's current
  concentrations; very large intervals reduce the resolution of any in-step
  dynamics that a sibling never observes (the emitted delta is still exact for
  the closed reaction system over that interval).
- simbio's own `simbio.io.sbml` importer is **not** used (it is broken in the
  released package — see *Loading from Antimony*); this wrapper bridges the
  network into simbio's core directly.
- simbio requires Python ≥ 3.12.
