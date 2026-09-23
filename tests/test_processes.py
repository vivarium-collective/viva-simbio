"""Tests for SimbioProcess — the real simbio CRN bridge."""

import numpy as np
import pytest
from process_bigraph import Composite, allocate_core, gather_emitter_results

from viva_simbio import SimbioProcess, build_core, build_crn_model

REVERSIBLE = {
    "species": {"A": 1.0, "B": 2.0, "AB": 0.0},
    "reactions": [
        {"name": "bind", "reactants": ["A", "B"], "products": ["AB"], "rate": 1.0},
        {"name": "unbind", "reactants": ["AB"], "products": ["A", "B"], "rate": 0.3},
    ],
    "volume": 1.0,
}


def test_build_crn_model():
    model, reaction_names = build_crn_model(
        REVERSIBLE["species"], REVERSIBLE["reactions"], REVERSIBLE["volume"]
    )
    assert reaction_names == ["bind", "unbind"]
    # species are real attributes on the dynamically built Compartment
    for name in ("A", "B", "AB"):
        assert hasattr(model, name)


def test_ports_are_dicts():
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    assert isinstance(proc.inputs(), dict)
    assert isinstance(proc.outputs(), dict)
    assert set(proc.inputs()) == {"concentrations", "parameters"}
    assert set(proc.outputs()) == {"concentrations"}


def test_initial_state():
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    assert proc.initial_state()["concentrations"] == {"A": 1.0, "B": 2.0, "AB": 0.0}


def test_update_returns_deltas():
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    out = proc.update({"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}, interval=1.0)
    d = out["concentrations"]
    assert set(d) == {"A", "B", "AB"}
    # A + B -> AB consumes A and B, produces AB
    assert d["AB"] > 0
    assert d["A"] < 0 and d["B"] < 0


def test_mass_conservation():
    """A and B are conserved in (free + bound) form by the closed CRN."""
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    out = proc.update({"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}, interval=2.0)
    d = out["concentrations"]
    # d[AB] of AB formed must equal A consumed and B consumed
    assert abs(d["AB"] + d["A"]) < 1e-6
    assert abs(d["AB"] + d["B"]) < 1e-6


def test_rate_override_input():
    """A higher 'bind' rate fed via the parameters port produces more product."""
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    state = {"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}
    slow = proc.update({**state, "parameters": {"bind": 0.5}}, interval=1.0)
    fast = proc.update({**state, "parameters": {"bind": 5.0}}, interval=1.0)
    assert fast["concentrations"]["AB"] > slow["concentrations"]["AB"]


# --- Antimony path -----------------------------------------------------------

BRUSSELATOR_ANT = """\
model brusselator
  species X = 1, Y = 1;
  k1 = 1; k2 = 3; k3 = 1; k4 = 1;
  J1: -> X; k1;
  J2: X -> Y; k2 * X;
  J3: 2 X + Y -> 3 X; k3 * X^2 * Y;
  J4: X ->; k4 * X;
end
"""


def test_antimony_model_loads_and_steps():
    proc = SimbioProcess(config={"antimony": BRUSSELATOR_ANT}, core=allocate_core())
    # initial_state() triggers the lazy build (libantimony -> libSBML -> simbio)
    assert proc.initial_state()["concentrations"] == {"X": 1.0, "Y": 1.0}
    assert proc._species_names == ["X", "Y"]
    assert proc._param_names == ["k1", "k2", "k3", "k4"]
    out = proc.update({"concentrations": {"X": 1.0, "Y": 1.0}}, interval=0.5)
    assert set(out["concentrations"]) == {"X", "Y"}


def test_antimony_parameter_override_changes_dynamics():
    proc = SimbioProcess(config={"antimony": BRUSSELATOR_ANT}, core=allocate_core())
    state = {"concentrations": {"X": 1.0, "Y": 1.0}}
    base = proc.update(state, interval=0.5)["concentrations"]
    hot = proc.update({**state, "parameters": {"k2": 8.0}}, interval=0.5)["concentrations"]
    assert base != hot


def test_antimony_hill_kinetics_supported():
    """Non-mass-action (Hill) rate laws load through libSBML -> simbio."""
    hill = """\
model hill
  species S = 5, P = 0;
  vmax = 2; K = 1; n = 4;
  r: S -> P; vmax * S^n / (K^n + S^n);
end
"""
    proc = SimbioProcess(config={"antimony": hill}, core=allocate_core())
    out = proc.update({"concentrations": {"S": 5.0, "P": 0.0}}, interval=1.0)
    assert out["concentrations"]["P"] > 0
    assert out["concentrations"]["S"] < 0


def test_composite_run_accumulates_in_store():
    """Deltas accumulate in the shared store across a real Composite run."""
    # build_core() registers SimbioProcess so `local:SimbioProcess` resolves.
    core = build_core()
    doc = {
        "simbio": {
            "_type": "process",
            "address": "local:SimbioProcess",
            "config": REVERSIBLE,
            "interval": 0.5,
            "inputs": {"concentrations": ["stores", "concentrations"]},
            "outputs": {"concentrations": ["stores", "concentrations"]},
        },
        "stores": {"concentrations": dict(REVERSIBLE["species"])},
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {"concentrations": "map[string,float]", "time": "float"}},
            "inputs": {
                "concentrations": ["stores", "concentrations"],
                "time": ["global_time"],
            },
        },
    }
    sim = Composite({"state": doc}, core=core)
    sim.run(5.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert len(rows) > 1
    final = rows[-1]["concentrations"]
    # AB grows toward steady state; A+AB and B+AB stay conserved
    assert final["AB"] > 0.5
    assert abs((final["A"] + final["AB"]) - 1.0) < 1e-6
    assert abs((final["B"] + final["AB"]) - 2.0) < 1e-6


def test_stoichiometry():
    """2A -> A2 consumes A twice as fast as it makes A2."""
    cfg = {
        "species": {"A": 2.0, "A2": 0.0},
        "reactions": [
            {"name": "dimerize", "reactants": [["A", 2]], "products": ["A2"], "rate": 1.0}
        ],
        "volume": 1.0,
    }
    proc = SimbioProcess(config=cfg, core=allocate_core())
    out = proc.update({"concentrations": {"A": 2.0, "A2": 0.0}}, interval=1.0)["concentrations"]
    assert abs(out["A"] + 2 * out["A2"]) < 1e-6  # 2 A consumed per A2 formed
