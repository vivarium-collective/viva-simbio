"""Tests for the simbio (Antimony) composite generators and registration."""

import numpy as np
from process_bigraph import Composite, gather_emitter_results

from viva_simbio import build_core
from viva_simbio.composites.crn import brusselator, lotka_volterra, repressilator

GENERATORS = ["simbio_brusselator", "simbio_lotka_volterra", "simbio_repressilator"]


def _run(doc, total_time):
    # build_core() registers SimbioProcess so the composite's `local:SimbioProcess`
    # address resolves (a bare allocate_core() does not see the editable workspace pkg).
    core = build_core()
    sim = Composite({"state": doc}, core=core)
    sim.run(total_time)
    return gather_emitter_results(sim)[("emitter",)]


def test_generators_registered():
    from viva_superpowers.composite_generator import _REGISTRY

    for gen in GENERATORS:
        matches = [eid for eid in _REGISTRY if eid.endswith("." + gen)]
        assert matches, f"{gen} missing; have {list(_REGISTRY)[:5]}"


def test_document_wires_both_inputs():
    """Both process input ports must be wired (concentrations + parameters)."""
    doc = brusselator()
    inputs = doc["simbio"]["inputs"]
    assert inputs["concentrations"] == ["stores", "concentrations"]
    assert inputs["parameters"] == ["stores", "parameters"]
    # the parameters store is populated with the model's rate constants
    assert set(doc["stores"]["parameters"]) == {"k1", "k2", "k3", "k4"}


def test_brusselator_oscillates():
    rows = _run(brusselator(interval=0.25), total_time=20.0)
    xs = [r["concentrations"]["X"] for r in rows]
    # an oscillator swings across a wide range rather than settling
    assert max(xs) - min(xs) > 1.0


def test_lotka_volterra_runs():
    rows = _run(lotka_volterra(interval=0.25), total_time=20.0)
    prey = [r["concentrations"]["prey"] for r in rows]
    assert max(prey) - min(prey) > 1.0
    assert min(prey) >= 0.0  # populations stay non-negative


def test_repressilator_hill_kinetics_runs():
    rows = _run(repressilator(interval=0.5), total_time=40.0)
    p1 = [r["concentrations"]["p1"] for r in rows]
    # protein 1 is produced and oscillates under Hill repression
    assert max(p1) > 1.0


def test_parameter_store_perturbs_dynamics():
    """Changing the parameters store changes the trajectory (input wire is live)."""
    base = brusselator(k2=3.0, interval=0.25)
    hot = brusselator(k2=8.0, interval=0.25)
    b_last = _run(base, 6.0)[-1]["concentrations"]
    h_last = _run(hot, 6.0)[-1]["concentrations"]
    assert b_last != h_last
