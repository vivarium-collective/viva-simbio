"""Tests for the canonical model-source Steps/Process (parity with copasi/tellurium)."""

import textwrap

import pytest
from process_bigraph import allocate_core

from viva_simbio import (
    SimbioSteadyStateStep,
    SimbioUTCProcess,
    SimbioUTCStep,
    register_simbio_types,
)

# A small Antimony model written to a temp .ant file exercises the file-path
# branch of load_model_source without needing a bundled SBML fixture.
MODEL = textwrap.dedent("""\
    model decay
      species A = 10, B = 0;
      k = 0.5;
      r: A -> B; k * A;
    end
""")


@pytest.fixture
def core():
    return register_simbio_types(allocate_core())


@pytest.fixture
def model_file(tmp_path):
    path = tmp_path / "decay.ant"
    path.write_text(MODEL)
    return str(path)


def test_utc_step_returns_numeric_result(core, model_file):
    step = SimbioUTCStep(
        config={"model_source": model_file, "time": 10.0, "n_points": 11}, core=core
    )
    out = step.update({})
    assert "result" in out
    r = out["result"]
    assert set(r.keys()) >= {"time", "columns", "values"}
    assert len(r["time"]) == 11
    assert len(r["values"]) == 11
    assert all(len(row) == len(r["columns"]) for row in r["values"])
    # A decays into B
    a_idx = r["columns"].index("A")
    assert r["values"][-1][a_idx] < r["values"][0][a_idx]


def test_utc_step_accepts_string_source(core):
    step = SimbioUTCStep(config={"model_source": MODEL, "time": 5.0, "n_points": 6}, core=core)
    r = step.update({})["result"]
    assert len(r["time"]) == 6


def test_utc_step_rejects_n_points_below_two(core, model_file):
    step = SimbioUTCStep(config={"model_source": model_file, "time": 10.0, "n_points": 1}, core=core)
    with pytest.raises(ValueError, match="n_points must be >= 2"):
        step.update({})


def test_utc_step_species_override_couples(core, model_file):
    step = SimbioUTCStep(config={"model_source": model_file, "time": 5.0, "n_points": 2}, core=core)
    base = step.update({})["result"]
    perturbed = step.update({"species_counts": {"A": 100.0}})["result"]
    a = base["columns"].index("A")
    assert perturbed["values"][0][a] == 100.0
    assert perturbed["values"][-1][a] > base["values"][-1][a]


def test_steady_state_step(core, model_file):
    step = SimbioSteadyStateStep(config={"model_source": model_file}, core=core)
    out = step.update({})
    ss = out["steady_state_concentrations"]
    # Irreversible decay: A -> 0, B -> 10 at steady state
    assert ss["A"] < 1e-3
    assert abs(ss["B"] - 10.0) < 1e-2


def test_utc_process_steps_in_time(core, model_file):
    proc = SimbioUTCProcess(config={"model_source": model_file}, core=core)
    init = proc.initial_state()
    assert init["time"] == 0.0
    assert init["species_concentrations"]["A"] == 10.0
    out = proc.update({}, 5.0)
    assert out["time"] == 5.0
    assert out["species_concentrations"]["A"] < 10.0


def test_missing_model_source_raises(core):
    with pytest.raises(ValueError, match="model_source"):
        SimbioUTCStep(config={"model_source": ""}, core=core).update({})
