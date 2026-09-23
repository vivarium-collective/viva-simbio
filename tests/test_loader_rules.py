"""Regression tests for SBML/Antimony loader handling of rules & boundary species.

These mirror features common in real BioModels that the from-scratch
SBML->simbio bridge must handle: assignment rules (derived quantities used in
rate laws), boundary species (referenced but not consumed), and rate laws that
would reduce to 0/0 at build time if locals were inlined as float literals.
"""

import numpy as np
import pytest

from viva_simbio.antimony_loader import model_from_antimony


def _solve(model):
    from simbio import Simulator

    return Simulator(model).solve(save_at=np.linspace(0, 5, 6))


def test_assignment_rule_is_not_a_state_species():
    """An assignment rule (`Total := A + B`) is inlined, not integrated."""
    ant = """model m
      species A = 10, B = 0;
      Total := A + B;
      k = 0.5;
      r: A -> B; k * A / Total;
    end"""
    model, species, _params = model_from_antimony(ant)
    assert "Total" not in species
    assert set(species) == {"A", "B"}
    res = _solve(model)
    assert "Total" not in res.data_vars


def test_boundary_species_not_consumed():
    """A boundary species ($S) drives a rate law but is not consumed."""
    ant = """model m
      species A = 10, B = 0;
      $S = 2;
      k = 0.5;
      r: A -> B; k * A * S;
    end"""
    model, species, _ = model_from_antimony(ant)
    res = _solve(model)
    # A is consumed, B produced; the closed A<->B mass is conserved
    a = float(np.asarray(res["A"])[-1])
    b = float(np.asarray(res["B"])[-1])
    assert a < 10.0 and b > 0.0
    assert abs((a + b) - 10.0) < 1e-6


def test_division_rate_law_with_assignment_denominator():
    """`X / CT` where CT is an assignment rule must not raise 0/0 at build."""
    ant = """model m
      species C2 = 1, M = 0;
      CT := C2 + M;
      k = 1.0;
      r: C2 -> M; k * M / CT;
    end"""
    # CT = C2 + M = 1 at t0 (nonzero), so M/CT is well defined — and crucially
    # building the RateLaw must not evaluate a literal 0.0/0.0.
    model, species, _ = model_from_antimony(ant)
    res = _solve(model)
    assert not np.any(np.isnan(np.asarray(res["M"])))

# A non-unit compartment with an initialAmount species + an extensive kinetic
# law (the SBML rate explicitly multiplies by the compartment). COPASI/Tellurium
# track concentration = amount/volume and the comp factor cancels in dConc/dt.
# Regression guard for the amount<->concentration volume-scaling bug that made
# BIOMD1/2/9 diverge (nRMSE ~0.4): the loader must report A's initial as
# amount/volume and decay it at rate k (volume-independent), not the raw amount
# decaying at rate k*volume.
_VOL = 1.0e-3
_AMOUNT0 = 5.0
_K = 0.7
_SBML_AMOUNT_NONUNIT_VOL = f"""<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model id="decay_amount">
    <listOfCompartments>
      <compartment id="c" size="{_VOL}"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="A" compartment="c" initialAmount="{_AMOUNT0}" hasOnlySubstanceUnits="false"/>
      <species id="B" compartment="c" initialAmount="0" hasOnlySubstanceUnits="false"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="{_K}"/>
    </listOfParameters>
    <listOfReactions>
      <reaction id="r" reversible="false">
        <listOfReactants><speciesReference species="A"/></listOfReactants>
        <listOfProducts><speciesReference species="B"/></listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>c</ci><apply><times/><ci>k</ci><ci>A</ci></apply></apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"""


def test_initial_amount_converted_to_concentration_by_volume():
    """A's reported initial value is amount/volume, not the raw amount."""
    from viva_simbio.antimony_loader import model_from_sbml
    from simbio import Simulator

    model, species, _ = model_from_sbml(_SBML_AMOUNT_NONUNIT_VOL)
    assert set(species) == {"A", "B"}
    res = Simulator(model).solve(save_at=np.linspace(0, 2.0, 21))
    a0 = float(np.asarray(res["A"])[0])
    assert a0 == pytest.approx(_AMOUNT0 / _VOL, rel=1e-6)  # 5/0.001 = 5000


def test_extensive_rate_law_decays_independent_of_volume():
    """dConc/dt = -k*[A]: the comp factor cancels, so A decays at rate k
    regardless of compartment volume (it did NOT before the fix)."""
    from viva_simbio.antimony_loader import model_from_sbml
    from simbio import Simulator

    model, _species, _ = model_from_sbml(_SBML_AMOUNT_NONUNIT_VOL)
    t = np.linspace(0, 2.0, 41)
    res = Simulator(model).solve(save_at=t)
    a = np.asarray(res["A"])
    expected = (_AMOUNT0 / _VOL) * np.exp(-_K * t)
    # nRMSE against the analytic concentration trajectory.
    rng = expected.max() - expected.min()
    nrmse = float(np.sqrt(np.mean((a - expected) ** 2)) / rng)
    assert nrmse < 1e-3
