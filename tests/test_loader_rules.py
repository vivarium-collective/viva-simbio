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


def _run(model, t):
    import simbio

    return simbio.Simulator(model).solve(save_at=t)


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


# --- Rate rules (explicit ODEs dX/dt = f) --------------------------------
#
# Many BioModels express their dynamics via `rateRule`s on species declared
# boundaryCondition=true (and sometimes on parameters) rather than reactions
# (e.g. BIOMD0000000006 / Goldbeter 1991: u' = ..., v' = ...). Before the fix
# the loader ignored rateRules entirely, so those targets were left as fixed
# constants and the model emitted no dynamics (empty/flat output). These tests
# assert a rateRule target is INTEGRATED — it varies over time and is finite.


def test_rate_rule_species_is_dynamic():
    """A boundary species governed only by a rate rule must be integrated,
    not left as a flat constant. `X' = a - X` relaxes from 0 toward a=2."""
    ant = """model m
      species $X = 0;
      a = 2;
      X' = a - X;
    end"""
    model, species, _params = model_from_antimony(ant)
    assert species == ["X"]
    import simbio

    t = np.linspace(0, 5, 11)
    res = simbio.Simulator(model).solve(save_at=t)
    x = np.asarray(res["X"])
    assert not np.any(np.isnan(x))  # finite, integrator did not fail
    assert x.max() - x.min() > 1e-3  # actually varies, not a flat constant
    # Analytic: X(t) = a*(1 - e^{-t}); check the approach to the a=2 asymptote.
    assert x[-1] == pytest.approx(2.0 * (1 - np.exp(-5.0)), rel=1e-2)


def test_rate_rule_species_not_consumed_by_reaction():
    """A rate-rule species is driven ONLY by its rule; a reaction referencing
    it does not additionally change it (SBML rate-rule semantics)."""
    # X has a rate rule (X' = 1, so X == t), and also appears as a reactant in
    # a reaction producing P. X must follow its rule (reach ~5 at t=5), and P
    # must still be produced from X's trajectory.
    ant = """model m
      species $X = 0, P = 0;
      X' = 1;
      kf = 0.3;
      r: X -> P; kf * X;
    end"""
    model, _species, _ = model_from_antimony(ant)
    import simbio

    res = simbio.Simulator(model).solve(save_at=np.linspace(0, 5, 51))
    x = np.asarray(res["X"])
    p = np.asarray(res["P"])
    assert x[-1] == pytest.approx(5.0, rel=1e-3)  # X == t, unaffected by reaction
    assert p[-1] > 0.0  # P still produced from X


def test_rate_rule_parameter_becomes_integrated_state():
    """A parameter that is a rate-rule target must become an integrable state
    variable (via `initial`), not a fixed constant. `p' = 2` from p0=1."""
    ant = """model m
      species $X = 0;
      p = 1;
      p' = 2;
      X' = p;
    end"""
    model, species, params = model_from_antimony(ant)
    assert "p" in params
    import simbio

    t = np.linspace(0, 3, 4)
    res = simbio.Simulator(model).solve(save_at=t)
    p = np.asarray(res["p"])
    x = np.asarray(res["X"])
    assert not np.any(np.isnan(p))
    # p(t) = 1 + 2t -> [1, 3, 5, 7]
    assert np.allclose(p, 1.0 + 2.0 * t, rtol=1e-4)
    # X(t) = integral of p = t + t^2 -> driven by the evolving parameter
    assert np.allclose(x, t + t**2, atol=1e-3)


def test_rate_rule_coupled_with_assignment_rule():
    """Rate rules whose RHS references an assignment-rule variable must resolve
    (assignment rules are inlined before rate-rule RHS is built). Mirrors the
    BIOMD6 shape: y := a + b (assignment), x' = y (rate)."""
    ant = """model m
      species $x = 0;
      a = 1; b = 2;
      y := a + b;
      x' = y;
    end"""
    model, species, _ = model_from_antimony(ant)
    assert "y" not in species  # assignment var, inlined not integrated
    import simbio

    t = np.linspace(0, 4, 9)
    res = simbio.Simulator(model).solve(save_at=t)
    x = np.asarray(res["x"])
    assert not np.any(np.isnan(x))
    # y == 3 constant, so x' = 3 -> x == 3t
    assert np.allclose(x, 3.0 * t, atol=1e-3)


# --- Math builtins, functionDefinitions, piecewise/logical -----------------
#
# libSBML L3 rate laws reference math functions (exp, ln, sin, sqrt, ...),
# SBML functionDefinitions (reusable lambdas), the SBML `piecewise`, and
# logical/relational operators. Before the fix these were absent from the eval
# scope, so hundreds of BioModels failed with `NameError: name 'exp' is not
# defined` / `name 'function_1' is not defined`, or `SyntaxError` on `&&`/
# piecewise. These now evaluate symbolically and integrate.


def test_math_builtin_exp_in_rate_law():
    """`exp(...)` in a kinetic law evaluates symbolically and integrates."""
    ant = """model m
      species X = 1;
      k = 0.5;
      X' = -k * exp(-X);
    end"""
    model, species, _ = model_from_antimony(ant)
    assert species == ["X"]
    x = np.asarray(_run(model, np.linspace(0, 3, 7))["X"])
    assert np.all(np.isfinite(x))
    assert x[-1] < x[0]  # decays (dX/dt = -k e^{-X} < 0), actual dynamics


def test_math_builtin_log_and_trig_supported():
    """ln / sqrt / sin build without NameError and integrate to finite values."""
    ant = """model m
      species X = 2;
      X' = -sqrt(X) - ln(X) + sin(X) * 0;
    end"""
    model, _species, _ = model_from_antimony(ant)
    x = np.asarray(_run(model, np.linspace(0, 1, 6))["X"])
    assert np.all(np.isfinite(x))
    assert x.max() - x.min() > 1e-6


def test_function_definition_is_substituted():
    """An SBML functionDefinition (a reusable lambda) is applied in a rate law.
    Mirrors BIOMD0000000089 (NameError: name 'function_1' before the fix)."""
    ant = """function mm(s, vmax, km)
      vmax * s / (km + s)
    end
    model m
      species S = 10, P = 0;
      vm = 1; k = 2;
      r: S -> P; mm(S, vm, k);
    end"""
    model, species, _ = model_from_antimony(ant)
    assert set(species) == {"S", "P"}
    res = _run(model, np.linspace(0, 5, 6))
    s = np.asarray(res["S"])
    p = np.asarray(res["P"])
    assert s[-1] < s[0] and p[-1] > 0.0  # MM consumption, dynamics present
    assert abs((s[-1] + p[-1]) - 10.0) < 1e-6  # mass conserved


def test_function_definition_referencing_math_builtin():
    """A functionDefinition body may itself call a math builtin (exp)."""
    ant = """function g(x)
      exp(-x)
    end
    model m
      species S = 2, P = 0;
      k = 1;
      r: S -> P; k * g(S) * S;
    end"""
    model, _species, _ = model_from_antimony(ant)
    res = _run(model, np.linspace(0, 3, 4))
    s = np.asarray(res["S"])
    assert np.all(np.isfinite(s)) and s[-1] < s[0]


def test_piecewise_selects_branch():
    """`piecewise(1, X < 5, 0)`: X rises at rate 1 while X<5, then stops.
    Mirrors the piecewise construct in BIOMD0000000162 (SyntaxError before)."""
    ant = """model m
      species X = 0;
      X' = piecewise(1, X < 5, 0);
    end"""
    model, _species, _ = model_from_antimony(ant)
    x = np.asarray(_run(model, np.linspace(0, 10, 11))["X"])
    assert np.all(np.isfinite(x))
    assert 4.0 < x[-1] < 6.0  # capped near the X==5 threshold


def test_piecewise_logical_and_builds_and_is_correct_in_window():
    """`&&` inside a piecewise condition no longer SyntaxErrors; the AND is
    rewritten to a product of Real indicators. Started inside the (1,4) window,
    the rate is the in-window value (2), so Y rises from 2 toward 3."""
    ant = """model m
      species Y = 2;
      Y' = piecewise(2, (Y > 1) && (Y < 4), 0.5);
    end"""
    model, _species, _ = model_from_antimony(ant)
    y = np.asarray(_run(model, np.linspace(0, 0.5, 2))["Y"])
    assert np.all(np.isfinite(y))
    assert y[-1] == pytest.approx(3.0, abs=1e-2)  # rate 2 over dt=0.5


def test_min_max_supported():
    """SBML min/max (built via the abs identity) evaluate and integrate.
    `X' = max(0, 1 - X)` relaxes X from 0 toward 1 (= 1 - e^{-t} while X<1)."""
    ant = """model m
      species X = 0;
      X' = max(0, 1 - X);
    end"""
    model, _species, _ = model_from_antimony(ant)
    t = np.linspace(0, 5, 11)
    x = np.asarray(_run(model, t)["X"])
    assert np.all(np.isfinite(x))
    assert x.max() - x.min() > 1e-3  # varies
    assert x[-1] == pytest.approx(1 - np.exp(-5.0), rel=1e-2)  # approaches 1


def test_unsupported_construct_fails_cleanly():
    """A genuinely unparseable rate-law fragment raises a clear ValueError
    naming the construct, not a bare SyntaxError leaking from eval."""
    from viva_simbio.antimony_loader import _eval_rate_law

    with pytest.raises(ValueError, match="Unsupported construct"):
        _eval_rate_law("1 <", {})
