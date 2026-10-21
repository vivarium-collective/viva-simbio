"""Load a genuine simbio model from an Antimony string.

Why not ``simbio.io.sbml.loads``? simbio 1.1.0's bundled SBML/MathML importer
is broken against the ``poincare``/``symbolite`` versions it resolves on PyPI
(``symbolite.core.as_function`` was removed; ``add_species`` builds a bare
``Variable`` instead of a ``Species``; ``get_species_reference`` calls the
keyword-only ``Species`` positionally as if it were a ``Reactant``; and the
``MathMLSymbol`` rate-law type fails poincare's unit translation). Repairing
all of that means maintaining a fork of simbio's importer.

Instead this module keeps the bridge **fully real** without touching that
importer: the model text is parsed by **libantimony** (real), the reaction
network is extracted with **libSBML** (real), and the model is rebuilt with
simbio's own working core constructors (`Species`, `Parameter`, `RateLaw`) —
so simbio still assembles and integrates the ODEs. Each reaction's kinetic law
is evaluated into a genuine simbio/symbolite expression, so arbitrary rate laws
(mass-action, Hill, Michaelis–Menten, ...) are supported, not just mass-action.
"""

from __future__ import annotations

import ast
import math
import re
from typing import Any


def antimony_to_sbml(antimony: str, *, model_name: str | None = None) -> tuple[str, str]:
    """Compile an Antimony string to SBML via libantimony. Returns (sbml, name)."""
    import antimony as sb

    sb.clearPreviousLoads()
    code = sb.loadAntimonyString(antimony)
    if code < 0:
        raise ValueError(f"Antimony parse error: {sb.getLastError()}")
    name = model_name or sb.getMainModuleName()
    return sb.getSBMLString(name), name


# --- Math builtins & logical handling for rate laws ----------------------
#
# libSBML L3 formulas reference math functions (exp, ln, sin, sqrt, ...), the
# SBML `piecewise`, and logical/relational operators that are NOT ordinary
# Python-evaluable names/operators over symbolic Species/Parameter objects.
# Without these in the eval scope the loader raised `NameError: name 'exp' is
# not defined` (~hundreds of BioModels) or `SyntaxError` on `&&`/`||`/piecewise.
# The functions below operate SYMBOLICALLY on simbio/symbolite expressions, so
# they translate through poincare to the solver (they do not just handle floats).

_MATH_SCOPE_CACHE: dict[str, Any] | None = None


def _piecewise(*args: Any) -> Any:
    """SBML ``piecewise(v1, c1, v2, c2, ..., [default])``.

    Conditions arrive as Real 0/1 indicators (AST-rewritten below), so the
    selection is built arithmetically — ``c*v + (1-c)*rest`` — which simbio's
    backend can translate. Earlier conditions take precedence (first true wins).
    """
    has_default = len(args) % 2 == 1
    default = args[-1] if has_default else 0.0
    n_pairs = (len(args) - (1 if has_default else 0)) // 2
    pairs = [(args[2 * i], args[2 * i + 1]) for i in range(n_pairs)]
    result: Any = default
    for value, cond in reversed(pairs):
        result = cond * value + (1.0 - cond) * result
    return result


def _pw_indicator(cond: Any) -> Any:
    """Coerce a (symbolic) Boolean comparison to a Real 0/1 indicator.

    simbio's backend cannot translate symbolic boolean AND/OR, but it can
    translate a comparison multiplied into Real arithmetic. An all-literal
    comparison folds to a Python ``bool`` during eval; map it to 0.0/1.0.
    """
    if isinstance(cond, bool):
        return 1.0 if cond else 0.0
    from symbolite.abstract import real as sr

    return sr.mul(cond, 1.0)


def _pw_or(a: Any, b: Any) -> Any:
    """Logical OR over Real 0/1 indicators: a + b - a*b (stays in {0,1})."""
    return a + b - a * b


def _math_scope() -> dict[str, Any]:
    """Build (once) the symbolic math-function scope shared by every rate law."""
    global _MATH_SCOPE_CACHE
    if _MATH_SCOPE_CACHE is not None:
        return _MATH_SCOPE_CACHE
    from symbolite.abstract import real as sr

    def _log(*a: Any) -> Any:
        # libSBML L3 emits `log10(x)` for a base-less log and `log(base, x)`
        # for an explicit base; support both. (A single-arg `log` is base-10.)
        if len(a) == 1:
            return sr.log10(a[0])
        base, x = a
        return sr.log(x) / sr.log(base)

    def _root(n: Any, x: Any) -> Any:  # root(degree, radicand)
        return x ** (1.0 / n)

    def _pow(a: Any, b: Any) -> Any:
        return a**b

    def _min(a: Any, b: Any) -> Any:  # symbolic min via abs identity
        return (a + b - sr.abs(a - b)) / 2

    def _max(a: Any, b: Any) -> Any:
        return (a + b + sr.abs(a - b)) / 2

    scope: dict[str, Any] = {
        "exp": sr.exp, "ln": sr.log, "log": _log, "log10": sr.log10, "log2": sr.log2,
        "sqrt": sr.sqrt, "abs": sr.abs, "floor": sr.floor,
        "ceil": sr.ceil, "ceiling": sr.ceil, "factorial": sr.factorial,
        "root": _root, "pow": _pow, "min": _min, "max": _max,
        "sin": sr.sin, "cos": sr.cos, "tan": sr.tan,
        "asin": sr.asin, "acos": sr.acos, "atan": sr.atan, "atan2": sr.atan2,
        "sinh": sr.sinh, "cosh": sr.cosh, "tanh": sr.tanh,
        "asinh": sr.asinh, "acosh": sr.acosh, "atanh": sr.atanh,
        "piecewise": _piecewise,
        "pi": sr.pi, "exponentiale": sr.e, "e": sr.e,
        "true": 1.0, "false": 0.0, "True": 1.0, "False": 0.0,
    }
    _MATH_SCOPE_CACHE = scope
    return scope


class _BoolToArith(ast.NodeTransformer):
    """Rewrite comparisons / boolean ops into arithmetic on Real 0/1 indicators.

    ``a < b``           -> ``__pw_ind__(a < b)``        (Boolean -> Real 0/1)
    ``x and y``         -> ``x * y``                    (AND  = product)
    ``x or y``          -> ``__pw_or__(x, y)``          (OR   = x+y-x*y)
    ``not x``           -> ``1.0 - x``                  (NOT  = complement)

    This avoids symbolite's symbolic boolean operators, which simbio's backend
    cannot translate, while preserving the numeric semantics the solver needs.
    """

    def visit_Compare(self, node: ast.Compare) -> ast.AST:
        self.generic_visit(node)
        return ast.Call(
            func=ast.Name(id="__pw_ind__", ctx=ast.Load()), args=[node], keywords=[]
        )

    def visit_BoolOp(self, node: ast.BoolOp) -> ast.AST:
        self.generic_visit(node)
        values = node.values
        if isinstance(node.op, ast.And):
            expr: ast.AST = values[0]
            for val in values[1:]:
                expr = ast.BinOp(left=expr, op=ast.Mult(), right=val)
            return expr
        expr = values[0]
        for val in values[1:]:
            expr = ast.Call(
                func=ast.Name(id="__pw_or__", ctx=ast.Load()),
                args=[expr, val],
                keywords=[],
            )
        return expr

    def visit_UnaryOp(self, node: ast.UnaryOp) -> ast.AST:
        self.generic_visit(node)
        if isinstance(node.op, ast.Not):
            return ast.BinOp(left=ast.Constant(1.0), op=ast.Sub(), right=node.operand)
        return node


# Formulas containing any of these need the logical/piecewise path; everything
# else uses the fast direct-eval path (no behaviour change, no regressions).
_NEEDS_BOOL = re.compile(r"piecewise|[<>!&|=]")


def _eval_rate_law(formula: str, names: dict[str, Any]):
    """Evaluate a libSBML L3 infix formula into a symbolite expression.

    ``names`` maps every species / parameter / compartment id to its simbio
    object (or a float for a fixed compartment size). simbio's `Species` and
    `Parameter` overload the arithmetic operators, so ordinary Python
    evaluation of the formula builds the symbolic rate law. Math builtins,
    `piecewise`, and user-defined functions (passed in ``names``) extend that
    scope; logical/relational operators are rewritten into Real arithmetic.
    """
    python_expr = formula.replace("^", "**")
    scope = {**_math_scope(), **names}
    if not _NEEDS_BOOL.search(python_expr):
        return eval(python_expr, {"__builtins__": {}}, scope)

    # Logical / relational / piecewise path. SBML logical operators (`&&`, `||`,
    # `!`) are not valid Python, so translate them to keywords, then rewrite
    # comparisons and boolean ops into arithmetic over Real 0/1 indicators.
    translated = python_expr.replace("&&", " and ").replace("||", " or ")
    translated = re.sub(r"!(?!=)", " not ", translated)
    try:
        tree: Any = ast.parse(translated, mode="eval")
    except SyntaxError as exc:
        raise ValueError(
            f"Unsupported construct in SBML rate law {formula!r}: {exc}"
        ) from exc
    tree = _BoolToArith().visit(tree)
    ast.fix_missing_locations(tree)
    code = compile(tree, "<sbml-rate-law>", "eval")
    scope = {**scope, "__pw_ind__": _pw_indicator, "__pw_or__": _pw_or}
    return eval(code, {"__builtins__": {}}, scope)


def model_from_sbml(sbml: str, name: str = "model"):
    """Build a simbio `Compartment` subclass from an SBML string.

    Returns ``(Model, species_names, parameter_names)``.
    """
    import libsbml
    from poincare import Independent
    from poincare.reactions import Reactant, RateLaw
    from simbio import Compartment, Parameter, Volume, assign
    from simbio import initial as initial_field  # aliased: `initial` is reused as
    # a local float in the species loop below (amount/concentration value).
    from simbio.core import amount, concentration, volume as volume_field

    document = libsbml.readSBMLFromString(sbml)
    model = document.getModel()
    if model is None:
        raise ValueError("libSBML could not parse the model")

    compartment_size = {
        c.getId(): (c.getSize() if c.isSetSize() else 1.0)
        for c in model.getListOfCompartments()
    }

    # Assignment rules define a quantity as a formula of others (e.g.
    # CT = C2 + CP + M + pM). They are NOT state variables; we inline their
    # expression wherever they're referenced rather than integrating them.
    # Without this, such a "species" stays a free variable at initial 0 and a
    # rate law like M/CT evaluates 0/0 at compile time (BIOMD0000000005).
    assignment_rules = {
        rule.getVariable(): rule.getMath()
        for rule in model.getListOfRules()
        if rule.isAssignment()
    }
    assignment_vars = set(assignment_rules)

    # Rate rules define an explicit ODE dX/dt = f(...) for a species, parameter,
    # or compartment. These are the model's dynamics for many BioModels (species
    # declared boundaryCondition=true whose change is given entirely by a
    # rateRule, e.g. BIOMD0000000006 / Goldbeter 1991). The target must become an
    # INTEGRATED state variable driven by this derivative — not a fixed constant
    # and not a reaction. Without handling these, such models emit flat constants
    # (or the integrator fails on "no dynamics"). Resolved below, after
    # assignment rules are inlined so RHS references resolve.
    rate_rules = {
        rule.getVariable(): rule.getMath()
        for rule in model.getListOfRules()
        if rule.isRate()
    }
    rate_rule_vars = set(rate_rules)

    namespace: dict[str, Any] = {"__annotations__": {}}
    volume = next(iter(compartment_size.values()), 1.0)
    namespace["volume"] = volume_field(default=float(volume))
    namespace["__annotations__"]["volume"] = Volume

    eval_names: dict[str, Any] = {}
    species_objs: dict[str, Any] = {}
    species_names: list[str] = []
    # Per-species compartment volume + whether the species is tracked as a
    # concentration. SBML rate laws are written in substance/time and (unless
    # `hasOnlySubstanceUnits`) reference species by CONCENTRATION, so we track
    # such species as concentrations and convert an `initialAmount` to a
    # concentration via the compartment volume. The extensive rate law is then
    # divided by that volume per reaction (below) to yield dConc/dt — without
    # this, models whose compartment volume != 1 diverge from COPASI/Tellurium
    # (the amount↔concentration scaling bug; BIOMD1/2/9 had nRMSE ~0.4).
    species_volume: dict[str, float] = {}
    concentration_species: set[str] = set()
    # Boundary / constant species are referenced by rate laws but are NOT
    # consumed by reactions (their value is held fixed unless a rule changes it).
    fixed_species: set[str] = set()
    for s in model.getListOfSpecies():
        sid = s.getId()
        if sid in assignment_vars:
            continue  # derived; inlined below
        vol_s = float(compartment_size.get(s.getCompartment(), volume))
        species_volume[sid] = vol_s
        substance_units = s.getHasOnlySubstanceUnits()
        if s.isSetInitialConcentration():
            conc0 = s.getInitialConcentration()
            amt0 = conc0 * vol_s
        elif s.isSetInitialAmount():
            amt0 = s.getInitialAmount()
            conc0 = (amt0 / vol_s) if vol_s else amt0
        else:
            conc0 = amt0 = 0.0
        if substance_units:
            initial, factory = amt0, amount
        else:
            initial, factory = conc0, concentration
            concentration_species.add(sid)
        if initial is None or (isinstance(initial, float) and math.isnan(initial)):
            initial = 0.0
        obj = factory(default=float(initial))
        namespace[sid] = obj
        namespace["__annotations__"][sid] = type(obj)
        species_objs[sid] = obj
        eval_names[sid] = obj
        species_names.append(sid)
        # Boundary/constant species, and any species governed by a rate rule,
        # are not consumed by reactions (SBML: a rate-rule species' change is
        # given entirely by its rule). Rate-rule species still get an explicit
        # derivative equation below — they must not be left as free constants.
        if s.getBoundaryCondition() or s.getConstant() or sid in rate_rule_vars:
            fixed_species.add(sid)

    parameter_names: list[str] = []
    for p in model.getListOfParameters():
        pid = p.getId()
        if pid in assignment_vars:
            continue  # derived; inlined below
        value = p.getValue() if p.isSetValue() else 0.0
        if pid in rate_rule_vars:
            # A parameter driven by a rate rule is a state variable, not a
            # constant: build it with `initial(...)` so it is integrable and
            # attach its derivative equation below.
            param = initial_field(default=float(value))
            namespace[pid] = param
            namespace["__annotations__"][pid] = type(param)
        else:
            param = assign(default=float(value))
            namespace[pid] = param
            namespace["__annotations__"][pid] = Parameter
        eval_names[pid] = param
        parameter_names.append(pid)

    for cid, size in compartment_size.items():
        eval_names.setdefault(cid, float(size))

    # The SBML `time` csymbol (used by time-dependent rate laws / piecewise
    # pulses) maps to poincare's independent variable. Only register it if no
    # model entity already claims the name.
    if "time" not in eval_names:
        time_var = Independent()
        namespace["time"] = time_var
        namespace["__annotations__"]["time"] = Independent
        eval_names["time"] = time_var

    # SBML functionDefinitions are reusable MathML lambdas referenced by name in
    # kinetic laws and rules (e.g. `function_1(x, y)`, `Constant_flux__irreversible(v)`).
    # Build a Python callable per definition that substitutes the call arguments
    # into the body expression. Bodies may reference math builtins and OTHER
    # functionDefinitions, so each body is evaluated against the shared
    # `func_defs` dict — late binding handles forward references and ordering.
    # (SBML lambdas reference only their own arguments, not model variables.)
    func_defs: dict[str, Any] = {}
    for fd in model.getListOfFunctionDefinitions():
        body = fd.getBody()
        if body is None:
            continue
        arg_names = [fd.getArgument(i).getName() for i in range(fd.getNumArguments())]
        body_formula = libsbml.formulaToL3String(body)

        def _make(arg_names: list[str], body_formula: str):
            def _fn(*call_args: Any) -> Any:
                local = dict(zip(arg_names, call_args))
                return _eval_rate_law(body_formula, {**func_defs, **local})

            return _fn

        func_defs[fd.getId()] = _make(arg_names, body_formula)
    eval_names.update(func_defs)

    # Resolve assignment-rule expressions into symbolic expressions, inlining
    # them into eval_names. Iterate so rules that reference other rules resolve
    # once their dependencies are available.
    pending = dict(assignment_rules)
    while pending:
        progressed = False
        for var, math_node in list(pending.items()):
            try:
                expr = _eval_rate_law(libsbml.formulaToL3String(math_node), eval_names)
            except (KeyError, NameError):
                continue  # depends on an as-yet-unresolved assignment var
            eval_names[var] = expr
            del pending[var]
            progressed = True
        if not progressed:
            # Cyclic or references something unknown: fall back to 0.0 so the
            # build doesn't hang/crash (these vars are not integrated anyway).
            for var in pending:
                eval_names.setdefault(var, 0.0)
            break

    # Attach explicit ODEs for rate-rule targets. The RHS is resolved in the
    # same scope used for kinetic laws (species Variables, parameters, inlined
    # assignment-rule expressions, and compartment sizes) — crucially AFTER the
    # assignment rules above are inlined into eval_names, so references resolve.
    # A rate rule gives the derivative directly in the target's own units
    # (concentration or amount), so — unlike a reaction kinetic law — NO
    # volume division is applied here.
    for var, math_node in rate_rules.items():
        target = eval_names.get(var)
        if not hasattr(target, "derive"):
            # Not yet an integrable Variable (e.g. a compartment tracked as a
            # float size). Promote it to a state Variable so it can be driven.
            if isinstance(target, (int, float)):
                default = float(target)
            else:
                default = float(compartment_size.get(var, 0.0))
            target = initial_field(default=default)
            namespace[var] = target
            namespace["__annotations__"][var] = type(target)
            eval_names[var] = target
        rhs = _eval_rate_law(libsbml.formulaToL3String(math_node), eval_names)
        namespace[f"_rate_{var}"] = target.derive() << rhs

    for r in model.getListOfReactions():
        kinetic = r.getKineticLaw()
        if kinetic is None:
            raise NotImplementedError(f"Reaction {r.getId()} has no kinetic law")
        # Register local (reaction-scoped) parameters as model Parameters under
        # namespaced names, so they stay symbolic in the rate law. Substituting
        # them as plain floats would let an all-literal sub-expression evaluate
        # at build time — e.g. a kinetic law that reduces to 0.0/0.0 raises
        # ZeroDivisionError during compilation (BIOMD0000000005).
        local: dict[str, Any] = {}
        for p in kinetic.getListOfParameters():
            pid = p.getId()
            namespaced = f"{r.getId()}__{pid}"
            param = assign(default=float(p.getValue()) if p.isSetValue() else 0.0)
            namespace[namespaced] = param
            namespace["__annotations__"][namespaced] = Parameter
            local[pid] = param
        scope = {**eval_names, **local}
        rate_law = _eval_rate_law(libsbml.formulaToL3String(kinetic.getMath()), scope)

        def _reactants(refs):
            # Only true state species become reactants. Boundary / constant /
            # assignment-rule species are not consumed by reactions (SBML
            # boundaryCondition semantics), and skipping them avoids KeyErrors
            # for species that were never registered as state variables.
            out = []
            for sr in refs:
                sid = sr.getSpecies()
                if sid in species_objs and sid not in fixed_species:
                    out.append(Reactant(species_objs[sid], sr.getStoichiometry()))
            return out

        reactants = _reactants(r.getListOfReactants())
        products = _reactants(r.getListOfProducts())

        # SBML kinetic laws give the reaction rate in substance/time. simbio
        # applies `rate_law` to a concentration species as dConc/dt, so convert
        # by dividing the extensive rate by the species' compartment volume. We
        # can only do this with a single scalar divisor, so apply it when every
        # state species the reaction touches is a concentration sharing one
        # volume (the single-compartment case). Mixed-volume or amount-tracked
        # reactions keep the raw rate (volume==1 makes the division a no-op for
        # the common case, so nothing regresses there).
        touched = [sr.getSpecies() for sr in
                   list(r.getListOfReactants()) + list(r.getListOfProducts())
                   if sr.getSpecies() in species_objs
                   and sr.getSpecies() not in fixed_species]
        vols = {species_volume[sid] for sid in touched}
        all_conc = all(sid in concentration_species for sid in touched)
        if touched and all_conc and len(vols) == 1:
            divisor = next(iter(vols))
            if divisor not in (0.0, 1.0):
                rate_law = rate_law / float(divisor)

        namespace[f"_rxn_{r.getId()}"] = RateLaw(
            reactants=reactants, products=products, rate_law=rate_law
        )

    model_cls = type(name, (Compartment,), namespace)
    return model_cls, species_names, parameter_names


def model_from_antimony(antimony: str, *, model_name: str | None = None):
    """Compile Antimony -> SBML -> simbio model.

    Returns ``(Model, species_names, parameter_names)``.
    """
    sbml, name = antimony_to_sbml(antimony, model_name=model_name)
    return model_from_sbml(sbml, name=name)


def _looks_like_sbml(text: str) -> bool:
    head = text.lstrip()[:200].lower()
    return head.startswith("<?xml") or "<sbml" in head


def load_model_source(model_source: str, *, model_format: str = "auto"):
    """Build a simbio model from an SBML/Antimony file path, URL, or string.

    ``model_format`` is ``"auto"`` (sniff), ``"sbml"``, or ``"antimony"``. This
    is the entry point the canonical model-source Steps/Processes use, mirroring
    pbg-copasi's ``model_source`` and pbg-tellurium's ``model``/``model_file``.

    Returns ``(Model, species_names, parameter_names)``.
    """
    import os
    from urllib.request import urlopen

    text = model_source
    if model_source.startswith(("http://", "https://")):
        with urlopen(model_source) as response:  # noqa: S310 - user-supplied model URL
            text = response.read().decode("utf-8")
    elif os.path.exists(model_source):
        with open(model_source, encoding="utf-8") as handle:
            text = handle.read()

    fmt = model_format
    if fmt == "auto":
        fmt = "sbml" if _looks_like_sbml(text) else "antimony"

    if fmt == "sbml":
        return model_from_sbml(text)
    return model_from_antimony(text)
