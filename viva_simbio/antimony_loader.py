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

import math
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


def _eval_rate_law(formula: str, names: dict[str, Any]):
    """Evaluate a libSBML L3 infix formula into a symbolite expression.

    ``names`` maps every species / parameter / compartment id to its simbio
    object (or a float for a fixed compartment size). simbio's `Species` and
    `Parameter` overload the arithmetic operators, so ordinary Python
    evaluation of the formula builds the symbolic rate law.
    """
    python_expr = formula.replace("^", "**")
    return eval(python_expr, {"__builtins__": {}}, names)


def model_from_sbml(sbml: str, name: str = "model"):
    """Build a simbio `Compartment` subclass from an SBML string.

    Returns ``(Model, species_names, parameter_names)``.
    """
    import libsbml
    from poincare.reactions import Reactant, RateLaw
    from simbio import Compartment, Parameter, Volume, assign
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
        if s.getBoundaryCondition() or s.getConstant():
            fixed_species.add(sid)

    parameter_names: list[str] = []
    for p in model.getListOfParameters():
        pid = p.getId()
        if pid in assignment_vars:
            continue  # derived; inlined below
        value = p.getValue() if p.isSetValue() else 0.0
        param = assign(default=float(value))
        namespace[pid] = param
        namespace["__annotations__"][pid] = Parameter
        eval_names[pid] = param
        parameter_names.append(pid)

    for cid, size in compartment_size.items():
        eval_names.setdefault(cid, float(size))

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
