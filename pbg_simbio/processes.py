"""Process-bigraph wrapper for `simbio` Chemical Reaction Network simulation.

This is a **real bridge**: ``SimbioProcess.update`` builds a genuine
``simbio.Compartment`` from the reaction spec in its config and integrates it
with simbio's LSODA solver (via `poincare`). No reaction kinetics are
reimplemented here — the ODEs are assembled and solved by simbio itself.

simbio: https://github.com/dyscolab/simbio  (extends `poincare`)
"""

from __future__ import annotations

from typing import Any

import numpy as np
from process_bigraph import Process, Step


def _normalize_reactant(entry: Any) -> tuple[str, int]:
    """Accept ``"A"``, ``["A", 2]`` or ``{"species": "A", "stoichiometry": 2}``.

    Returns ``(species_name, stoichiometry)``.
    """
    if isinstance(entry, str):
        return entry, 1
    if isinstance(entry, dict):
        return entry["species"], int(entry.get("stoichiometry", 1))
    if isinstance(entry, (list, tuple)):
        name = entry[0]
        stoich = int(entry[1]) if len(entry) > 1 else 1
        return name, stoich
    raise TypeError(f"Cannot interpret reactant {entry!r}")


def build_crn_model(
    species: dict[str, float],
    reactions: list[dict[str, Any]],
    volume: float = 1.0,
):
    """Assemble a `simbio.Compartment` subclass from a declarative CRN spec.

    Parameters
    ----------
    species
        Mapping of species name -> initial concentration.
    reactions
        Each reaction is a dict with ``reactants``, ``products`` (lists of
        species names, ``[name, stoich]`` pairs, or ``{"species", "stoichiometry"}``
        dicts), a numeric ``rate`` (mass-action rate constant), and an optional
        ``name`` (defaults to ``r0``, ``r1``, ...).
    volume
        Compartment volume.

    Returns
    -------
    (Model, reaction_names)
        The dynamically-created `Compartment` subclass and the ordered list of
        reaction names (whose rate parameters can be overridden at solve time).
    """
    from simbio import Compartment, MassAction, Parameter, Species, Volume, assign
    from simbio.core import concentration, volume as volume_field

    namespace: dict[str, Any] = {"__annotations__": {}}
    namespace["volume"] = volume_field(default=float(volume))
    namespace["__annotations__"]["volume"] = Volume

    species_objs: dict[str, Any] = {}
    for name, initial in species.items():
        obj = concentration(default=float(initial))
        namespace[name] = obj
        namespace["__annotations__"][name] = Species
        species_objs[name] = obj

    reaction_names: list[str] = []
    for index, rxn in enumerate(reactions):
        name = rxn.get("name", f"r{index}")
        reaction_names.append(name)

        rate_param = f"_rate_{name}"
        namespace[rate_param] = assign(default=float(rxn.get("rate", 1.0)))
        namespace["__annotations__"][rate_param] = Parameter

        def _terms(entries):
            terms = []
            for entry in entries:
                sp_name, stoich = _normalize_reactant(entry)
                if sp_name not in species_objs:
                    raise KeyError(
                        f"Reaction references unknown species {sp_name!r}; "
                        f"declare it in `species`."
                    )
                term = species_objs[sp_name]
                terms.append(stoich * term if stoich != 1 else term)
            return terms

        namespace[f"eq_{name}"] = MassAction(
            reactants=_terms(rxn.get("reactants", [])),
            products=_terms(rxn.get("products", [])),
            rate=namespace[rate_param],
        )

    model = type("SimbioCRNModel", (Compartment,), namespace)
    return model, reaction_names


class SimbioProcess(Process):
    """Time-stepped bridge to a simbio Chemical Reaction Network.

    Two equivalent ways to define the network:

    * **Antimony** — set ``config["antimony"]`` to an Antimony model string. It
      is compiled to SBML (libantimony), extracted (libSBML), and rebuilt with
      simbio's core (see :mod:`pbg_simbio.antimony_loader`). Arbitrary rate
      laws (mass-action, Hill, ...) are supported.
    * **Reaction spec** — set ``config["species"]`` + ``config["reactions"]``
      for a mass-action network built directly via :func:`build_crn_model`.

    The surrounding bigraph owns the absolute species concentrations in a
    shared store. Each step the process reads those concentrations, integrates
    the network forward by ``interval`` with simbio, and emits the **change** in
    each concentration as a ``map[string,float]`` delta — so a sibling process
    (an influx, dilution, transport term, or another reaction module) writing
    the same store composes additively.

    Inputs
    ------
    concentrations : map[string,float]
        Current absolute concentration of each species (used as the integration
        initial condition). Missing species fall back to their model initial.
    parameters : map[string,float]
        Optional overrides for the model's parameters (rate constants, Hill
        coefficients, ...), keyed by parameter name. A controller, a
        temperature model, or a slow regulatory process can write these.

    Outputs
    -------
    concentrations : map[string,float]
        Per-species concentration delta over the interval (accumulates
        additively in the shared store).
    """

    config_schema = {
        # Antimony model string (primary path). If empty, species+reactions are used.
        "antimony": {"_type": "string", "_default": ""},
        "species": {"_type": "map[float]", "_default": {}},
        # List of reaction dicts; kept permissive so arbitrary CRN topologies
        # pass through bigraph-schema unchanged.
        "reactions": {"_type": "list", "_default": []},
        "volume": {"_type": "float", "_default": 1.0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        self._model = None
        self._simulator = None
        self._species_names: list[str] = []
        self._param_names: list[str] = []
        self._species_initials: dict[str, float] = {}
        # For the reaction-spec path, reaction names carry _rate_<name> params.
        self._reaction_names: list[str] = []
        # Pre-compute the initial state so initial_state() works before _build().
        if not self.config["antimony"]:
            self._species_initials = dict(self.config["species"])

    def _build(self):
        from simbio import Parameter, Simulator

        if self.config["antimony"]:
            from .antimony_loader import model_from_antimony

            self._model, self._species_names, self._param_names = model_from_antimony(
                self.config["antimony"]
            )
            self._species_initials = {
                name: float(getattr(self._model, name).initial or 0.0)
                for name in self._species_names
            }
        else:
            self._model, self._reaction_names = build_crn_model(
                species=self.config["species"],
                reactions=self.config["reactions"],
                volume=self.config["volume"],
            )
            self._species_names = list(self.config["species"])
            self._species_initials = dict(self.config["species"])
            # reaction-spec rate constants are exposed as _rate_<name> params
            self._param_names = list(self._reaction_names)

        self._simulator = Simulator(self._model)

    def _param_object(self, name):
        """Resolve a parameter name to its simbio Parameter on the model."""
        for attr in (name, f"_rate_{name}"):
            obj = getattr(self._model, attr, None)
            if obj is not None:
                return obj
        return None

    def inputs(self):
        return {
            "concentrations": "map[string,float]",
            "parameters": "map[string,float]",
        }

    def outputs(self):
        return {"concentrations": "map[string,float]"}

    def initial_state(self):
        if self._simulator is None:
            self._build()
        return {"concentrations": dict(self._species_initials)}

    def update(self, state, interval):
        if self._simulator is None:
            self._build()

        current = {**self._species_initials, **(state.get("concentrations") or {})}
        param_overrides = state.get("parameters") or {}

        values: dict[Any, float] = {}
        for name in self._species_names:
            values[getattr(self._model, name)] = float(current[name])
        for name, value in param_overrides.items():
            param = self._param_object(name)
            if param is not None:
                values[param] = float(value)

        result = self._simulator.with_values(values).solve(
            t_span=(0.0, float(interval)),
            save_at=[0.0, float(interval)],
        )

        delta = {}
        for name in self._species_names:
            if name not in result:  # e.g. assignment-rule / non-integrated species
                continue
            final = float(np.asarray(result[name])[-1])
            delta[name] = final - float(current[name])
        return {"concentrations": delta}


# ---------------------------------------------------------------------------
# Canonical model-source Steps / Process (parity with pbg-copasi & pbg-tellurium)
#
# These mirror the process/step "quartet" the other canonical wrappers expose:
#   - a base loader mixin              (BaseSimbioStep   ~ BaseTelluriumStep / BaseCopasi)
#   - a one-shot UTC trajectory Step   (SimbioUTCStep    ~ TelluriumUTCStep / CopasiUTCStep)
#   - a steady-state Step              (SimbioSteadyStateStep ~ *SteadyStateStep)
#   - a time-coupled Process           (SimbioUTCProcess ~ TelluriumProcess / CopasiUTCProcess)
#
# They take a `model_source` (SBML/Antimony file path, URL, or string) and emit
# the canonical `numeric_result` trajectory shape, so they drop straight into
# the pbg-biomodels compare-biomodel composite alongside COPASI and Tellurium.
# ---------------------------------------------------------------------------


class BaseSimbioStep(Step):
    """Abstract base for simbio model-source Steps.

    Loads an SBML/Antimony model (via :func:`load_model_source`) into a genuine
    simbio model + ``Simulator`` and caches species / parameter ids. Subclasses
    implement ``update()`` with the specific task they perform (UTC, steady
    state).
    """

    config_schema = {
        "model_source": {"_type": "string", "_default": ""},
        "model_format": {"_type": "string", "_default": "auto"},
        # Integration tolerances for the LSODA solver. simbio/poincare default
        # to rtol=1e-3, atol=1e-6 — ~1000x looser than COPASI/Tellurium CVODE,
        # which makes stiff models diverge badly. Default to CVODE-comparable
        # tolerances here; tighten further (e.g. rtol=1e-8) for very stiff CRNs.
        "rtol": {"_type": "float", "_default": 1.0e-6},
        "atol": {"_type": "float", "_default": 1.0e-9},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        self._model = None
        self._simulator = None
        self._species_ids: list[str] = []
        self._param_ids: list[str] = []

    def _solver(self):
        """Build the LSODA solver from the configured rtol/atol."""
        from poincare.solvers import LSODA

        return LSODA(rtol=float(self.config["rtol"]),
                     atol=float(self.config["atol"]))

    def _simbio_initialize(self):
        if self._simulator is not None:
            return
        from simbio import Simulator

        from .antimony_loader import load_model_source

        source = self.config["model_source"]
        if not source:
            raise ValueError(
                f"{type(self).__name__} requires a non-empty 'model_source'."
            )
        self._model, self._species_ids, self._param_ids = load_model_source(
            source, model_format=self.config["model_format"]
        )
        self._simulator = Simulator(self._model)

    def _initial_concentrations(self) -> dict[str, float]:
        return {
            sid: float(getattr(self._model, sid).initial or 0.0)
            for sid in self._species_ids
        }

    def _override_values(self, incoming: dict[str, float]) -> dict[Any, float]:
        """Map an incoming {species_id: value} dict to solver `values=` overrides."""
        values: dict[Any, float] = {}
        for sid, value in (incoming or {}).items():
            if sid in self._species_ids:
                values[getattr(self._model, sid)] = float(value)
        return values

    def inputs(self):
        return {}

    def initial_state(self):
        self._simbio_initialize()
        return {"species_concentrations": self._initial_concentrations()}


class SimbioUTCStep(BaseSimbioStep):
    """One-shot uniform-time-course Step returning a dense trajectory.

    Integrates ``model_source`` from 0 to ``time`` at ``n_points`` evenly
    spaced samples and returns the canonical ``numeric_result`` shape
    (``{time, columns, values}``) — matching ``CopasiUTCStep`` /
    ``TelluriumUTCStep``. Optional incoming ``species_concentrations`` /
    ``species_counts`` set the integration initial condition (for coupling).
    """

    config_schema = {
        **BaseSimbioStep.config_schema,
        "time": {"_type": "float", "_default": 10.0},
        "n_points": {"_type": "integer", "_default": 2},
    }

    def inputs(self):
        return {
            "species_concentrations": "map[float]",
            "species_counts": "map[float]",
        }

    def outputs(self):
        return {"result": "numeric_result"}

    def update(self, state):
        self._simbio_initialize()

        n_points = int(self.config["n_points"])
        if n_points < 2:
            raise ValueError("n_points must be >= 2")

        incoming = (
            (state or {}).get("species_counts")
            or (state or {}).get("species_concentrations")
            or {}
        )
        values = self._override_values(incoming)

        save_at = np.linspace(0.0, float(self.config["time"]), n_points)
        result = self._simulator.with_values(values).with_solver(self._solver()).solve(
            save_at=save_at
        )

        # Only species simbio actually integrated appear in the result. SBML
        # species defined by assignment rules (e.g. "Rum1Total"), or otherwise
        # not emitted, are skipped rather than raising KeyError.
        columns = [c for c in self._species_ids if c in result]
        arrays = {c: np.asarray(result[c]) for c in columns}
        # If the model has no integrable state (e.g. every species is boundary
        # or driven by an unsupported rate rule), simbio returns an empty
        # dataset with no `time` coord — fall back to the requested grid so the
        # step yields an empty-but-valid numeric_result instead of crashing.
        if "time" in result.coords:
            times = [float(t) for t in np.asarray(result.coords["time"])]
        else:
            times = [float(t) for t in save_at]
        values_rows = [
            [float(arrays[c][i]) for c in columns] for i in range(len(times))
        ]
        return {
            "result": {
                "time": times,
                "columns": columns,
                "values": values_rows,
            }
        }


class SimbioSteadyStateStep(BaseSimbioStep):
    """Approximate steady-state Step.

    simbio (an ODE integrator) has no algebraic steady-state solver, so this
    integrates to a long ``equilibration_time`` and reports the final
    concentrations as the steady state — honestly an approximation, not a
    Newton solve. Output matches ``TelluriumSteadyStateStep``.
    """

    config_schema = {
        **BaseSimbioStep.config_schema,
        "equilibration_time": {"_type": "float", "_default": 1.0e4},
    }

    def outputs(self):
        return {"steady_state_concentrations": "overwrite[map[float]]"}

    def update(self, state):
        self._simbio_initialize()

        incoming = (
            (state or {}).get("species_counts")
            or (state or {}).get("species_concentrations")
            or {}
        )
        values = self._override_values(incoming)

        t_end = float(self.config["equilibration_time"])
        result = self._simulator.with_values(values).with_solver(self._solver()).solve(
            t_span=(0.0, t_end), save_at=[0.0, t_end],
        )
        steady = {
            sid: float(np.asarray(result[sid])[-1])
            for sid in self._species_ids if sid in result
        }
        return {"steady_state_concentrations": steady}


class SimbioUTCProcess(Process):
    """Time-coupled bridge Process over a ``model_source`` model.

    Advances the integration by ``interval`` on each ``update()`` and reports
    absolute species concentrations and time — matching the role of
    ``TelluriumProcess`` / ``CopasiUTCProcess``. Optional incoming
    ``species_concentrations`` re-seed the state each step (coupling).

    For additive, delta-based composition prefer :class:`SimbioProcess`; this
    class exists for parity with the other canonical Process wrappers, which
    report absolute state via ``overwrite``.
    """

    config_schema = {
        "model_source": {"_type": "string", "_default": ""},
        "model_format": {"_type": "string", "_default": "auto"},
        # See BaseSimbioStep — CVODE-comparable defaults, override for stiff CRNs.
        "rtol": {"_type": "float", "_default": 1.0e-6},
        "atol": {"_type": "float", "_default": 1.0e-9},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        self._model = None
        self._simulator = None
        self._species_ids: list[str] = []
        self._param_ids: list[str] = []
        self._state: dict[str, float] = {}
        self._time = 0.0

    def _solver(self):
        from poincare.solvers import LSODA

        return LSODA(rtol=float(self.config["rtol"]),
                     atol=float(self.config["atol"]))

    def _build(self):
        if self._simulator is not None:
            return
        from simbio import Simulator

        from .antimony_loader import load_model_source

        source = self.config["model_source"]
        if not source:
            raise ValueError("SimbioUTCProcess requires a non-empty 'model_source'.")
        self._model, self._species_ids, self._param_ids = load_model_source(
            source, model_format=self.config["model_format"]
        )
        self._simulator = Simulator(self._model)
        self._state = {
            sid: float(getattr(self._model, sid).initial or 0.0)
            for sid in self._species_ids
        }

    def inputs(self):
        return {"species_concentrations": "maybe[map[float]]"}

    def outputs(self):
        return {
            "species_concentrations": "overwrite[map[float]]",
            "time": "overwrite[float]",
        }

    def initial_state(self):
        self._build()
        return {
            "species_concentrations": dict(self._state),
            "time": self._time,
        }

    def update(self, state, interval):
        self._build()

        incoming = (state or {}).get("species_concentrations")
        if incoming:
            for sid, value in incoming.items():
                if sid in self._state:
                    self._state[sid] = float(value)

        values = {getattr(self._model, sid): self._state[sid] for sid in self._species_ids}
        result = self._simulator.with_values(values).with_solver(self._solver()).solve(
            t_span=(0.0, float(interval)), save_at=[0.0, float(interval)],
        )
        self._state = {
            sid: (float(np.asarray(result[sid])[-1]) if sid in result
                  else self._state.get(sid, 0.0))
            for sid in self._species_ids
        }
        self._time += float(interval)
        return {
            "species_concentrations": dict(self._state),
            "time": self._time,
        }
