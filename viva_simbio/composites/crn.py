"""Composite generators for simbio Chemical Reaction Networks.

The headline models are defined as **Antimony** strings and loaded into genuine
simbio models (see :mod:`viva_simbio.antimony_loader`). Each generator wires a
:class:`SimbioProcess` to a shared ``concentrations`` store **and** a
``parameters`` store (so a sibling can perturb rate constants), plus a RAM
emitter — so the dashboard's Composites tab can run and sweep them.
"""

from __future__ import annotations

from viva_superpowers.composite_generator import composite_generator

# --- Antimony model definitions (shown in the demo report) -------------------

BRUSSELATOR = """\
model brusselator
  # Classic Prigogine oscillator: a chemical clock.
  species X = 1, Y = 1;
  k1 = 1; k2 = 3; k3 = 1; k4 = 1;
  J1:           -> X;   k1;          # constant production of X
  J2: X         -> Y;   k2 * X;      # X converts to Y
  J3: 2 X + Y   -> 3 X; k3 * X^2 * Y;  # autocatalytic feedback
  J4: X         ->;     k4 * X;      # X decays
end
"""

LOTKA_VOLTERRA = """\
model lotka_volterra
  # Predator-prey oscillations as a reaction network.
  # (alpha/beta/gamma/delta are reserved math-function names in Antimony,
  #  so the rate constants are named k_*.)
  species prey = 10, pred = 5;
  k_birth = 1.1; k_predation = 0.4; k_death = 0.4; k_repro = 0.1;
  birth:     -> prey;  k_birth * prey;          # prey reproduce
  predation: prey ->;  k_predation * prey * pred;  # predators eat prey
  death:     pred ->;  k_death * pred;          # predators die
  repro:     -> pred;  k_repro * prey * pred;   # predators reproduce from prey
end
"""

REPRESSILATOR = """\
model repressilator
  # Three genes repressing each other in a ring (Elowitz & Leibler 2000),
  # written with Hill-function transcription -- not mass-action.
  species m1 = 0.2, m2 = 0.1, m3 = 0.3, p1 = 0, p2 = 0, p3 = 0;
  a = 10; a0 = 0.1; n = 3; K = 1; beta = 1;
  tx1: -> m1; a / (1 + (p3/K)^n) + a0;
  tx2: -> m2; a / (1 + (p1/K)^n) + a0;
  tx3: -> m3; a / (1 + (p2/K)^n) + a0;
  dm1: m1 ->; m1;   dm2: m2 ->; m2;   dm3: m3 ->; m3;
  tl1: -> p1; beta * m1;   tl2: -> p2; beta * m2;   tl3: -> p3; beta * m3;
  dp1: p1 ->; beta * p1;   dp2: p2 ->; beta * p2;   dp3: p3 ->; beta * p3;
end
"""


def _antimony_document(*, antimony, species_init, params_init, interval):
    """Wire a SimbioProcess + emitter around shared concentration & parameter stores.

    Both process input ports are wired:
      * ``concentrations`` <-> the ``concentrations`` store (read + delta write)
      * ``parameters``      <- the ``parameters`` store (rate constants etc.)
    """
    return {
        "simbio": {
            "_type": "process",
            "address": "local:SimbioProcess",
            "config": {"antimony": antimony},
            "interval": interval,
            "inputs": {
                "concentrations": ["stores", "concentrations"],
                "parameters": ["stores", "parameters"],
            },
            "outputs": {"concentrations": ["stores", "concentrations"]},
        },
        "stores": {
            "concentrations": dict(species_init),
            "parameters": dict(params_init),
        },
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


@composite_generator(
    name="simbio_brusselator",
    description="Brusselator chemical oscillator (Antimony) integrated by simbio.",
    parameters={
        "k2": {"type": "float", "default": 3.0,
               "description": "Conversion rate X->Y; tune to cross the Hopf bifurcation"},
        "k3": {"type": "float", "default": 1.0,
               "description": "Autocatalytic feedback rate"},
        "interval": {"type": "float", "default": 0.25,
                     "description": "Composite tick size"},
    },
)
def brusselator(core=None, *, k2=3.0, k3=1.0, interval=0.25):
    species_init = {"X": 1.0, "Y": 1.0}
    params_init = {"k1": 1.0, "k2": k2, "k3": k3, "k4": 1.0}
    return _antimony_document(antimony=BRUSSELATOR, species_init=species_init,
                              params_init=params_init, interval=interval)


@composite_generator(
    name="simbio_lotka_volterra",
    description="Lotka-Volterra predator-prey oscillator (Antimony) integrated by simbio.",
    parameters={
        "k_birth": {"type": "float", "default": 1.1, "description": "Prey birth rate"},
        "k_predation": {"type": "float", "default": 0.4, "description": "Predation rate"},
        "k_death": {"type": "float", "default": 0.4, "description": "Predator death rate"},
        "k_repro": {"type": "float", "default": 0.1, "description": "Predator reproduction"},
        "interval": {"type": "float", "default": 0.25,
                     "description": "Composite tick size"},
    },
)
def lotka_volterra(core=None, *, k_birth=1.1, k_predation=0.4, k_death=0.4,
                   k_repro=0.1, interval=0.25):
    species_init = {"prey": 10.0, "pred": 5.0}
    params_init = {"k_birth": k_birth, "k_predation": k_predation,
                   "k_death": k_death, "k_repro": k_repro}
    return _antimony_document(antimony=LOTKA_VOLTERRA, species_init=species_init,
                              params_init=params_init, interval=interval)


@composite_generator(
    name="simbio_repressilator",
    description="Three-gene repressilator with Hill kinetics (Antimony) integrated by simbio.",
    parameters={
        "a": {"type": "float", "default": 10.0,
              "description": "Max transcription rate"},
        "n": {"type": "float", "default": 3.0, "description": "Hill coefficient"},
        "beta": {"type": "float", "default": 1.0,
                 "description": "Translation / degradation rate"},
        "interval": {"type": "float", "default": 0.5,
                     "description": "Composite tick size"},
    },
)
def repressilator(core=None, *, a=10.0, n=3.0, beta=1.0, interval=0.5):
    species_init = {"m1": 0.2, "m2": 0.1, "m3": 0.3, "p1": 0.0, "p2": 0.0, "p3": 0.0}
    params_init = {"a": a, "a0": 0.1, "n": n, "K": 1.0, "beta": beta}
    return _antimony_document(antimony=REPRESSILATOR, species_init=species_init,
                              params_init=params_init, interval=interval)
