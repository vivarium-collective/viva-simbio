"""viva-simbio: process-bigraph wrapper for the simbio CRN simulator."""

from .antimony_loader import load_model_source, model_from_antimony, model_from_sbml
from .composites import brusselator, lotka_volterra, repressilator
from .core import build_core, register_simbio_processes
from .processes import (
    BaseSimbioStep,
    SimbioProcess,
    SimbioSteadyStateStep,
    SimbioUTCProcess,
    SimbioUTCStep,
    build_crn_model,
)
from .types import register_simbio_types

__all__ = [
    # composable, Antimony-native CRN process
    "SimbioProcess",
    "build_crn_model",
    # canonical model-source Steps / Process (parity with copasi & tellurium)
    "BaseSimbioStep",
    "SimbioUTCStep",
    "SimbioSteadyStateStep",
    "SimbioUTCProcess",
    # loaders + types
    "model_from_antimony",
    "model_from_sbml",
    "load_model_source",
    "register_simbio_types",
    # workspace core (registers types + this workspace's own processes)
    "build_core",
    "register_simbio_processes",
    # composite generators
    "brusselator",
    "lotka_volterra",
    "repressilator",
]
