"""build_core() — the workspace core that the composites and tests run in.

`process_bigraph.allocate_core()` auto-discovers processes/types from *installed
distributions* that depend on bigraph-schema (it scans
``importlib.metadata.packages_distributions()``). That works for sibling
``pbg-*`` wrapper packages installed as regular wheels, but it does **not** see
this workspace's own package when it is installed *editable* (``pip install -e``):
an editable install records only a ``.pth`` shim in its ``RECORD``, so
``packages_distributions()`` never maps ``viva_simbio`` back to the ``viva-simbio``
distribution and the discovery scan skips it entirely. The result is that
``SimbioProcess`` is missing from the core's link registry and any composite that
references ``local:SimbioProcess`` fails with::

    Exception: no link found at address: {'protocol': 'local', 'data': 'SimbioProcess'}

So we register this workspace's *own* processes explicitly here. The
registrations are idempotent (skipped if discovery already provided them, e.g.
in a non-editable / Docker install), so ``build_core()`` is correct regardless
of how ``viva-simbio`` was installed.
"""

from __future__ import annotations

from process_bigraph import allocate_core

from .processes import (
    BaseSimbioStep,
    SimbioProcess,
    SimbioSteadyStateStep,
    SimbioUTCProcess,
    SimbioUTCStep,
)
from .types import register_simbio_types

# The process/step wrappers defined IN this workspace package. Keyed by the
# short name used in composite addresses (``local:<name>``).
SIMBIO_PROCESSES = {
    "SimbioProcess": SimbioProcess,
    "SimbioUTCProcess": SimbioUTCProcess,
    "BaseSimbioStep": BaseSimbioStep,
    "SimbioUTCStep": SimbioUTCStep,
    "SimbioSteadyStateStep": SimbioSteadyStateStep,
}


def register_simbio_processes(core):
    """Register viva-simbio's own process/step wrappers into a core.

    Idempotent: a wrapper already present (e.g. provided by auto-discovery in a
    non-editable install) is left untouched.
    """
    for name, cls in SIMBIO_PROCESSES.items():
        if name not in core.link_registry:
            core.register_link(name, cls)
    return core


def build_core(core=None):
    """Return a process-bigraph core with all viva-simbio types + processes registered.

    This is the canonical core for the workspace: composites that address
    ``local:SimbioProcess`` (and the test suite) must build their `Composite`
    against a core returned from here, not a bare ``allocate_core()``.
    """
    if core is None:
        core = allocate_core()
    register_simbio_types(core)
    register_simbio_processes(core)
    return core
