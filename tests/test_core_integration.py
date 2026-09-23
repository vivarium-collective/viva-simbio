"""Guard: build_core() must register this workspace's own processes.

These processes live INSIDE the editable-installed ``viva_simbio`` package, which
``process_bigraph.allocate_core()`` does not auto-discover (an editable install
is invisible to ``importlib.metadata.packages_distributions()``). If the explicit
registration in ``viva_simbio.core`` is ever dropped, every composite that
addresses ``local:SimbioProcess`` breaks deep inside Composite construction with
a cryptic ``no link found at address`` error. This test makes that regression
fail fast with a clear message instead.
"""

from viva_simbio import build_core
from viva_simbio.core import SIMBIO_PROCESSES


def test_build_core_succeeds():
    assert build_core() is not None


def test_workspace_processes_registered():
    core = build_core()
    missing = [name for name in SIMBIO_PROCESSES if name not in core.link_registry]
    assert not missing, (
        f"build_core() did not register workspace processes {missing}; "
        f"composites addressing `local:<name>` will fail to resolve"
    )


def test_simbio_process_resolves_by_address():
    """The exact failure mode the bug produced: local:SimbioProcess must resolve."""
    core = build_core()
    assert core.link_registry.get("SimbioProcess") is not None


def test_numeric_result_type_registered():
    core = build_core()
    assert core.access("numeric_result") is not None
