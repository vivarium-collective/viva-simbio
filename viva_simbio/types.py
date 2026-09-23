"""Custom bigraph-schema types for viva-simbio.

`numeric_result` is the canonical trajectory shape shared with pbg-copasi /
pbg-tellurium and the pbg-biomodels comparison flow:

    {"time": list[float], "columns": list[string], "values": list[list[float]]}

It is registered defensively (skipped if a host workspace, e.g. pbg-biomodels,
has already registered it) so importing viva-simbio into another workspace does
not clash.
"""

SIMBIO_TYPES = {
    "numeric_result": {
        "time": "list[float]",
        "columns": "list[string]",
        "values": "list[list[float]]",
    },
    "numeric_results": "map[numeric_result]",
}


def register_simbio_types(core):
    """Register viva-simbio bigraph-schema types into a process-bigraph core."""
    try:
        existing = set(core.types())
    except Exception:
        existing = set()
    to_register = {k: v for k, v in SIMBIO_TYPES.items() if k not in existing}
    if to_register:
        core.register_types(to_register)
    return core
