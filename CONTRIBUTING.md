# Contributing to viva-simbio

## Development setup

uv is required. Install with `brew install uv` or `pip install uv`.
simbio needs Python ≥ 3.12.

    uv venv --python 3.12 .venv
    source .venv/bin/activate
    uv pip install -e ".[dev]"
    pytest

## Releasing to PyPI

Tag a commit with `git tag v<VERSION>` and push the tag. The
`.github/workflows/release.yml` workflow publishes to PyPI automatically
using trusted publishing (no tokens needed after initial setup).

PyPI trusted publishing must be configured once per repo. See
https://docs.pypi.org/trusted-publishers/.
