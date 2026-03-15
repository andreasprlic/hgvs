# Tech Strategy - Golden Paths

This is the **SINGLE SOURCE OF TRUTH** for technology choices in the hgvs project.

## Compliance

1. **Follow This File**: Use the technologies listed in the Golden Paths below
2. **No Deviations**: Do not suggest alternatives unless explicitly instructed
3. **Latest Stable**: Always use the latest stable version unless pinned

## Project Overview

**hgvs** is a Python bioinformatics library for parsing, formatting, mapping, and validating genome variants using the HGVS nomenclature. Published to PyPI as part of the [biocommons](https://github.com/biocommons) organization.

## Python Stack

| Component | Choice |
|-----------|--------|
| Runtime | Python 3.11–3.13 (`>=3.11, <3.14`) |
| Package Manager | uv |
| Build System | setuptools + setuptools_scm |
| Linting & Formatting | Ruff (line length 100, double quotes, preview mode) |
| Type Checking | mypy |
| Dependency Audit | deptry |
| Testing | pytest + pytest-cov |
| Test Recording | VCRpy + pytest-recording |
| Test Environments | tox + tox-uv |
| Interactive Shell | IPython / Jupyter |

## Key Dependencies

| Dependency | Purpose |
|------------|---------|
| attrs | Data classes and object utilities |
| parsley | PEG parsing for HGVS expressions |
| psycopg2 | PostgreSQL adapter (UTA database) |
| biocommons.seqrepo | Biological sequence repository |
| bioutils | Shared biocommons utilities |
| uta-align | Optional: sequence alignment |

## Database

| Component | Choice |
|-----------|--------|
| Engine | PostgreSQL |
| Adapter | psycopg2 |
| Purpose | UTA (Universal Transcript Archive) — biological data lookup |
| Connection | `postgresql://user:pass@host:5432/uta/uta_20241220` |

## Documentation

| Component | Choice |
|-----------|--------|
| Tool | MkDocs |
| Theme | mkdocs-material |
| API Docs | mkdocstrings[python] |

## CI/CD

| Component | Choice |
|-----------|--------|
| Platform | GitHub Actions (reusable biocommons workflow) |
| Publishing | PyPI via `uv publish` (Trusted Publishing) |
| Trigger | All branches and tags |

## Pre-commit Hooks

- Ruff (check + format)
- Standard file hygiene (trailing whitespace, merge conflicts, JSON/YAML/TOML validation, private key detection)
- .gitignore sorting

## Test Markers

`network`, `slow`, `extra`, `issues`, `mapping`, `models`, `normalization`, `quick`, `regression`, `validation`
