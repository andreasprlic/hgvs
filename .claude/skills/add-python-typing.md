# Skill: Add Python Typing to an Existing Codebase

Add modern Python 3.11+ type annotations to an existing project, including
the PEP 561 marker, handling common pitfalls along the way.

---

## Trigger

Use this skill when asked to:
- "Add typing to this project"
- "Add type annotations"
- "Make this project typed / mypy-clean / pyright-clean"

---

## Workflow

### 1. Audit first

Use an Explore agent to survey all source files and report:
- Which functions/methods lack parameter or return type annotations
- Which files already have partial annotations
- What style is already in use (modern `X | Y` vs legacy `Optional[X]`)

### 2. Add annotations file-by-file

Work through files in dependency order (data model first, then consumers):
```
enums → location/edit/posedit → sequencevariant
→ dataproviders/interface → variantmapper → assemblymapper
→ validator → normalizer → projector
```

Use **Python 3.11+ syntax** throughout — the project requires Python ≥ 3.11:
| Old style | New style |
|-----------|-----------|
| `Optional[str]` | `str \| None` |
| `Union[X, Y]` | `X \| Y` |
| `List[str]` | `list[str]` |
| `Dict[str, int]` | `dict[str, int]` |
| `Tuple[X, Y]` | `tuple[X, Y]` |
| `from typing import List, Dict, ...` | remove; use builtins |

Only import from `typing` when needed: `Any`, `cast`, `TYPE_CHECKING`, etc.

### 3. Add the PEP 561 marker

```bash
touch src/hgvs/py.typed          # empty marker file
```

Then add it to `pyproject.toml` package-data so it ships in wheels:
```toml
[tool.setuptools.package-data]
"*" = ["_data/*", "py.typed"]
```

Without `py.typed`, type checkers silently ignore all annotations when the
package is installed as a dependency.

### 4. Run tests

```bash
uv run pytest tests/ -x -q --no-cov -m "not network and not slow and not extra"
```

All tests must pass before committing.

---

## Common Pitfalls & Fixes

### `attrs` classes — dynamic attributes unknown to Pyright

`@attr.s(slots=True)` generates `__init__` and `__attrs_attrs__` at runtime.
Pyright reports "No parameter named 'base'" or "Attribute '__attrs_attrs__'
is unknown". **These are pre-existing — do not try to fix them.**

### Config / `__getattr__` returning mixed types

A config class whose `__getattr__` returns `str | int | bool | None`
depending on the INI value will cause cascading errors everywhere a config
value is used as a function default. Fix: annotate `__getattr__` as `-> Any`.

```python
# config.py
from typing import Any

class ConfigGroup:
    def __getattr__(self, k: str) -> Any:   # ← not str, not bool
        return _val_xform(self.__dict__["_section"][k])
```

This prevents Pyright from widening every parameter that has a config default
to `int | Any | bool | None`.

### Circular imports — use `from __future__ import annotations`

When two modules reference each other's types (e.g. `variantmapper` ↔
`sequencevariant`), add this as the **very first** import:

```python
from __future__ import annotations
```

This makes all annotations strings at runtime, breaking the cycle. You can
then import the module (not the class) for isinstance checks.

### Data provider parameter — use Interface, not Any

The `hdp` parameter accepted by mappers, validators, and normalizers should
be typed as the concrete interface type, not `Any`:

```python
import hgvs.dataproviders.interface

def __init__(self, hdp: hgvs.dataproviders.interface.Interface, ...) -> None:
```

### Subclass overrides with fewer parameters (Liskov violations)

When a subclass intentionally drops parameters the base class has (e.g.
`AssemblyMapper` absorbs `alt_aln_method` into `self`), Pyright flags
`reportIncompatibleMethodOverride`. Suppress with a targeted ignore:

```python
def g_to_c(self, var_g, tx_ac, alt_aln_method=None) -> SequenceVariant:  # type: ignore[override]
```

Document *why* the override is intentionally narrower.

### Functions that accept `float("inf")` as a sentinel

Annotating `start: int` then calling with `float("inf") - 1` (which is
`float`) will error. Use `int | float`:

```python
def _fetch_bounded_seq(self, var, start: int | float, end: int | float, ...) -> str:
```

### Distinguish errors you introduced from pre-existing ones

When Pyright reports errors after your edits, sort them into two buckets:

**Pre-existing (ignore):**
- `reportMissingImports` for packages not in Pyright's env (`attr`, `bioutils`, `IPython`)
- `reportSelfClsParameterName` for comparison operators using `lhs` instead of `self`
- `__attrs_attrs__` / attrs constructor parameter errors
- `hgvs.edit.X` / `hgvs.sequencevariant.X` not known — implicit side-effect imports

**Your responsibility (fix):**
- Errors that appear only in files you annotated
- Type mismatches between your new annotations and actual usage
- Missing imports for modules you reference in annotations

---

## Commit strategy

Make small, focused commits:
1. Annotation pass across source files
2. `py.typed` + `pyproject.toml`
3. Targeted fixes for errors introduced by annotations

Each commit should leave the test suite green.
