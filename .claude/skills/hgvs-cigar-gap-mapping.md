---
name: hgvs-cigar-gap-mapping
description: How to diagnose and fix variant mapping bugs involving CIGAR I/D gap segments — including internal gaps, adjacent gaps, and the double-gap cancellation pattern. Use when g_to_c or g_to_n produces wrong output near an alignment discrepancy.
---

# HGVS CIGAR Gap Mapping

Reference for the three gap-related code paths in `VariantMapper` (`g_to_c` / `g_to_n`) and how to diagnose which one is needed.

---

## CIGAR Op Meanings (hgvs convention)

| Op | Which sequence advances | Meaning |
|----|------------------------|---------|
| `=` / `X` | both | normal match / mismatch |
| `I` | genome (ref) only | genome has a base **not present in transcript** |
| `D` | transcript (tgt) only | transcript has a base **not present in genome** |
| `N` | genome (ref) only | intron — genomic bases skipped in transcript |

> `I` and `D` are easy to confuse. The memory aid: **I = genome Inserts a base** (transcript lacks it). **D = transcript has an eXtra base** (genome lacks it, "D" is the odd one out).

---

## Three Gap Code Paths

### 1. `pos.uncertain = True` (endpoint inside gap)

Triggered when `mapper.g_to_c(pos)` returns `uncertain=True` — at least one variant endpoint maps to an `I` or `D` CIGAR segment.

**Code path**: `else` branch in `g_to_c` / `g_to_n`.

**Fix used**: `_get_altered_tx_sequence` with `pos_g = mapper.c_to_g(pos_c)`.

---

### 2. Internal gap (`_variant_has_internal_gap`)

Both endpoints are in `=` segments (`uncertain=False`), but an `I` or `D` segment lies **strictly inside** the variant interval. The genomic span is wider (or narrower) than the transcript span.

**Diagnosis**: `_variant_has_internal_gap(mapper, var_g)` returns `True`.

**Code path**: `elif _variant_has_internal_gap` in the `not pos.uncertain` branch.

**Fix used**: `_get_altered_tx_sequence` with `pos_g = var_g.posedit.pos` (no expansion needed).

**Example**: `g.119027726_119027728delinsTT` → `c.526_527delinsAA`
- I-segment at 119027727 is inside the 3-base interval
- Genomic: 3 bases, transcript: 2 bases (I-seg excluded)

---

### 3. Adjacent gap — double-gap cancellation (`_variant_spans_i_segment`)

Both endpoints are in `=` segments (`uncertain=False`), but an `I` segment sits **immediately adjacent** (its interbase start = the variant's interbase end, or vice versa). The alt allele's deletion effectively cancels with the adjacent I-segment, yielding a simpler transcript edit.

**Diagnosis**: `_variant_spans_i_segment(mapper, var_g)` returns `True` (uses interbase end, not base-1).

**Code path**: `elif _variant_spans_i_segment` in the `not pos.uncertain` branch.

**Fix used**: `_expand_pos_g_for_adjacent_gap` + `_get_altered_tx_sequence`.

**Example**: `NC_000011.10:g.119027721_119027726delinsTCACA` → `NM_001164277.1:c.532G>A`

```
3-way alignment (minus strand view):
  ref: T G G T G T G G T   (9 bases, includes I-seg at pos 119027727)
  alt: T G - T G T G A T   (1-base deletion in alt at same position as I-seg)
  tx : T - G T G T G G T   (I-seg absent from transcript)

  The alt's deletion and the I-seg "cancel" → only a G→A substitution remains at c.532
```

**Mechanics**:
1. `_expand_pos_g_for_adjacent_gap` extends `pos_g.end.base` by one to include the I-seg base in `seq`
2. `_get_altered_tx_sequence` builds:
   - `tx_ref` = seq excluding I-seg → `CCACAC` → RC → `GTGTGG`
   - `tx_alt` = alt_seq + I-seg base (suffix at `j == var_end`) → `TCACAC` → RC → `GTGTGA`
3. Trim common prefix (5 chars) → `G>A` at c.532

The key line in `_get_altered_tx_sequence` (delins case):
```python
# Include adjacent I-seg at exactly var_end in suffix_tx so the double-gap cancels
suffix_tx = [seq[j] for j in range(var_end, len(seq)) if j not in i_offsets or j == var_end]
```

---

## How to Identify Which Path is Needed

```python
from hgvs.alignmentmapper import AlignmentMapper
import hgvs.variantmapper as vm_module

mapper = AlignmentMapper(hdp, TX_AC, GENOME_AC, 'splign')
var_g.fill_ref(hdp)

pos_c = mapper.g_to_c(var_g.posedit.pos)
print('uncertain:', pos_c.uncertain)

vmap = hgvs.variantmapper.VariantMapper(hdp)
print('internal gap:', vmap._variant_has_internal_gap(mapper, var_g))
print('adjacent gap:', vmap._variant_spans_i_segment(mapper, var_g))
```

Decision:
```
uncertain=True                    → path 1 (existing uncertain logic)
uncertain=False, internal gap     → path 2 (_variant_has_internal_gap)
uncertain=False, adjacent gap     → path 3 (_variant_spans_i_segment + expand)
uncertain=False, neither          → simple strand-flip (no gap involved)
```

---

## Guard for Partially Uncertain Variants

Methods that call `.base` on `var_g.posedit.pos.start` / `.end` will fail for partially uncertain variants whose positions are `Interval` objects (e.g. `g.108337304_(108337428_?)del`). Always guard:

```python
try:
    start_offset = var_g.posedit.pos.start.base - 1 - gc_offset
    end_offset   = var_g.posedit.pos.end.base - gc_offset
except AttributeError:
    return False   # can't determine CIGAR op; skip the gap-fix path
```

This guard is in place in `_variant_spans_i_segment`. Apply the same pattern to any new helper that accesses `.base` on variant positions.

---

## Related Skills

- `hgvs-mapping-verification` — step-by-step verification of expected mapping outputs
- `hgvs-test-cache` — populate UTA cache after adding new test cases
