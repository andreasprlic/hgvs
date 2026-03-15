---
name: hgvs-mapping-verification
description: Verify HGVS variant mapping test expectations from first principles before implementing a fix. Use when a plan asserts an expected g_to_c/g_to_n result that you cannot independently confirm, or when a test fails and you need to know whether the code or the expectation is wrong.
---

# HGVS Mapping Verification

Use this skill **before implementing a fix** for a reported wrong mapping. Plans and bug reports sometimes contain incorrect expected outputs. Verify independently from the sequence data.

---

## Step 0 — Variant context investigation

Run this first for any unfamiliar variant. It identifies the gene, which exon the variant falls in, and whether any internal CIGAR gaps (D/I within an exon) are adjacent to the variant — the most common source of edge-case bugs.

```python
import hgvs.parser
import hgvs.dataproviders.uta as uta
import hgvs.assemblymapper

hp  = hgvs.parser.Parser()
hdp = uta.connect()
am  = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38', alt_aln_method='splign', replace_reference=True)

GENOME_AC = 'NC_000011.10'   # adjust as needed
TX_AC     = 'NM_XXXXXX.X'
VAR_G     = 'NC_000011.10:g.X_YdelinsZ'

var_g = hp.parse(VAR_G)

# 1. Full-stack mapping
var_c = am.g_to_c(var_g, TX_AC)
var_p = am.c_to_p(var_c)
print('c.:', var_c)
print('p.:', var_p)
print()

# 2. Transcript / gene info
tx_info = hdp.get_tx_info(TX_AC, GENOME_AC, 'splign')
print('gene:', tx_info['hgnc'])
print('cds_start_i / cds_end_i:', tx_info['cds_start_i'], tx_info['cds_end_i'])
print()

# 3. Exon structure — find which exon(s) the variant overlaps
exons = hdp.get_tx_exons(TX_AC, GENOME_AC, 'splign')
g_start = var_g.posedit.pos.start.base
g_end   = var_g.posedit.pos.end.base

print('Exons overlapping variant:')
for e in exons:
    if e['alt_start_i'] <= g_end and e['alt_end_i'] >= g_start:
        print(f"  exon {e['ord']}: genome {e['alt_start_i']}-{e['alt_end_i']}, "
              f"tx {e['tx_start_i']}-{e['tx_end_i']}, strand {e['alt_strand']}, "
              f"cigar {e['cigar']}")
print()

# 4. Internal CIGAR gap analysis — D/I within an exon (not intron N)
import re
cigar_re = re.compile(r'(\d+)([=DIMNX])')

print('Internal D/I gaps in overlapping exons:')
for e in exons:
    if not (e['alt_start_i'] <= g_end and e['alt_end_i'] >= g_start):
        continue
    strand    = e['alt_strand']
    alt_end   = e['alt_end_i']
    ops       = cigar_re.findall(e['cigar'])
    # walk CIGAR from genomic 3' end (alt_end) on minus strand
    g_cur = alt_end if strand == -1 else e['alt_start_i']
    for length_str, op in (reversed(ops) if strand == -1 else ops):
        length = int(length_str)
        if op in '=MX':
            if strand == -1:
                g_cur -= length
            else:
                g_cur += length
        elif op in 'DI':
            boundary = g_cur  # genomic position of the gap boundary
            dist_start = abs(g_start - boundary)
            dist_end   = abs(g_end   - boundary)
            print(f"  exon {e['ord']}: op={op} at genomic boundary ~{boundary}, "
                  f"dist from var start={dist_start}, dist from var end={dist_end}")
            # D doesn't consume genomic bases; I does
            if op == 'I' and strand == -1:
                g_cur -= length
            elif op == 'I':
                g_cur += length
```

**What to look for:**
- `D` in the exon CIGAR → transcript has bases with no genomic counterpart (variant lands *before* the gap)
- `I` in the exon CIGAR → genome has a base with no transcript counterpart (variant lands *after* the gap on the transcript)
- Distance 0–10 from a gap boundary → likely an edge case; the variant may straddle the gap

> **UTA CIGAR convention** (opposite of standard SAM/BAM): UTA CIGARs are *transcript-centric* — `D` = deletion in genome relative to transcript, `I` = insertion in genome relative to transcript. `build_tx_cigar` + `CIGARMapper` swap this so that genome=ref and transcript=tgt, but the raw `e['cigar']` values from `hdp.get_tx_exons()` follow the UTA convention. See the header of `src/hgvs/utils/cigarmapper.py` for the full explanation.

---

## Step 1 — Confirm what the code actually produces

```python
HGVS_CACHE_MODE=run uv run python -c "
import hgvs.dataproviders.uta as uta
import hgvs.parser, hgvs.assemblymapper, hgvs.variantmapper

hdp = uta.connect()
hp  = hgvs.parser.Parser()
am  = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38')
vm  = hgvs.variantmapper.VariantMapper(hdp)   # no normalisation

var_g = hp.parse_hgvs_variant('NC_000011.10:g.X_YdelinsZ')

# Raw (no normalisation)
print('raw :', vm.g_to_c(var_g, 'NM_XXXXXX.X'))
# Normalised
print('norm:', am.g_to_c(var_g, 'NM_XXXXXX.X'))
"
```

---

## Step 2 — Check whether the uncertain path is triggered

```python
import hgvs.dataproviders.uta as uta
from hgvs.alignmentmapper import AlignmentMapper

hdp    = uta.connect()
mapper = AlignmentMapper(hdp, 'NM_XXXXXX.X', 'NC_000011.10', 'splign')
var_g.fill_ref(hdp)

pos_c = mapper.g_to_c(var_g.posedit.pos)
print('pos_c         :', pos_c)
print('pos_c.uncertain:', pos_c.uncertain)   # True → uncertain path; False → direct mapping
```

If `uncertain=False`, the direct (non-uncertain) path runs and the result is **mechanically correct** — the plan may have the wrong expectation.

---

## Step 3 — Identify CIGAR segment for each variant endpoint

```python
mapper = AlignmentMapper(hdp, 'NM_XXXXXX.X', 'NC_000011.10', 'splign')
gc_offset = mapper.gc_offset
ref_pos   = mapper.cigarmapper.ref_pos
cigar_op  = mapper.cigarmapper.cigar_op

print('gc_offset:', gc_offset)
print('ref_pos  :', ref_pos)
print('cigar_op :', cigar_op)

for label, base in [('start', START_BASE), ('end', END_BASE)]:
    offset = base - 1 - gc_offset          # same formula as _extract_genomic_position
    for i, op in enumerate(cigar_op):
        if offset < ref_pos[i + 1]:
            print(f'{label}: 1-based={base}, cigar_offset={offset}, op={op!r}, '
                  f'segment=[{ref_pos[i]},{ref_pos[i+1]})')
            break
```

Key ops:
- `=` / `X` / `M` — normal match/mismatch → deterministic mapping
- `I` — genome-only base (no transcript equivalent) → may trigger `uncertain=True`
- `D` — transcript-only base (no genome equivalent) → may trigger `uncertain=True`
- `N` — intron

**The formula `base - 1 - gc_offset` is used for BOTH start and end positions** in `_extract_genomic_position`. The CIGAR mapper uses `end="start"` or `end="end"` for intronic offset handling, but the ref-position lookup uses the same arithmetic.

---

## Step 4 — Verify the expected transcript edit manually

For a **minus-strand** transcript:

```python
hdp = uta.connect()

# 1. Fetch genomic forward-strand sequence
g_seq = hdp.get_seq(GENOME_AC, START_BASE - 1, END_BASE)   # 0-based half-open
print('genomic ref (fwd):', g_seq)

# 2. Compute transcript ref and alt
from bioutils.sequences import reverse_complement
tx_ref = reverse_complement(g_seq)
tx_alt = reverse_complement(GENOMIC_ALT)
print('tx ref:', tx_ref)
print('tx alt:', tx_alt)

# 3. Trim common prefix/suffix
n_prefix = 0
while n_prefix < min(len(tx_ref), len(tx_alt)) and tx_ref[n_prefix] == tx_alt[n_prefix]:
    n_prefix += 1
n_suffix = 0
while (n_suffix < min(len(tx_ref) - n_prefix, len(tx_alt) - n_prefix)
       and tx_ref[-(n_suffix+1)] == tx_alt[-(n_suffix+1)]):
    n_suffix += 1
print('trimmed ref:', tx_ref[n_prefix: len(tx_ref)-n_suffix or None])
print('trimmed alt:', tx_alt[n_prefix: len(tx_alt)-n_suffix or None])
print('pos shift   : start +', n_prefix, '/ end -', n_suffix)
```

For a **plus-strand** transcript: `tx_ref = g_seq`, `tx_alt = GENOMIC_ALT` directly.

---

## Step 5 — Fetch transcript sequence to double-check ref

```python
from hgvs.alignmentmapper import AlignmentMapper
mapper = AlignmentMapper(hdp, 'NM_XXXXXX.X', 'NC_000011.10', 'splign')
cds_start_i = mapper.cds_start_i          # 0-based start of CDS in transcript

# c.527 → n-pos (0-based) = cds_start_i + 527 - 1
n_start = cds_start_i + C_START - 1
n_end   = cds_start_i + C_END
tx_seq  = hdp.get_seq('NM_XXXXXX.X', n_start, n_end)
print('transcript ref at c.{}_{}: {}'.format(C_START, C_END, tx_seq))
```

If `tx_seq` matches the trimmed `tx_ref` from Step 4, the expected edit is correct. If they disagree, the expected output in the plan is wrong.

---

## Decision tree

```
Step 0: any D/I gap within ≤10 bases of the variant?
  ├─ YES → variant may straddle or be adjacent to an internal gap;
  │         pay close attention to Steps 2-3 and the uncertain flag.
  └─ NO  → standard case; proceed to Step 1.

uncertain=False AND _variant_has_internal_gap=True?
  └─ YES → internal gap path runs (_get_altered_tx_sequence with original pos_g)

uncertain=False AND _variant_spans_i_segment=True?
  └─ YES → adjacent gap / double-gap cancellation path
            (_expand_pos_g_for_adjacent_gap + _get_altered_tx_sequence)
            See: hgvs-cigar-gap-mapping skill for full explanation

uncertain=False, no gap flags?
  └─ direct strand-flip path; result is mechanically correct.
     If plan says it's wrong, verify manually (Steps 4-5).

uncertain=True?
  └─ uncertain path runs; investigate why:
      ├─ Both endpoints in "I"?  → existing logic handles it
      ├─ One endpoint in "I", one in "="? → I-segment boundary bug
      └─ Endpoint in "D" or "N"? → intron/D-segment case
```

---

## Common pitfalls

| Pitfall | How to catch it |
|---------|-----------------|
| Internal D/I gap in the exon goes unnoticed | Step 0: always check exon CIGARs, not just N (introns) |
| Plan assumes `uncertain=True` but code gives `uncertain=False` | Step 2 |
| Expected c. output is a missense but genomic alt has net indel → always a frameshift | Steps 4-5 |
| Off-by-one: using `base - gc_offset` vs `base - 1 - gc_offset` | Step 3 comment |
| I-segment is ADJACENT to variant, not INSIDE it | Step 3: check if I-op is within `[start_offset, end_offset]` |
| Test expectation set from plan without independent verification | Always run Steps 1-5 before writing fix |
| `exon_info[i]` index trap in `_get_boundary()` | `exon_starts`/`exon_ends` are **sorted independently** of `exon_info`. Index `i` from the sorted arrays does NOT index back into `exon_info`. Additionally, the sentinel `exon_starts.append(exon_ends[-1])` makes `i` reachable at `len(exon_info)`. Fix: look up the exon by value — `next(e for e in exon_info if e["tx_start_i"] == exon_starts[i], None)` — and guard `if exon is not None`. |
| Roundtrip not possible around D/I positions| Don't expect this and don't try to correct it |
| NC_000011.10:g.119027721_119027726delinsTCACA is a difficult variant | The only true cdot is NM_001164277.1:c.532G>A |

---

## Related skills

- `hgvs-test-cache` — populate UTA cache after adding new test cases
- `debugging` — general debugging workflow
