---
name: hgvs-test-cache
description: How to populate the hgvs UTA test cache when adding new test cases that require UTA data not already in the local cache file.
---

# hgvs Test Cache Management

Tests in hgvs use a local cache (`tests/data/cache-py3.hdp`) so they can run offline without a live UTA database connection. When adding a new test case that requires data not already cached, you must populate the cache first.

## Cache Modes

| Mode string | Constant | Behavior |
|-------------|----------|----------|
| `'learn'`   | `LEARN=1` | Reads from cache; on miss, fetches from live UTA **and writes to cache** |
| `'run'`     | `RUN=2`   | Reads from cache only; raises `HGVSDataNotAvailableError` on miss |
| `'verify'`  | `VERIFY=3`| Always fetches live; raises if result differs from cache |

> **Common mistake**: Using `mode='write'` silently falls back to `mode=None` (no caching at all). Always use `'learn'` to populate the cache.

## Workflow: Adding a New Test Case

### 1. Identify what data is needed

Run the mapping with `mode='learn'` against a live UTA instance. This fetches missing data and writes it to the cache automatically:

```python
import sys
sys.path.insert(0, 'tests')
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser

UTA_URL = "postgresql://anonymous:anonymous@uta.biocommons.org:5432/uta/uta_20241220"
hdp = hgvs.dataproviders.uta.connect(UTA_URL, mode='learn', cache='tests/data/cache-py3.hdp')
am = hgvs.assemblymapper.AssemblyMapper(hdp)
hp = hgvs.parser.Parser()

var_g = hp.parse_hgvs_variant("NC_000002.12:g.73385901_73385903del")
var_c = am.g_to_c(var_g, "NM_015120.4")
print(str(var_c))  # capture this as your expected value
```

### 2. Verify the cache is sufficient

Re-run with `mode='run'` to confirm no live network calls are needed:

```python
hdp = hgvs.dataproviders.uta.connect(mode='run', cache='tests/data/cache-py3.hdp')
am = hgvs.assemblymapper.AssemblyMapper(hdp)
# run the same mapping — must succeed without errors
```

### 3. Write the test

```python
def test_projection_with_d_in_cigar(self):
    hgvs_g = "NC_000002.12:g.73385901_73385903del"
    hgvs_c = "NM_015120.4:c.34_36del"
    var_g = self.hp.parse_hgvs_variant(hgvs_g)
    var_c = self.am.g_to_c(var_g, "NM_015120.4")
    assert str(var_c) == hgvs_c
```

### 4. Run the test with cache

```bash
HGVS_CACHE_MODE=run uv run pytest tests/test_hgvs_assemblymapper.py -v
```

### 5. Commit the updated cache

The cache file `tests/data/cache-py3.hdp` is a binary shelve file checked into the repo. Always commit it alongside the new test:

```bash
git add tests/data/cache-py3.hdp tests/test_hgvs_assemblymapper.py
git commit -m "Add test for <description>; update UTA cache"
```

## Understanding CIGAR Alignment Strings

UTA stores transcript-to-genome alignments as CIGAR strings. When a test exercises a variant near an alignment discrepancy, the cache must include data for that region.

| Op | Meaning in hgvs (`CIGARMapper`) |
|----|----------------------------------|
| `=` | Match: transcript and genome agree |
| `X` | Mismatch: transcript and genome differ (1 base each) |
| `I` | Genome has an extra base **not present in the transcript** (advances genome/ref only) |
| `D` | Transcript has an extra base **not present in the genome** (advances transcript/tgt only) |
| `N` | Intron: genomic bases skipped in the transcript |

> **Convention note**: Raw exon CIGARs from `hdp.get_tx_exons()` (UTA) use the same I/D meanings as above after `build_tx_cigar` applies the transcript→genome swap. Always verify against `CIGARMapper.ref_pos` / `tgt_pos` rather than hand-counting.

Example — `NM_015120.4` exon 1 CIGAR `146=3D289=`:
- First 146 bases match
- `3D`: transcript has 3 extra bases with no genomic counterpart
- Next 289 bases match
- A deletion spanning across the `3D` boundary (e.g. `g.73385901_73385903del`) maps to `c.34_36del`

## Troubleshooting

**`HGVSDataNotAvailableError` during test run**: Cache is missing data. Re-run the mapping with `mode='learn'` to populate it.

**Test passes locally but fails in CI**: CI runs in `mode='run'`; the cache must be committed with all required entries.

**`mode='write'` doesn't persist data**: This mode string is not recognized — use `'learn'` instead.
