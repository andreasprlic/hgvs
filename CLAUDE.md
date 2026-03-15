# Claude Agentic Framework — HGVS Project

**Project**: `hgvs` — Python library for parsing, formatting, validating, normalizing, and mapping genomic sequence variants per [HGVS nomenclature](https://hgvs-nomenclature.org/).
**Package manager**: `uv` | **Python**: 3.11–3.13 | **License**: Apache-2.0
**Repo**: https://github.com/biocommons/hgvs/ | **Docs**: https://hgvs.readthedocs.io/

## Quick Reference

```bash
# Setup (one-time)
make devready                   # create venv, install deps, pre-commit hooks
source .venv/bin/activate

# Quality gates — run before every commit
make cqa                        # lock check, lint (ruff), reformat, dependency check

# Testing
make test                       # pytest with coverage (VCR cassette cache, offline-safe)
make test-learn                 # add new network responses to VCR test cache
make test-relearn               # destroy and rebuild test cache
pytest -m quick                 # run only quick tests
pytest -m "not network"         # skip tests that need a live UTA/SeqRepo connection
pytest tests/issues/            # regression tests for specific reported issues

# Docs
make docs-serve                 # build and serve docs locally (MkDocs + Material theme)
make docs-test                  # test docs build with zero warnings

# Build & publish
make build                      # build distribution packages
make publish                    # publish to PyPI (requires UV_PUBLISH_TOKEN env var)

# CLI
hgvs-shell                      # interactive HGVS Python shell
```

## Project Overview

`hgvs` is a production-grade bioinformatics library used in clinical genomics. It handles:

- **Parsing**: PEG grammar-based parsing of HGVS variant strings (e.g. `NM_000518.4:c.315+1G>A`)
- **Formatting**: Convert internal objects back to standard HGVS strings
- **Validation**: Verify variants against reference sequences (requires SeqRepo)
- **Normalization**: Rewrite variants in canonical form
- **Mapping**: Project variants between genome (`g.`), transcript (`c.`, `n.`, `r.`), and protein (`p.`) coordinate systems

**External data dependencies** (required for non-cached tests):
- **UTA** (Universal Transcript Archive): PostgreSQL DB — default `postgresql://anonymous@uta.biocommons.org:5432/uta/uta_20241220`; set via `UTA_DB_URL` env var
- **SeqRepo**: Local or remote sequence repository; set via `HGVS_SEQREPO_DIR` env var
- Docker Compose setup in `./misc/docker-compose.yml` for full local test environment

## Tech Stack

Defined in `.claude/rules/tech-strategy.md` — auto-loaded for every session.

## Source Layout

```
src/hgvs/
  parser.py           # HGVS string parsing (PEG grammar via parsley)
  sequencevariant.py  # Core SequenceVariant class
  variantmapper.py    # g/c/n/r/p coordinate mapping
  alignmentmapper.py  # sequence alignment handling
  assemblymapper.py   # cross-assembly mapping
  validator.py        # variant validation
  normalizer.py       # canonical form normalization
  projector.py        # variant projection between sequences
  easy.py             # high-level convenience API
  config.py           # global configuration
  dataproviders/      # pluggable UTA, NCBI, SeqFetcher backends
  extras/babelfish.py # format conversion utilities
  utils/              # CIGAR mapping, position helpers, etc.
  generated/          # auto-generated grammar parser (do not edit)
  pretty/             # pretty-printing utilities

tests/
  test_hgvs_parser.py
  test_hgvs_variantmapper*.py
  test_hgvs_alignmentmapper.py
  test_hgvs_assemblymapper.py
  test_hgvs_dataproviders_uta.py
  issues/             # regression tests for GitHub issues
  cassettes/          # VCR HTTP cassettes for reproducible offline tests
```

## Core Principles

These seven principles distill every rule, skill, and standard in this framework. Follow them and everything else follows.

### 1. Understand First
Read before writing; grep before creating; verify APIs via docs before assuming training data is current.

### 2. Prove It Works
Write tests first, run quality gates (tests, linter, types, build) before every commit, and add a regression test for every bug fix.

### 3. Keep It Safe
No secrets in code, validate all input, use parameterized queries, apply least privilege, and flag vulnerabilities immediately.

### 4. Keep It Simple
Single responsibility, no premature abstraction, delete dead code, avoid `any` types, fix warnings before committing.

### 5. Don't Repeat Yourself
Check `.claude/skills/` before generating ad-hoc solutions; maintain a single source of truth for business logic.

### 6. Ship It
Work on a branch, commit iteratively, and push to remote — work isn't done until `git push` succeeds.

### 7. Leave a Trail
Artifacts in `./artifacts/`, track work with Beads (`bd` CLI), document decisions in ADRs, name things clearly.

Full details in `.claude/rules/` (auto-loaded).

## Workflow

**Branching**: Always branch from `main`. Never commit directly to `main`.

**Planning flow**: PR-FAQ → PRD → ADR → Design Spec → Plan → Implementation Beads

**Artifacts**: All planning docs stored in `./artifacts/`:

| Type | Pattern | Example |
|------|---------|---------|
| Vision | `pr_faq_[feature].md` | `pr_faq_user_auth.md` |
| Requirements | `prd_[feature].md` | `prd_user_auth.md` |
| Architecture | `adr_[topic].md` | `adr_database_choice.md` |
| System Design | `system_design_[component].md` | `system_design_api.md` |
| Design | `design_spec_[component].md` | `design_spec_login_form.md` |
| Roadmap | `roadmap_[project].md` | `roadmap_mvp.md` |
| Plan | `plan_[task].md` | `plan_api_refactor.md` |
| Security Audit | `security_audit_[date].md` | `security_audit_2025-01.md` |
| Post-Mortem | `postmortem_[incident-id].md` | `postmortem_inc-2025-001.md` |

## Working Directories

| Directory | Purpose | Lifecycle |
|-----------|---------|-----------|
| `./artifacts/` | Durable documents (plans, ADRs, PRDs, design specs) | Committed to repo |
| `./scratchpad/` | Ephemeral working notes, exploration output, draft content | Gitignored, disposable |


## MCP Tools

| Tool | Use For |
|------|---------|
| Sequential Thinking | Complex analysis, trade-off evaluation |
| Context7 | Library documentation lookup |
| Filesystem | File system operations beyond workspace |

## Skills

Check `.claude/skills/` before ad-hoc generation. Skills are auto-suggested based on context via `.claude/skills/skill-rules.json`.
