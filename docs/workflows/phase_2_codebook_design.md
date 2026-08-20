# Phase 2: Codebook design

Phase 2 defines the codebook JSON mapping each target to the bit IDs that identify it in the SOLAR (splint/padlock, STARmap-style) assay.

For CLI-first workflows, generate the codebook with `mkprobes make-codebook`, then validate it with the checks below before running probe generation.

## Generating a codebook

```bash
mkprobes make-codebook <dataset> genes.tss.txt
```

This picks the smallest MHD code (minimum-Hamming-distance; matrices vendored in the package) whose capacity exceeds the target count by ≥5%, assigns each target three readout bits by seeded shuffle, fills spare capacity with `Blank-N` decoy codewords, and swaps any codeword that would perfectly confound imaging rounds onto a blank. Output defaults to `genes.codebook.json`; the codebook hash is logged for provenance.

**Optional — expression-informed assignment.** If per-target expression data is available (it is not required), the assignment can instead be optimized to balance total expression load across readout bits, so no bit's fluorescence is dominated by a few highly expressed genes:

```bash
# name of an annotation table registered at ingest, or a parquet/csv/tsv path
mkprobes make-codebook <dataset> genes.tss.txt --expression fpkm
mkprobes make-codebook <dataset> genes.tss.txt --expression expression.tsv --expression-column tpm
```

The table needs a `transcript_id` and/or `gene_id` column; targets missing from it are filled with the table median (with a warning). `--iterations` (default 200) controls how many assignments are tried; the one with the highest per-bit load entropy wins. For single-cell-derived optimization, the Python API offers `CodebookPickerSingleCell.find_optimalish` (percentile-based per-cell load balancing).

To extend an existing panel with non-overlapping bits:

```bash
mkprobes make-codebook <dataset> new_genes.txt --existing-codebook panel_v1.json
```

## Why this phase matters

The codebook is the contract between biological targets and encoding bits. Even if sequence design is perfect, a malformed or misaligned codebook causes panel-level decoding ambiguity. This phase prevents those failures before heavy compute begins.

## Inputs

- Validated gene list (`.txt`, one per line).
- Optional transcript-resolved target list.

## 2A. Validate genes/transcripts

### Check and normalize names (reference datasets)

```bash
mkprobes chkgenes data/mouse panel_a/genes.txt
```

Possible outputs:

- `panel_a/genes.converted.txt`
- `panel_a/genes.mapping.json`

**Intention:** ensure target names resolve to canonical identifiers before designing probes.

**Key arguments explained**

- `PATH`:
  - dataset root used for annotation lookup.
- `genes.txt`:
  - raw target list to validate; expected one target per line.

### Resolve transcripts

```bash
mkprobes convert-to-transcripts data/mouse panel_a/genes.converted.txt --mode canonical
```

Mode options include `canonical`, `gencode`, `ensembl`, `appris`, `apprisalt`, `longest`, `all`.

**Intention:** convert gene-centric lists into explicit transcript targets when required by your design policy.

**Key argument explained**

- `--mode`:
  - chooses transcript selection policy.
  - `canonical` is a practical default for stable panel design on reference (human/mouse) datasets.
  - transcript-selection policy should be fixed before phase 3 to avoid reproducibility drift.

### Custom datasets: offline transcript selection

On custom datasets (from `mkprobes ingest` or `mkprobes create-dataset`), transcript selection is fully offline:

- `--longest` picks the longest transcript sequence per gene; `--all` keeps every isoform (also available as `-m longest` / `-m all` for `convert-to-transcripts`).
- `--canonical` falls back to `--longest` on custom datasets.
- Input tokens can be transcript IDs (passed through unchanged), gene IDs/names, or symbols resolved case-insensitively through annotation tables registered on the dataset (e.g. an orthologs table).
- Unresolved tokens error out with close-match suggestions.

The Ensembl/mygene/APPRIS modes remain reference-only (human/mouse).

## 2B. Prepare codebook JSON

Create `panel_a/codebook.json` with shape:

```json
{
  "TargetA": [1, 9, 17],
  "TargetB": [2, 10, 18]
}
```

## Required schema

Use this schema for CLI probe-design workflows:

1. Top-level type: JSON object (`dict`).
2. Keys: target names (`str`) matching phase-3 `--gene` values.
3. Values: arrays of exactly 3 integers.
4. Each target's 3 integers must be distinct.
5. No duplicate 3-bit code tuples across targets.

Example valid entry:

```json
"Gad1": [2, 10, 18]
```

## Validation checklist

Run these checks before phase 3.

### Check JSON structure and duplicate code tuples

Save as `check_codebook_schema.py` and run with `python check_codebook_schema.py`:

```python
import json
from pathlib import Path

codebook = json.loads(Path("panel_a/codebook.json").read_text())
if not isinstance(codebook, dict) or not codebook:
    raise SystemExit("codebook.json must be a non-empty JSON object")

seen = set()
for target, bits in codebook.items():
    if not isinstance(target, str) or not target.strip():
        raise SystemExit(f"Invalid target key: {target!r}")
    if not isinstance(bits, list) or len(bits) != 3:
        raise SystemExit(f"{target}: expected exactly 3 bit IDs")
    if any((not isinstance(b, int)) for b in bits):
        raise SystemExit(f"{target}: all bit IDs must be integers")
    if len(set(bits)) != 3:
        raise SystemExit(f"{target}: bit IDs must be distinct within target")
    tup = tuple(sorted(bits))
    if tup in seen:
        raise SystemExit(f"Duplicate code tuple found: {tup}")
    seen.add(tup)

print(f"Codebook validation OK ({len(codebook)} targets)")
```

### Check target-name mismatch vs run list

Save as `check_codebook_targets.py` and run with `python check_codebook_targets.py`:

```python
import json
from pathlib import Path

targets = {x.strip() for x in Path("panel_a/genes.converted.txt").read_text().splitlines() if x.strip()}
codebook = set(json.loads(Path("panel_a/codebook.json").read_text()))

missing_in_codebook = sorted(targets - codebook)
extra_in_codebook = sorted(codebook - targets)
if missing_in_codebook:
    raise SystemExit(f"Missing targets in codebook: {missing_in_codebook}")
if extra_in_codebook:
    raise SystemExit(f"Extra targets in codebook: {extra_in_codebook}")
print("Target-name alignment OK")
```

### Practical acceptance criteria

You can proceed to phase 3 only when all checks pass:

1. schema check passes.
2. target-name alignment check passes.
3. codebook hash recorded in run metadata/logs.

## 2C. Track codebook identity

Compute a stable hash for run logs/provenance:

```bash
mkprobes hash panel_a/codebook.json
```

**Intention:** create a stable identifier so you can reliably match output files to the exact codebook used.

## HPC add-on guidance

- Freeze codebook JSON before compute fan-out.
- Keep codebook files under source control with dataset/version metadata.
- If using multiple panels, enforce explicit non-overlapping bit ranges in your panel design process.

## Failure modes

- Target names in codebook do not match validated gene/transcript names:
  - normalize target list first (`chkgenes`, `convert-to-transcripts`) and regenerate codebook.
- Reused/conflicting bits:
  - rerun your codebook generation process with explicit collision checks.
- Silent target-list drift after codebook creation:
  - regenerate and re-hash the codebook whenever target lists change.
