# Getting Started

This guide gives a practical first pass through `mkprobes`, focused on CLI-first design of SOLAR (splint/padlock, STARmap-style) probesets.

## What you will do

1. Prepare a probe-design dataset.
2. Validate your target list and transcript mapping.
3. Create a codebook JSON (input file).
4. Generate candidate/screened/final probes with `mkprobes` commands.
5. Run panel-level QC on probe counts.

## Minimal first-run checklist

- [ ] Installation complete ({doc}`installation`); `mkprobes --help` works.
- [ ] You have `bowtie2`, `jellyfish`, and `gffread` in `PATH`.
- [ ] You have a writable project folder with enough disk for intermediate parquet outputs.

## Suggested project layout

```text
project/
├── data/
│   └── mouse/                      # created by mkprobes prepare (reference mode)
├── panel_a/
│   ├── genes.txt
│   ├── genes.converted.txt         # optional from chkgenes
│   ├── codebook.json               # mapping: gene/transcript -> bit IDs
│   └── output/                     # candidate/screen/final parquet outputs
```

## Quick start: toy CLI probe-design run

The sequence below is intentionally ordered so each step validates prerequisites for the next one.

### 1. Create and prepare a working area

```bash
mkdir -p project/panel_a/output
cd project
```

Reference dataset (mouse/human):

```bash
mkprobes prepare data --species mouse --threads 16
```

This creates/updates `data/mouse` with required GTF/FASTA, bowtie index, and k-mer files.

For any other species, build a dataset from a genome FASTA + GTF/GFF3 with `mkprobes ingest` — see the runbook {doc}`workflows/solar_new_species`.

Why this first: every downstream `mkprobes` command depends on this dataset being complete and consistent.

### 2. Create a toy gene list

Create `panel_a/genes.txt` with one target per line, for example:

```text
Pcp4
Gad1
Slc17a7
Camk2a
Rbfox3
```

Validate names and emit converted names if needed:

```bash
mkprobes chkgenes data/mouse panel_a/genes.txt
```

If generated, prefer `panel_a/genes.converted.txt` as canonical target input.

Why this step: target-name normalization prevents silent mismatches between biological intent and command inputs.

### 3. Create a toy codebook JSON

Prepare `panel_a/codebook.json`:

```json
{
  "Camk2a": [1, 9, 17],
  "Gad1": [2, 10, 18],
  "Pcp4": [3, 11, 19],
  "Rbfox3": [4, 12, 20],
  "Slc17a7": [5, 13, 21]
}
```

Record a stable identifier for this codebook:

```bash
mkprobes hash panel_a/codebook.json
```

Run strict codebook validation before generation. Save as `check_codebook.py` and run with `python check_codebook.py`:

```python
import json
from pathlib import Path

codebook = json.loads(Path("panel_a/codebook.json").read_text())
targets = {x.strip() for x in Path("panel_a/genes.converted.txt").read_text().splitlines() if x.strip()}

if set(codebook) != targets:
    raise SystemExit(f"Target mismatch. codebook_only={sorted(set(codebook)-targets)} targets_only={sorted(targets-set(codebook))}")
for t, bits in codebook.items():
    if not isinstance(bits, list) or len(bits) != 3 or any((not isinstance(b, int)) for b in bits):
        raise SystemExit(f"{t}: invalid bits (expected list of 3 ints)")
    if len(set(bits)) != 3:
        raise SystemExit(f"{t}: duplicate bit IDs within target")
print("Codebook validation OK")
```

For the full validation checklist, see {doc}`workflows/phase_2_codebook_design`.

Why this step: the codebook is the panel contract; validating now avoids wasting compute in phase 3.

### 4. Run candidate -> screen -> construct with CLI

Run per target:

```bash
while read -r gene; do
  mkprobes candidates data/mouse --gene "$gene" --output panel_a/output
  mkprobes screen panel_a/output "$gene" --minimum 60 --maxoverlap 20 --restriction BamHI,KpnI
  mkprobes construct data/mouse panel_a/output --gene "$gene" --codebook panel_a/codebook.json --target_probes 72 --restriction BamHI --restriction KpnI
done < panel_a/genes.converted.txt
```

Why this exact order:

1. `candidates` explores possible probes.
2. `screen` removes poor/off-target options.
3. `construct` emits final encoded outputs tied to your codebook.

### 5. Panel-level QC

Filter for targets meeting minimum probe count:

```bash
mkprobes filter-genes panel_a/output --genes panel_a/genes.converted.txt --min-probes 48 --out panel_a/genes.pass.txt
```

Why this step: panel quality should be decided at the panel level, not per-target in isolation.

## Expected intermediate outputs

For each gene/transcript, expect files like:

- `*_all.parquet`
- `*_bowtie.parquet`
- `*_crawled.parquet`
- `*_screened_ol*.parquet`
- `*_final_BamHIKpnI_<bits>.parquet`

## HPC quick adaptation

- Keep `data/mouse` (or `data/human`) on shared read-mostly storage.
- Write per-job `panel/output` to node-local scratch.
- Copy back final parquet outputs and QC summaries.
- Run one-gene smoke tests before full panel fan-out.

## Next steps

- New-species end-to-end runbook: {doc}`workflows/solar_new_species`
- Assay context + internal algorithm view: {doc}`workflows/assay_and_under_the_hood`
- Input/output flow map: {doc}`workflows/dataflow_map`
- Command details: {doc}`reference/cli`
