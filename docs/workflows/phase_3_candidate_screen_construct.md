# Phase 3: Candidate -> screen -> construct

Phase 3 builds actual probe candidates for each target and constructs the final SOLAR (splint/padlock, STARmap-style) sequence encodings.

This is the main `mkprobes` CLI execution stage.

## Why this phase matters

This phase is where design intent becomes concrete sequence output. It deliberately separates:

1. **candidate generation** (broad search),
2. **screening** (quality and off-target filtering),
3. **construction** (final encoded sequences).

Keeping these distinct makes debugging and parameter tuning tractable.

## Inputs

- Dataset path (`data/mouse`, `data/human`, or custom dataset).
- Codebook JSON (phase 2).
- Target list (`genes.converted.txt` or transcript list).

## Single-target debug path

Use this path first when tuning a new panel. Single-target runs make parameter effects obvious and reduce wasted compute.

### Step 1: candidates

```bash
mkprobes candidates data/mouse --gene Pcp4 --output panel_a/output
```

**What this step is doing**

- Enumerates candidate probe regions for one target transcript/gene.
- Computes off-target alignment context for downstream filtering.

Useful options:

- `--gene`:
  - target identifier to design against.
- `--output`:
  - directory where candidate and alignment files are written.
- `--allow gene1,gene2`:
  - allow specified related binders instead of auto-rejecting them.
- `--disallow gene1,gene2`:
  - explicitly reject specified binders.
- `--ignore-revcomp`:
  - disables reverse-complement matching; use only with a clear biological rationale.
- `--overwrite`:
  - forces regeneration when previous outputs exist.

Custom dataset notes (datasets from `ingest`/`create-dataset`):

- `--allow`/`--disallow` take transcript IDs (the FASTA record IDs), not gene names.
- Sibling isoforms of the target's gene are auto-allowed (derived from the dataset's GTF).
- The rRNA/tRNA blocklist is enforced automatically when the dataset carries one (`blocklist15.jf`).

### Step 2: screen

```bash
mkprobes screen panel_a/output Pcp4 --minimum 60 --maxoverlap 20 --restriction BamHI,KpnI --overwrite
```

**What this step is doing**

- Applies filtering logic to candidate probes.
- Enforces overlap/coverage strategy and optional restriction-site exclusions.

Key arguments explained:

- `DATA_DIR` and `GENE`:
  - identify which candidate files are screened.
- `--minimum`:
  - minimum probes to keep for this target; drives adaptive overlap behavior.
- `--maxoverlap`:
  - upper bound for overlap search when trying to satisfy `--minimum`.
- `--restriction`:
  - comma-separated enzymes whose recognition sites are filtered out.
- `--overwrite`:
  - recompute screening outputs even if previous outputs exist.

### Step 3: construct

```bash
mkprobes construct data/mouse panel_a/output --gene Pcp4 --codebook panel_a/codebook.json --target_probes 72 --restriction BamHI --restriction KpnI
```

**What this step is doing**

- Combines screened probes with codebook bits.
- Builds final encoded sequence payloads for downstream panel usage.

Key arguments explained:

- first positional `PATH`:
  - dataset root used for reference context.
- second positional `OUTPUT_PATH`:
  - location of screened/candidate files and final outputs.
- `--gene`:
  - target to construct.
- `--codebook`:
  - target-to-bit mapping JSON.
- `--target_probes`:
  - maximum probes retained/constructed per target.
- `--restriction`:
  - restriction-site policy propagated to final selection/filenames.

## Panel batch path (recommended)

```bash
mkprobes run-panel data/mouse panel_a/codebook.json
```

`run-panel` applies the command trio to every target in the codebook across
parallel workers (16 by default), with production settings
(`--minimum 60 --maxoverlap 0 --restriction BamHI,KpnI --target-probes 48`,
all overridable). Behavior:

- finished targets (final parquet present) are skipped; `--overwrite` redoes them;
- a `<codebook>.acceptable.json` allow-list (or `--allow-file`) is passed to
  `candidates --allow` and forces a re-screen/re-construct for those genes;
- failures are isolated per gene, logged to `output/<gene>.log`, collected in
  `<codebook>.failed.txt`, and the command exits non-zero;
- `mkprobes run-panel <dataset> <codebook> <gene>` re-runs one target (forced);
- `--list-failed` / `--list-failed-all` triage targets without final outputs.

**Intention:** apply exactly the same command trio to each target for
consistent panel assembly and easier HPC execution — one command per panel.

## Output files

Per gene/transcript, common outputs in `output/`:

- `<target>_all.parquet`
- `<target>_bowtie.parquet`
- `<target>_crawled.parquet`
- `<target>_crawled.stats.json`
- `<target>_screened_ol*.parquet`
- `<target>_screened_ol*.stats.json`
- `<target>_final_<restriction>_<bits>.parquet`

## Quality controls

- Validate final counts against intended minimum.
- Inspect `_crawled.stats.json` for dominant off-target binders.
- Re-run selected targets with `--overwrite` after tuning parameters.

Suggested command-level QC loop:

1. run one target with conservative defaults;
2. inspect outputs/stats;
3. lock parameter policy;
4. scale to full batch.

## HPC add-on guidance

- Shard target lists and run multiple batch jobs.
- Keep `output/` on local scratch while running, then archive final parquet and logs.
- Match job CPU allocation to per-job command parallelism.

## Failure modes

- No probes after match/off-target filtering:
  - verify transcript selection and allow-list strategy.
- Existing outputs silently skipped:
  - rerun with `--overwrite`.
- Sparse final probe counts:
  - inspect overlap strategy and accepted off-target controls.
- Inconsistent behavior across targets:
  - verify a single shared codebook and consistent command arguments across the full loop.
