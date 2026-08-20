# Phase 5: Manifest assembly and export

Phase 5 aggregates per-target final SOLAR (splint/padlock, STARmap-style) probe files and emits deliverables for ordering/downstream use.

Assembly is driven by `scripts/probegen/2_assemble_manifest.py`, run from the repository checkout.

## Inputs

- `manifest.json` (list of `ProbeSet` entries).
- Completed `output/*_final_*.parquet` files for targets in each codebook.

## Commands

### Check short/underperforming targets

```bash
uv run python scripts/probegen/2_assemble_manifest.py panel_a/manifest.json short 12
```

This reports targets with fewer than threshold probes and may trigger interactive accept-list flow when off-target tables are available.

### Generate final assembled outputs

```bash
uv run python scripts/probegen/2_assemble_manifest.py panel_a/manifest.json gen
```

The header/footer table and MHD matrices are vendored inside the package (loaded via `importlib.resources`); pass `--headerfooter` to override the header/footer table.

For non-model organisms, `gen` accepts:

- `--rm-species '<taxon>'`: passed verbatim as RepeatMasker `-species` (any taxon RepeatMasker's library supports).
- `--skip-repeatmasker`: skip repeat masking (and silence the skip warning).

## What assembly does

For each manifest probeset:

1. Loads codebook and per-target final parquet files.
2. Selects top probes per target (count depends on `n_probes` and species logic).
3. Optionally runs RepeatMasker (built-in mouse/human mapping, or `--rm-species`).
4. Builds final splint/padlock constructs.
5. Writes parquet + FASTA outputs.
6. Emits `<panel>_final.txt` and cumulative `_allout<timestamp>.txt`.
7. Writes `<name>.provenance.json` per panel (timestamp, mkprobes version, codebook hash, probeset config, RepeatMasker status).

## Typical outputs

Under `<panel_root>/generated/`:

- `<panel>.parquet`
- `<panel>_pad.fasta`
- `<panel>_splint.fasta`
- `<panel>_final.txt`
- `<panel>.provenance.json`
- `_allout<ISO_TIMESTAMP>.txt`

## Manual acceptance workflow

When `short` finds low-count targets and off-target files exist, the script can prompt for acceptable off-target transcripts using `questionary`.

Accepted targets are stored in:

- `<codebook_stem>.acceptable.json`

Then rerun phase 3 for affected targets.

## HPC add-on guidance

- Run `short` first as a lightweight QC gate.
- Run `gen` after freezing all acceptable overrides.
- Archive manifest, codebook, acceptable JSON, and generated outputs together as a single release bundle.

## Failure modes

- Missing final parquet files:
  - ensure phase 3 completed for every target in codebook.
- RepeatMasker unavailable:
  - install it (optional dependency; only needed here) or run with `--skip-repeatmasker`.
- Mismatch in manifest paths:
  - validate that `codebook` paths in manifest are relative to manifest location as expected by scripts.
