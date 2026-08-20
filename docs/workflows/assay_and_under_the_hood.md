# Assay context and under-the-hood logic

This page explains two things:

1. Why this probe-design workflow supports SOLAR (splint/padlock, STARmap-style), a combinatorial spatial transcriptomics assay.
2. What the `mkprobes` pipeline is doing internally at each stage.

The descriptions below are based on repository implementation in `src/mkprobes/*`.

For a step-by-step execution view, see {doc}`dataflow_map`.

## Spatial transcriptomics context

The design workflow produces transcript-targeting probe sequences that can be decoded in tissue-space rather than bulk lysate.

In practical terms:

1. Each target transcript is assigned a code (bit IDs) in a codebook.
2. Multiple probes are designed per target to increase detection robustness.
3. Probe sequences are filtered for specificity and thermodynamic feasibility.
4. Final constructs include encoded readout content tied to target identity.

Why this matters for a spatial assay:

- **Identity:** bit assignments separate transcripts through combinatorial coding.
- **Specificity:** off-target and sequence-quality filters reduce false signal.
- **Robustness:** multi-probe-per-target design mitigates local dropout.
- **Compatibility:** restriction/homopolymer checks reduce downstream chemistry risk.

## End-to-end internal flow

## 1. Candidate enumeration (`mkprobes candidates`)

Implementation anchors:

- `src/mkprobes/candidates.py`
- `src/mkprobes/utils/_crawler.py`
- `src/mkprobes/starmap/starmap.py`

Internal logic:

1. Fetch target transcript sequence from dataset reference.
2. Slide windows and grow candidate subsequences with constraints:
   - length bounds (default around 43–55 nt)
   - GC bounds
   - Tm bounds
   - hairpin threshold
3. Reject homopolymer-heavy sequences.
4. Split candidate into assay-facing components (`splint`/`padlock` style split) using target Tm logic.
5. Align generated sequences against reference index (Bowtie-based path) to estimate off-target behavior.
6. Compute per-sequence quality features (GC, hairpin, Tm, homopolymer-derived flags, etc.).

Design intention:

- Generate a wide but constrained search space, then annotate each candidate with enough metrics for strict downstream filtering.

## 2. Specificity and quality scoring

Implementation anchors:

- `SAMFrame.filter_by_match(...)` usage in `candidates.py`
- `agg_tm_offtarget(...)` usage in `candidates.py`
- `PROBE_CRITERIA` in `src/mkprobes/utils/_filtration.py`

Internal logic:

1. Keep candidates that match acceptable transcript set with minimum match criteria.
2. Aggregate off-target thermodynamic summaries (`max_tm_offtarget` style metric).
3. Apply sequence-rule criteria (quadruplet/homopolymer/composition checks).
4. Compute an aggregate quality score (`oks`) from boolean criteria.

Design intention:

- Convert alignment and sequence features into a consistent quality ranking before overlap optimization.

## 3. Screening and overlap optimization (`mkprobes screen`)

Implementation anchors:

- `src/mkprobes/screen.py`
- `src/mkprobes/utils/_filtration.py`

Internal logic:

1. Optionally remove candidates containing specified restriction enzyme sites.
2. Apply tiered filtering criteria (strict to progressively relaxed thresholds).
3. Use overlap-aware selection (`find_overlap` / weighted variant) to choose a non-redundant probe set with coverage.
4. Reshape selected split components into paired records and write screened parquet outputs.

Design intention:

- Balance specificity and physical transcript coverage while capping redundancy.

## 4. Construction and encoding (`mkprobes construct`)

Implementation anchors:

- `src/mkprobes/codebook/finalconstruct.py`

Internal logic:

1. Load screened probes and target bit IDs from codebook.
2. Map bit IDs to readout sequences (`readout_ref_filtered.csv` lookup).
3. Stitch readout content with padlock payload.
4. Skip constructions that violate hard constraints (e.g., homopolymers, restricted motifs like BamHI site in stitching path).
5. Emit final per-target parquet with encoded sequence columns.

Design intention:

- Produce assay-ready encoded probe constructs tied deterministically to the codebook.

## 5. Panel QC (`mkprobes filter-genes`)

Implementation anchors:

- `filter_genes` in `src/mkprobes/codebook/finalconstruct.py`

Internal logic:

1. Count available screened probes per target.
2. Flag targets below threshold.
3. Optionally output a pass-list for downstream panel release.

Design intention:

- Enforce a panel-level minimum viability criterion before assay execution.

## Why these choices help spatial assays specifically

1. **High-specificity targeting in-place:** alignment + off-target thermodynamic checks lower cross-binding risk in tissue.
2. **Combinatorial identity encoding:** bit-to-target mapping enables multiplexing without one-channel-per-target scaling.
3. **Coverage-aware selection:** overlap optimization spreads probes along transcript regions, improving detection stability.
4. **Chemistry-aware constraints:** homopolymer/restriction filtering reduces assay-fragile constructs.

## Practical interpretation notes

1. This pipeline is designed as a conservative filter stack: broad generation followed by increasingly strict biochemical and specificity constraints.
2. Codebook integrity is as important as sequence quality; treat codebook validation as a hard gate.
3. Panel quality should be evaluated globally (`filter-genes`), not from single-target success alone.
