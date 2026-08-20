# mkprobes

`mkprobes` designs **SOLAR** probesets — the lab's splint/padlock,
STARmap-style combinatorial FISH assay — from a species' transcriptome. It
covers the full workflow: reference/dataset preparation (including genome +
annotation ingestion for non-traditional model species), target selection,
per-gene candidate generation and off-target screening, probe construction
against a codebook, panel QC, and final assembly into an orderable oligo
pool.

The workflow at a glance:

```text
1. dataset    mkprobes prepare | ingest      reference (mouse/human) | any species
2. targets    mkprobes chkgenes / convert-to-transcripts
3. codebook   mkprobes make-codebook
4. probes     mkprobes candidates -> screen -> construct   (per target)
5. panel QC   mkprobes filter-genes
6. assembly   scripts/probegen/2_assemble_manifest.py      -> orderable oligos
```

New species? Start with the flagship runbook:
[SOLAR probesets for a new species](workflows/solar_new_species.md).

```{toctree}
:maxdepth: 1
:caption: Start here

installation
getting_started
troubleshooting
```

```{toctree}
:maxdepth: 1
:caption: Workflows

workflows/solar_new_species
workflows/phase_1_dataset_prep
workflows/phase_2_codebook_design
workflows/phase_3_candidate_screen_construct
workflows/phase_4_panel_qc_export
workflows/phase_5_manifest_assembly
workflows/assay_and_under_the_hood
workflows/dataflow_map
```

```{toctree}
:maxdepth: 1
:caption: Reference

reference/cli
reference/file_formats
```
