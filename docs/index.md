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
0. project    mkprobes init                  scaffold a project
1. dataset    mkprobes prepare | ingest      reference (mouse/human) | any species
2. targets    mkprobes chkgenes / convert-to-transcripts
3. codebook   mkprobes make-codebook
4. probes     mkprobes run-panel             candidates -> screen -> construct, all targets in parallel
5. panel QC   mkprobes filter-genes
6. assembly   mkprobes assemble              -> orderable oligos
```

## Where to start

**Designing a panel?** [Getting started](getting_started.md) is the
walkthrough — the whole workflow, start to finish. Everything else is a
companion to it.

**New to the assay?** Read
[What SOLAR is and what these probes do](what_is_solar.md) first: what these
probes physically are and why panels are encoded combinatorially.

**Setting up a machine?** [Before you start](before_you_start.md) covers the
external programs, the reference files to download, and how much disk, memory
and time each step needs.

**Not mouse or human?** Work through [Getting started](getting_started.md)
with [SOLAR probesets for a new species](workflows/solar_new_species.md)
open alongside it — that page covers only what differs.

**Looking up a flag?** [CLI reference](reference/cli.md) is generated from the
code, so it is always current.

```{toctree}
:maxdepth: 1
:caption: Start here

what_is_solar
before_you_start
installation
getting_started
troubleshooting
```

```{toctree}
:maxdepth: 1
:caption: How-to, step by step

workflows/build_a_dataset
workflows/choose_your_targets
workflows/design_the_codebook
workflows/design_probes
workflows/qc_your_panel
workflows/order_your_oligos
```

```{toctree}
:maxdepth: 1
:caption: Going further

workflows/solar_new_species
workflows/assay_and_under_the_hood
workflows/dataflow_map
```

```{toctree}
:maxdepth: 1
:caption: Reference

reference/cli
reference/columns
reference/file_formats
reference/fidelity
glossary
```
