# Design the codebook

**Step 3 of the {doc}`../getting_started` workflow.**

The codebook assigns each target three readout bits — the combination that
identifies it during imaging. It is the contract between your biological
targets and the encoding, and everything downstream is built against it. A
codebook that is wrong produces probes that are individually perfect and
collectively undecodable, so it is worth getting right before the expensive
step.

If the idea of combinatorial bits is new, {doc}`../what_is_solar` explains
what the bits physically are and why panels are encoded this way.

## Generate it

```bash
mkprobes make-codebook data/mouse panel_a/genes.converted.tss.txt -o panel_a/codebook.json
```

**Pass `-o` explicitly.** Without it, the output name is derived from the
input list rather than being a fixed default: `genes.converted.tss.txt`
produces `genes.converted.tss.codebook.json`. That is rarely what anyone
expects, and the manifest you generated with `mkprobes init` refers to
`codebook.json`.

What the command does:

1. Picks the smallest vendored MHD (minimum-Hamming-distance) code whose
   capacity exceeds your target count by at least 5%.
2. Assigns each target three readout bits.
3. Fills the spare capacity with `Blank-N` decoy codewords — these are how you
   measure the false-positive rate on real data.
4. Swaps any codeword that would perfectly confound imaging rounds onto a
   blank.
5. Logs a hash of the finished codebook.

`--n-bits` overrides the automatic sizing; you rarely want to.

## Optional: balance the readout load

If you have per-target expression data, the assignment can be optimized so no
single readout bit is dominated by a handful of very highly expressed genes:

```bash
# a table registered on the dataset...
mkprobes make-codebook data/mouse panel_a/genes.converted.tss.txt -o panel_a/codebook.json \
    --expression fpkm

# ...or any file on disk
mkprobes make-codebook data/mouse panel_a/genes.converted.tss.txt -o panel_a/codebook.json \
    --expression expression.tsv --expression-column tpm
```

This is genuinely optional — omit it and you get a plain seeded assignment
(`--seed`), which is a perfectly good panel.

The table needs a `transcript_id` and/or `gene_id` column. Targets missing
from it are filled with the table median, with a warning. `--iterations`
(default 200) sets how many assignments are tried; the one with the most even
per-bit load wins. `--expression-column` is only needed when the value column
is ambiguous.

For single-cell-derived optimization, the Python API offers
`CodebookPickerSingleCell.find_optimalish`, which balances per-cell load by
percentile.

## Extending an existing panel

To add targets to a panel you have already ordered, without reusing bits:

```bash
mkprobes make-codebook data/mouse new_genes.tss.txt -o panel_b/codebook.json \
    --existing-codebook panel_a/codebook.json
```

This derives the bit offset from the old codebook and refuses gene or bit
overlap rather than silently colliding. `--offset` sets the offset by hand and
is mutually exclusive with `--existing-codebook`.

## Which codebook produced a given file

Every codebook has a short hash — a stable identifier for that exact set of
targets and bit assignments. You do not have to record it anywhere: the tool
does.

`make-codebook` writes it beside the codebook, so it outlives the terminal:

```text
panel_a/codebook.json
panel_a/codebook.hash     <- 25dd20
```

It is also stamped into the provenance of every file designed against that
codebook — each per-target `_final_*.parquet`, and the assembled pool. So the
question "which codebook produced this?" is answered by the file itself:

```bash
mkprobes provenance panel_a/output/Sox2-201_final_BamHIKpnI_2,10,18.parquet
```

```text
{
  "codebook_hash": "25dd20",
  "bits": [2, 10, 18],
  "stage": "construct",
  ...
}
```

Match that against `codebook.hash` to confirm an output came from the codebook
you think it did. This matters when a codebook is regenerated: a different seed
or an edited target list produces different bit assignments, and the outputs
are otherwise indistinguishable — same target, same file name, different panel.

To print the hash of any codebook directly:

```bash
mkprobes hash panel_a/codebook.json
```

:::{note}
The hash covers the codebook as written, Blank codes included. Re-serialising
it — different indentation, different key order — does not change it, because
hashing sorts the keys first. Changing any target's bits does.
:::

## What a codebook looks like

```json
{
  "Sox2-201": [2, 10, 18],
  "Pax6-201": [1, 9, 17],
  "Blank-1":  [3, 11, 19]
}
```

The rules `mkprobes` enforces, and which you should preserve if you ever edit
one by hand:

1. A JSON object at the top level.
2. Keys are target names matching your target list exactly.
3. Values are arrays of exactly three integers.
4. The three integers within a target are distinct.
5. No two targets share the same three-bit code.

Generally: do not edit it by hand. Regenerate it. The file is cheap to
produce and hand edits are how rules 4 and 5 get broken.

## Before you spend hours on probes

Check that the codebook and the target list still agree — a target list edited
after the codebook was generated is a common and expensive mistake:

```bash
python - <<'EOF'
import json
from pathlib import Path

targets = {
    line.split("#", 1)[0].strip()
    for line in Path("panel_a/genes.converted.tss.txt").read_text().splitlines()
    if line.split("#", 1)[0].strip()
}
codebook = {k for k in json.loads(Path("panel_a/codebook.json").read_text()) if not k.startswith("Blank")}

if missing := sorted(targets - codebook):
    raise SystemExit(f"Missing from codebook: {missing}")
if extra := sorted(codebook - targets):
    raise SystemExit(f"In codebook but not in target list: {extra}")
print(f"OK - {len(targets)} targets aligned")
EOF
```

If that fails, regenerate the codebook from the current target list rather
than patching either file.

## When it goes wrong

- **Codebook names do not match the target list** — regenerate from
  `genes.converted.tss.txt`, do not rename keys.
- **You changed the target list after generating the codebook** — regenerate
  and re-hash. Any probes already designed against the old codebook are still
  valid for the targets that kept their bits, but `run-panel` works from the
  codebook, so the two must agree.
- **Two panels imaged together reuse bits** — generate the second with
  `--existing-codebook` pointing at the first.

---

Next: {doc}`design_probes`.
