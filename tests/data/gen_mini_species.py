"""
Generates the committed mini-species fixtures in this directory.

One synthetic genome (3 scaffolds) + the same 7-transcript/6-gene structure
annotated in four de novo styles (AUGUSTUS, StringTie, MAKER, NCBI-like).
Deterministic (seeded); outputs are committed, this script documents their
provenance and regenerates them if ever needed:

    uv run python tests/data/gen_mini_species.py
"""

import random
from pathlib import Path

HERE = Path(__file__).parent
rng = random.Random(42)

SCAFFOLDS = {"scaffold_1": 20_000, "scaffold_2": 15_000, "scaffold_3": 8_000}

# (gene, transcript, scaffold, strand, exons, biotype)
MODELS: list[tuple[str, str, str, str, list[tuple[int, int]], str]] = [
    ("g1", "g1.t1", "scaffold_1", "+", [(1001, 1800), (2501, 3400), (4101, 5100)], "protein_coding"),
    ("g1", "g1.t2", "scaffold_1", "+", [(1001, 1800), (4101, 5100)], "protein_coding"),
    ("g2", "g2.t1", "scaffold_1", "-", [(8001, 9200), (10001, 11000)], "protein_coding"),
    ("g3", "g3.t1", "scaffold_2", "+", [(501, 2500)], "protein_coding"),
    ("g4", "g4.t1", "scaffold_2", "-", [(5001, 5600), (6301, 6900), (7601, 8200), (9001, 9800)], "protein_coding"),
    ("g5", "g5.t1", "scaffold_3", "+", [(1001, 1700), (2501, 3300)], "protein_coding"),
    ("r1", "r1.t1", "scaffold_3", "+", [(5001, 5120)], "rRNA"),
]

# Same structures, different naming conventions per annotator style.
# Unique numeric index per gene (avoids ID collisions between g1 and r1).
GENE_IDX = {"g1": 1, "g2": 2, "g3": 3, "g4": 4, "g5": 5, "r1": 6}

STYLES = {
    "augustus": {
        "source": "AUGUSTUS",
        "gene_id": lambda g: g,
        "tx_id": lambda g, t: t,
        "extra": lambda g, t, b: "",
    },
    "stringtie": {
        "source": "StringTie",
        "gene_id": lambda g: f"STRG.{GENE_IDX[g]}",
        "tx_id": lambda g, t: f"STRG.{GENE_IDX[g]}.{t.split('.t')[1]}",
        "extra": lambda g, t, b: f' ref_gene_id "REF_{g.upper()}"; FPKM "12.5";',
    },
    "maker": {
        "source": "maker",
        "gene_id": lambda g: f"maker-scaffold1-snap-gene-0.{GENE_IDX[g]}",
        "tx_id": lambda g, t: f"maker-scaffold1-snap-gene-0.{GENE_IDX[g]}-mRNA-{t.split('.t')[1]}",
        "extra": lambda g, t, b: "",
    },
    "ncbi_style": {
        "source": "Gnomon",
        "gene_id": lambda g: f"gene-LOC10{GENE_IDX[g]}",
        "tx_id": lambda g, t: f"rna-XM_00{GENE_IDX[g]}00{t.split('.t')[1]}.1",
        "extra": lambda g, t, b: f' gene "SYM{g.upper()}"; gene_biotype "{b}";',
    },
}


def main() -> None:
    seqs = {
        name: "".join(rng.choice("ACGT") for _ in range(length)) for name, length in SCAFFOLDS.items()
    }
    with open(HERE / "mini_genome.fa", "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")

    for style, cfg in STYLES.items():
        lines: list[str] = ["#!genome-build mini_species_v1"]
        for gene, tx, scaffold, strand, exons, biotype in MODELS:
            gid, tid = cfg["gene_id"](gene), cfg["tx_id"](gene, tx)
            extra = cfg["extra"](gene, tx, biotype)
            start, end = exons[0][0], exons[-1][1]
            attrs = f'gene_id "{gid}"; transcript_id "{tid}";{extra}'
            lines.append(
                f"{scaffold}\t{cfg['source']}\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
            )
            for i, (es, ee) in enumerate(exons, 1):
                lines.append(
                    f"{scaffold}\t{cfg['source']}\texon\t{es}\t{ee}\t.\t{strand}\t.\t{attrs} "
                    f'exon_number "{i}";'
                )
        (HERE / f"{style}.gtf").write_text("\n".join(lines) + "\n")

    print(f"Wrote mini_genome.fa ({sum(SCAFFOLDS.values())} bp) and {len(STYLES)} GTF styles to {HERE}")


if __name__ == "__main__":
    main()
