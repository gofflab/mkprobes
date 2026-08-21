"""
Genome + annotation ingestion for SOLAR probe design datasets.

`mkprobes ingest` turns a genome FASTA + GTF (typically a de novo annotation:
MAKER, AUGUSTUS/BRAKER, StringTie, NCBI Gnomon) into a ready-to-use probe
design dataset:

1. Check external tools (gffread, bowtie2-build, jellyfish).
2. Normalize the annotation (GFF3 -> GTF via `gffread -T`; decompress .gz).
3. Validate the annotation strictly against the genome and emit a report
   (`validation_report.json`) - the report is the product as much as the
   dataset: every failure names what was checked, the offending values, and
   the fix.
4. Extract transcript sequences with gffread (`-w` spliced transcripts by
   default, `-x` for CDS).
5. Build the dataset via `Dataset.from_components` (bowtie2 index, jellyfish
   k-mers, optional rRNA/tRNA blocklist, optional annotation tables).
6. Round-trip validate (every GTF transcript resolvable in the FASTA).
7. Write a `solar_intake.yaml` provenance manifest.
"""

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Literal

import click
import polars as pl
import pyfastx
import yaml
from loguru import logger
from pydantic import BaseModel, Field
from rich.console import Console
from rich.table import Table

from .dataset import Dataset
from .external_data import ExternalData

EXTERNAL_TOOLS: dict[str, str] = {
    "gffread": "conda install -c bioconda gffread",
    "bowtie2": "conda install -c bioconda bowtie2",
    "bowtie2-build": "conda install -c bioconda bowtie2",
    "jellyfish": "conda install -c bioconda kmer-jellyfish",
}

#: Tools probe design needs, as opposed to dataset construction.
DESIGN_TOOLS: tuple[str, ...] = ("bowtie2",)

_TOOL_VERSION_ARGS: dict[str, list[str]] = {
    "gffread": ["--version"],
    "bowtie2": ["--version"],
    "bowtie2-build": ["--version"],
    "jellyfish": ["--version"],
}

RESERVED_SPECIES = ("human", "mouse")


def check_external_tools(
    required: tuple[str, ...] = ("gffread", "bowtie2-build", "jellyfish"),
) -> dict[str, str]:
    """
    Verifies required external tools are on PATH and captures their versions.

    Raises a single aggregated error naming every missing tool with an
    install hint, so the user fixes the environment once, not tool-by-tool.
    """
    missing: list[str] = []
    versions: dict[str, str] = {}
    for tool in required:
        if shutil.which(tool) is None:
            missing.append(f"  {tool}: not found. Install with: {EXTERNAL_TOOLS.get(tool, '?')}")
            continue
        try:
            out = subprocess.run(
                [tool, *_TOOL_VERSION_ARGS.get(tool, ["--version"])],
                capture_output=True,
                text=True,
                timeout=30,
            )
            first_line = (out.stdout or out.stderr).strip().splitlines()
            versions[tool] = first_line[0] if first_line else "unknown"
        except Exception:  # noqa: BLE001 - version capture is best-effort
            versions[tool] = "unknown"
    if missing:
        raise RuntimeError("Missing required external tools:\n" + "\n".join(missing))
    return versions


def detect_annotation_format(path: Path) -> Literal["gtf", "gff3"]:
    """Peeks at the first data line to distinguish GTF (`key "value";`) from GFF3 (`key=value`)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attr = fields[8]
            if "=" in attr.split(";")[0] and ' "' not in attr:
                return "gff3"
            return "gtf"
    raise ValueError(f"No data lines found in annotation file {path}.")


class ValidationIssue(BaseModel):
    code: str
    severity: Literal["error", "warning"]
    message: str
    fix: str | None = None


class IngestValidationReport(BaseModel):
    """Validation results for a genome + annotation pair, serialized to validation_report.json."""

    genome: str
    annotation: str
    n_genes: int = 0
    n_transcripts: int = 0
    feature_counts: dict[str, int] = Field(default_factory=dict)
    transcripts_per_gene_median: float | None = None
    transcripts_per_gene_max: int | None = None
    gene_name_source: str = "gene_name"
    n_unstranded_transcripts: int = 0
    issues: list[ValidationIssue] = Field(default_factory=list)

    #: IDs of transcripts with no strand. Excluded from the JSON report - there
    #: can be tens of thousands - and written to `unstranded_transcripts.txt`
    #: instead, so targets can be checked against them.
    unstranded_transcript_ids: list[str] = Field(default_factory=list, exclude=True)

    def add(self, code: str, severity: Literal["error", "warning"], message: str, fix: str | None = None):
        self.issues.append(ValidationIssue(code=code, severity=severity, message=message, fix=fix))

    @property
    def errors(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors


def _truncate(values: list[str], n: int = 10) -> str:
    shown = ", ".join(map(str, values[:n]))
    return shown + (f", ... ({len(values) - n} more)" if len(values) > n else "")


def validate_gtf(gtf_path: Path, genome_fasta: Path) -> IngestValidationReport:
    """
    Validates a (already GTF-format, uncompressed-or-gz) annotation against its genome.

    Collects issues rather than raising: the caller decides whether errors
    abort. See the module docstring for the checks performed.
    """
    report = IngestValidationReport(genome=str(genome_fasta), annotation=str(gtf_path))

    try:
        df = ExternalData.parse_gtf(gtf_path, filters=None, strip_version=False)
    except ValueError as e:
        report.add(
            "PARSE_FAILED",
            "error",
            f"Annotation could not be parsed as GTF: {e}",
            fix="Normalize with `gffread <in> -T -o <out.gtf>` and retry.",
        )
        return report

    report.feature_counts = dict(
        df.group_by("feature").len().sort("feature").iter_rows()
    )

    tx = df.filter(pl.col("feature") == "transcript")
    exons = df.filter(pl.col("feature") == "exon")

    if tx.is_empty():
        report.add(
            "NO_TRANSCRIPT_ROWS",
            "error",
            f"No `transcript` feature rows found (features present: {sorted(report.feature_counts)}).",
            fix="Normalize with `gffread <in> -T -o <out.gtf>`, which emits transcript rows.",
        )
        return report

    # Required attributes on transcript/exon rows.
    for col, feat_df, feat in (("transcript_id", tx, "transcript"), ("gene_id", tx, "transcript")):
        n_null = feat_df[col].null_count()
        if n_null:
            report.add(
                "REQUIRED_ATTR_MISSING",
                "error",
                f"{n_null}/{feat_df.height} `{feat}` rows lack `{col}`.",
                fix="Every transcript row must carry gene_id and transcript_id (GTF2.2 spec).",
            )
    if not exons.is_empty() and (n_null := exons["transcript_id"].null_count()):
        report.add(
            "REQUIRED_ATTR_MISSING",
            "error",
            f"{n_null}/{exons.height} `exon` rows lack `transcript_id`.",
            fix="Exon rows must reference their transcript (GTF2.2 spec).",
        )

    # Duplicate transcript IDs.
    dup = (
        tx.group_by("transcript_id").len().filter(pl.col("len") > 1)["transcript_id"].drop_nulls().to_list()
    )
    if dup:
        report.add(
            "DUPLICATE_TRANSCRIPT_ID",
            "error",
            f"{len(dup)} transcript_id values appear on multiple transcript rows: {_truncate(dup)}.",
            fix="Transcript IDs must be unique; deduplicate the annotation.",
        )

    # Coordinate sanity.
    bad_coords = df.filter(pl.col("start") > pl.col("end")).height
    if bad_coords:
        report.add(
            "COORDINATES_INVERTED",
            "error",
            f"{bad_coords} rows have start > end.",
            fix="GTF coordinates are 1-based inclusive with start <= end regardless of strand.",
        )

    # Strand. Counted in transcripts, not rows: a row count mixes transcripts
    # with their exons and roughly doubles the figure, which reads as noise
    # rather than as a fifth of the annotation.
    transcripts = df.filter(pl.col("feature") == "transcript") if "feature" in df.columns else df
    unstranded = transcripts.filter(~pl.col("strand").is_in(["+", "-"]))
    ids = (
        unstranded["transcript_id"].drop_nulls().unique(maintain_order=True).to_list()
        if "transcript_id" in unstranded.columns
        else []
    )
    report.unstranded_transcript_ids = ids
    report.n_unstranded_transcripts = len(ids)
    if ids:
        # Counted here rather than from report.n_transcripts, which is not
        # populated until later in this function.
        total = transcripts.height
        share = f" ({len(ids) / total:.1%})" if total else ""
        report.add(
            "STRAND_MISSING",
            "warning",
            f"{len(ids)} transcript(s){share} have no strand (+/-). They are kept, but gffread "
            "extracts the plus-strand sequence for them, so probes for any that are really on "
            "the minus strand will be antisense and will not bind. Nothing downstream catches "
            "this: the sequences pass every filter and fail only at the bench.",
            fix="IDs are written to unstranded_transcripts.txt. Check your targets against it "
            "before designing, and drop or replace any that appear. Single-exon transcripts "
            "carry no splice signal, so strand cannot be recovered from the annotation alone.",
        )

    # seqname vs genome contigs - the classic chr1-vs-1/scaffold mismatch that
    # silently yields an EMPTY gffread output.
    try:
        genome_contigs = set(pyfastx.Fasta(str(genome_fasta)).keys())
    except Exception as e:  # noqa: BLE001
        report.add(
            "GENOME_UNREADABLE",
            "error",
            f"Could not read genome FASTA {genome_fasta}: {e}",
            fix="Provide a valid (optionally gzipped) genome FASTA.",
        )
        return report
    gtf_seqnames = set(df["seqname"].unique().to_list())
    missing_contigs = sorted(gtf_seqnames - genome_contigs)
    if missing_contigs:
        severity: Literal["error", "warning"] = (
            "error" if len(missing_contigs) == len(gtf_seqnames) else "warning"
        )
        report.add(
            "SEQNAME_MISMATCH",
            severity,
            f"{len(missing_contigs)}/{len(gtf_seqnames)} GTF seqnames absent from the genome: "
            f"{_truncate(missing_contigs)}. Genome contigs look like: "
            f"{_truncate(sorted(genome_contigs)[:5], 5)}.",
            fix="Use the annotation built FOR this assembly, or rename seqnames to match "
            "(e.g. chr1 vs 1 vs scaffold_1). Transcripts on missing contigs are silently dropped.",
        )

    # Transcripts with no exon rows (gffread emits nothing for them with -w).
    if not exons.is_empty():
        tx_ids = set(tx["transcript_id"].drop_nulls().to_list())
        exon_tx = set(exons["transcript_id"].drop_nulls().to_list())
        no_exons = sorted(tx_ids - exon_tx)
        if no_exons:
            report.add(
                "TRANSCRIPT_WITHOUT_EXONS",
                "warning",
                f"{len(no_exons)} transcripts have no exon rows and will be dropped by gffread: "
                f"{_truncate(no_exons)}.",
            )
    else:
        report.add(
            "NO_EXON_ROWS",
            "error",
            "No `exon` feature rows found; gffread cannot extract any transcript sequence.",
            fix="Normalize with `gffread <in> -T -o <out.gtf>` or use an annotation with exon rows.",
        )

    # ID hygiene: characters that break probe naming or FASTA keying.
    ids = tx["transcript_id"].drop_nulls().to_list()
    bad_ids = [i for i in ids if any(c in i for c in ('|', '"'))]
    if bad_ids:
        report.add(
            "ID_FORBIDDEN_CHARS",
            "error",
            f"{len(bad_ids)} transcript IDs contain forbidden characters (| or \"): {_truncate(bad_ids)}.",
            fix="Rename IDs before ingest; these characters break downstream naming/parsing.",
        )
    colon_ids = [i for i in ids if ":" in i]
    if colon_ids:
        report.add(
            "ID_COLON",
            "warning",
            f"{len(colon_ids)} transcript IDs contain ':' (e.g. {_truncate(colon_ids, 3)}); "
            "supported, but worth double-checking probe names downstream.",
        )

    # Counts + isoform distribution.
    report.n_transcripts = tx.height
    report.n_genes = tx["gene_id"].n_unique()
    per_gene = tx.group_by("gene_id").len()["len"]
    report.transcripts_per_gene_median = float(per_gene.median()) if per_gene.len() else None
    report.transcripts_per_gene_max = int(per_gene.max()) if per_gene.len() else None

    # Which gene-name fallback applies (informational).
    if "gene_name" in df.columns and tx["gene_name"].null_count() == 0:
        # parse_gtf always creates gene_name; find its actual source.
        raw_has_gene_name = tx["attribute"].str.contains(r'gene_name "').any()
        raw_has_name = tx["attribute"].str.contains(r'Name "').any()
        raw_has_gene = tx["attribute"].str.contains(r'(?:^|; ?)gene "').any()
        report.gene_name_source = (
            "gene_name"
            if raw_has_gene_name
            else "Name" if raw_has_name else "gene" if raw_has_gene else "gene_id (fallback)"
        )
        if report.gene_name_source == "gene_id (fallback)":
            report.add(
                "GENE_NAME_FALLBACK",
                "warning",
                "No gene-naming attribute (gene_name/Name/gene) found; gene names fall back to gene_id. "
                "Targets must then be selected by ID (or register an ortholog/alias annotation table).",
            )

    # gffread's own validation, folded in when available (banner lines skipped).
    if shutil.which("gffread"):
        res = subprocess.run(
            ["gffread", "-E", str(gtf_path), "-o", "/dev/null"],
            capture_output=True,
            text=True,
        )
        banner = ("Command line was:", "gffread ", "..loaded", ".. loaded")
        for line in dict.fromkeys(res.stderr.strip().splitlines()):  # dedupe, keep order
            line = line.strip()
            if line and not line.startswith(banner):
                report.add("GFFREAD_WARNING", "warning", f"gffread -E: {line}")

    return report


def print_report(report: IngestValidationReport, console: Console | None = None) -> None:
    """Pretty-prints a validation report."""
    console = console or Console()
    table = Table(title="SOLAR ingest validation", show_lines=False)
    table.add_column("Metric")
    table.add_column("Value")
    table.add_row("Genome", report.genome)
    table.add_row("Annotation", report.annotation)
    table.add_row("Genes", str(report.n_genes))
    table.add_row("Transcripts", str(report.n_transcripts))
    table.add_row("Features", ", ".join(f"{k}:{v}" for k, v in report.feature_counts.items()))
    table.add_row(
        "Isoforms/gene",
        f"median {report.transcripts_per_gene_median}, max {report.transcripts_per_gene_max}",
    )
    table.add_row("Gene-name source", report.gene_name_source)
    console.print(table)
    for issue in report.issues:
        style = "bold red" if issue.severity == "error" else "yellow"
        console.print(f"[{style}]{issue.severity.upper()} [{issue.code}][/]: {issue.message}")
        if issue.fix:
            console.print(f"  [dim]fix:[/] {issue.fix}")
    verdict = "[bold green]PASSED[/]" if report.ok else "[bold red]FAILED[/]"
    console.print(
        f"Validation {verdict}: {len(report.errors)} error(s), {len(report.warnings)} warning(s)."
    )


def _materialize_plain(path: Path, dest: Path) -> Path:
    """Copies `path` to `dest`, decompressing when gzipped. Returns dest."""
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as src, open(dest, "wb") as out:
            shutil.copyfileobj(src, out)
    elif path.resolve() != dest.resolve():
        shutil.copy(path, dest)
    return dest


def normalize_annotation(annotation: Path, dest_gtf: Path) -> Literal["gtf", "gff3"]:
    """
    Materializes the annotation as a plain-text GTF at `dest_gtf`.

    GFF3 inputs are converted with `gffread -T`; gzipped inputs are
    decompressed. Returns the detected source format.
    """
    fmt = detect_annotation_format(annotation)
    if fmt == "gff3":
        plain = dest_gtf.with_suffix(".input.gff3")
        _materialize_plain(annotation, plain)
        logger.info(f"Converting GFF3 -> GTF with gffread -T: {annotation.name} -> {dest_gtf.name}")
        res = subprocess.run(
            ["gffread", str(plain), "-T", "-o", str(dest_gtf)], capture_output=True, text=True
        )
        plain.unlink(missing_ok=True)
        if res.returncode != 0 or not dest_gtf.exists():
            raise RuntimeError(f"gffread -T failed converting {annotation}: {res.stderr.strip()[-2000:]}")
    else:
        _materialize_plain(annotation, dest_gtf)
    return fmt


def run_gffread(
    genome: Path, gtf: Path, out_fasta: Path, mode: Literal["transcripts", "cds"] = "transcripts"
) -> None:
    """
    Extracts transcript (spliced exons, `-w`) or CDS (`-x`) sequences.

    gffread writes FASTA headers carrying the exact transcript_id from the
    GTF, so with identity keys (strip_version=False) the GTF<->FASTA
    agreement holds by construction. The genome must be plain text (gffread
    needs random access); the caller decompresses .gz first.
    """
    flag = "-w" if mode == "transcripts" else "-x"
    cmd = ["gffread", flag, str(out_fasta), "-g", str(genome), str(gtf)]
    logger.info(f"Extracting {mode}: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"gffread extraction failed: {res.stderr.strip()[-2000:]}")
    if not out_fasta.exists() or out_fasta.stat().st_size == 0:
        raise RuntimeError(
            f"gffread produced no sequences at {out_fasta}. This usually means the GTF seqnames "
            "do not match the genome contig names - run with --validate-only for details."
        )


def extract_biotype_blocklist(gtf_df: pl.DataFrame, transcripts_fasta: Path, biotypes: list[str], out: Path) -> Path | None:
    """
    Writes a FASTA of transcripts whose biotype matches `biotypes` (e.g. rRNA, tRNA, snoRNA).

    Looks for the first biotype-ish column present (gene_biotype,
    transcript_biotype, gene_type, transcript_type). Returns None (with a
    warning) if none exists or nothing matches.
    """
    biotype_col = next(
        (c for c in ("gene_biotype", "transcript_biotype", "gene_type", "transcript_type") if c in gtf_df.columns),
        None,
    )
    if biotype_col is None:
        logger.warning(
            f"--blocklist-biotypes requested ({biotypes}) but the GTF has no biotype column; skipping."
        )
        return None
    ids = set(
        gtf_df.filter(pl.col("feature") == "transcript")
        .filter(pl.col(biotype_col).is_in(biotypes))["transcript_id"]
        .drop_nulls()
        .to_list()
    )
    if not ids:
        logger.warning(f"No transcripts with {biotype_col} in {biotypes}; skipping biotype blocklist.")
        return None
    fa = pyfastx.Fasta(str(transcripts_fasta))
    written = 0
    with open(out, "w") as fh:
        for name in ids:
            try:
                seq = fa[name].seq
            except KeyError:
                continue
            fh.write(f">{name}\n{seq}\n")
            written += 1
    if not written:
        logger.warning("Biotype blocklist transcripts not found in extracted FASTA; skipping.")
        out.unlink(missing_ok=True)
        return None
    logger.info(f"Biotype blocklist: {written} sequences ({biotype_col} in {biotypes}) -> {out.name}")
    return out


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_intake_manifest(
    dataset_dir: Path,
    *,
    species: str,
    genome: Path,
    annotation: Path,
    annotation_format: str,
    extract_mode: str,
    fasta_key_regex: str,
    strip_version: bool,
    tool_versions: dict[str, str],
    report: IngestValidationReport,
    blocklist_files: list[str],
    annotation_tables: dict[str, str],
    argv: list[str],
    gene_name_column: str | None = None,
) -> Path:
    """Writes solar_intake.yaml: auto-filled provenance + QC; operator completes the stubs."""
    manifest = {
        "manifest_version": 2,
        "species": {"display_name": species, "internal_id": species, "taxonomy_id": None},
        "inputs": {
            "genome_fasta": str(genome),
            "genome_sha256": _sha256(genome),
            "annotation": str(annotation),
            "annotation_sha256": _sha256(annotation),
            "annotation_format": annotation_format,
            "blocklist_fasta": blocklist_files,
            "annotation_tables": annotation_tables,
            "gene_name_column": gene_name_column,
        },
        "provenance": {
            "assembly_name": "",
            "assembly_source": "",
            "assembly_release_date": "",
            "annotation_method": "",
            "annotation_build_date": "",
            "data_owner": "",
        },
        "processing": {
            "ingest_command": " ".join(argv),
            "extract_mode": extract_mode,
            "fasta_key_regex": fasta_key_regex,
            "strip_version": strip_version,
            "software_versions": tool_versions,
        },
        "quality_control": {
            "n_genes": report.n_genes,
            "n_transcripts": report.n_transcripts,
            "transcripts_per_gene_median": report.transcripts_per_gene_median,
            "gene_name_source": report.gene_name_source,
            "validation_errors": len(report.errors),
            "validation_warnings": len(report.warnings),
            "passed_for_probe_generation": report.ok,
        },
        "review": {"reviewer": "", "review_date": "", "notes": ""},
    }
    out = dataset_dir / "solar_intake.yaml"
    out.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return out


def validate_dataset(ds: Dataset, gtf_transcript_ids: list[str]) -> None:
    """
    Round-trip validation of a freshly built dataset.

    Every transcript_id in the GTF must resolve in the extracted FASTA, and
    vice versa (both directions reported). Also confirms the k-mer file and
    bowtie2 index exist and parse.
    """
    fasta_keys = set(ds.data.fa.keys())
    gtf_ids = set(gtf_transcript_ids)
    in_gtf_not_fasta = sorted(gtf_ids - fasta_keys)
    in_fasta_not_gtf = sorted(fasta_keys - gtf_ids)
    if in_gtf_not_fasta:
        raise ValueError(
            f"{len(in_gtf_not_fasta)}/{len(gtf_ids)} GTF transcripts missing from the extracted FASTA "
            f"(dropped by gffread - no exons, bad contig, or strand issues): "
            f"{_truncate(in_gtf_not_fasta)}. See validation_report.json warnings."
        )
    if in_fasta_not_gtf:
        logger.warning(
            f"{len(in_fasta_not_gtf)} FASTA records not present as GTF transcripts "
            f"(harmless but unexpected): {_truncate(in_fasta_not_gtf)}"
        )
    # Spot-check sequence retrieval through the normal lookup path.
    for tid in list(gtf_ids)[:5]:
        seq = ds.data.get_seq(tid, convert=False)
        if not seq:
            raise ValueError(f"Empty sequence for transcript {tid}.")
    if ds.kmer18 is None:
        raise ValueError("18-mer jellyfish file is missing.")
    if ds.kmer18.is_empty():
        logger.warning(
            "18-mer file has no k-mers above threshold (normal for small transcriptomes); "
            "repeat-region filtering will be a no-op."
        )
    _ = ds.data.bowtie2_index  # raises FileNotFoundError if missing
    logger.info(
        f"Dataset validation passed: {len(fasta_keys)} transcripts, "
        f"k-mers and bowtie2 index present{', blocklist active' if ds.trna_rna_kmers else ''}."
    )


@click.command()
@click.argument("dataset_dir", type=click.Path(file_okay=False, path_type=Path))
@click.option(
    "--genome",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Genome FASTA (optionally .gz).",
)
@click.option(
    "--gtf",
    "annotation",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Annotation in GTF or GFF3 (optionally .gz). GFF3 is converted with gffread -T.",
)
@click.option("--species", "-s", required=True, type=str, help="Species name (metadata; not 'human'/'mouse').")
@click.option(
    "--extract",
    "extract_mode",
    type=click.Choice(["transcripts", "cds"]),
    default="transcripts",
    show_default=True,
    help="What gffread extracts: spliced transcripts with UTRs (-w) or CDS only (-x).",
)
@click.option(
    "--rrna-fasta",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="rRNA FASTA for the 15-mer blocklist. Repeatable.",
)
@click.option(
    "--trna-fasta",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="tRNA FASTA for the 15-mer blocklist. Repeatable.",
)
@click.option(
    "--blocklist-biotypes",
    default=None,
    help="Comma-separated biotypes to auto-blocklist from the GTF's biotype column "
    "(e.g. 'rRNA,tRNA,snoRNA'). Requires a biotype attribute in the GTF.",
)
@click.option(
    "--gene-name-column",
    "gene_name_column",
    type=str,
    default=None,
    metavar="COLUMN",
    help="Column of a registered annotation table holding the gene names you want to write "
    "in target lists, e.g. --gene-name-column Hsapiens_gene_name. Name lookup then uses only "
    "that column. Cells holding a comma-separated list count as one name per entry.",
)
@click.option(
    "--annotation-table",
    "annotation_tables",
    multiple=True,
    metavar="NAME=PATH",
    help="Register a named annotation table (parquet/csv/tsv with transcript_id and/or "
    "gene_id column), e.g. --annotation-table orthologs=orthologs.tsv. Repeatable.",
)
@click.option(
    "--keep-genome",
    is_flag=True,
    help="Copy the genome FASTA into the dataset directory (for genome-mode simulation). "
    "Off by default: genomes are large; provenance (path+sha256) is always recorded.",
)
@click.option("--fasta-key-regex", default=r"^(\S+)", show_default=True, help="FASTA header key regex.")
@click.option(
    "--strip-version/--no-strip-version",
    default=False,
    show_default=True,
    help="Strip trailing .N from IDs. Default OFF for ingest: gffread headers match the GTF "
    "verbatim, and stripping merges StringTie-style isoforms (STRG.1.1/STRG.1.2).",
)
@click.option("--validate-only", is_flag=True, help="Run validation and write the report; build nothing.")
@click.option("--overwrite", is_flag=True, help="Overwrite existing dataset artifacts.")
def ingest(
    dataset_dir: Path,
    genome: Path,
    annotation: Path,
    species: str,
    extract_mode: Literal["transcripts", "cds"],
    rrna_fasta: tuple[Path, ...],
    trna_fasta: tuple[Path, ...],
    blocklist_biotypes: str | None,
    annotation_tables: tuple[str, ...],
    gene_name_column: str | None,
    keep_genome: bool,
    fasta_key_regex: str,
    strip_version: bool,
    validate_only: bool,
    overwrite: bool,
):
    """Ingest a genome FASTA + GTF/GFF3 into a SOLAR probe-design dataset."""
    console = Console()
    if species.lower() in RESERVED_SPECIES:
        raise click.BadParameter(
            f"Species {species!r} is reserved for reference datasets (mkprobes prepare). "
            "Use a distinct name for a custom dataset.",
            param_hint="--species",
        )
    tables: dict[str, str | Path] = {}
    for spec in annotation_tables:
        name, sep, table_path = spec.partition("=")
        if not sep or not name or not table_path:
            raise click.BadParameter(f"Expected NAME=PATH, got {spec!r}", param_hint="--annotation-table")
        tables[name] = Path(table_path)

    tool_versions = check_external_tools()
    dataset_dir.mkdir(parents=True, exist_ok=True)

    # 1. Materialize a plain-text GTF (converts GFF3, decompresses .gz).
    #
    # The guard is on `dataset.json`, the marker a build completed, not on the
    # annotation copy. `--validate-only` also writes annotation.gtf, so guarding
    # on that made the documented flow - validate, then ingest - fail on its own
    # leftovers.
    dest_gtf = dataset_dir / "annotation.gtf"
    built_marker = dataset_dir / "dataset.json"
    if built_marker.exists() and not overwrite and annotation.resolve() != dest_gtf.resolve():
        raise click.ClickException(
            f"{dataset_dir} already holds a built dataset. Pass --overwrite to replace it."
        )
    annotation_format = normalize_annotation(annotation, dest_gtf)

    # 2. Genome must be plain text for gffread's random access.
    if genome.suffix == ".gz":
        plain_genome = dataset_dir / genome.with_suffix("").name if keep_genome else (
            dataset_dir / f".tmp_{genome.with_suffix('').name}"
        )
        logger.info(f"Decompressing genome to {plain_genome}")
        _materialize_plain(genome, plain_genome)
    else:
        plain_genome = genome

    try:
        # 3. Validate and report.
        report = validate_gtf(dest_gtf, plain_genome)
        (dataset_dir / "validation_report.json").write_text(report.model_dump_json(indent=2))
        # Written even with --validate-only: knowing which targets are unsafe is
        # the point of validating before committing to a build.
        if report.unstranded_transcript_ids:
            unstranded_path = dataset_dir / "unstranded_transcripts.txt"
            unstranded_path.write_text("\n".join(report.unstranded_transcript_ids) + "\n")
            logger.info(
                f"{len(report.unstranded_transcript_ids)} unstranded transcript ID(s) "
                f"written to {unstranded_path}."
            )
        print_report(report, console)
        if not report.ok:
            raise click.ClickException(
                f"Validation failed with {len(report.errors)} error(s); dataset not built. "
                f"Report: {dataset_dir / 'validation_report.json'}"
            )
        if validate_only:
            console.print("[green]--validate-only: stopping after validation.[/]")
            return

        # 4. Extract transcript sequences.
        transcripts_fasta = dataset_dir / "transcripts.fasta"
        if transcripts_fasta.exists() and not overwrite:
            logger.info(f"{transcripts_fasta} exists; skipping extraction (pass --overwrite to redo).")
        else:
            run_gffread(plain_genome, dest_gtf, transcripts_fasta, extract_mode)

        # 5. Optional biotype-derived blocklist.
        blocklist_files = [*rrna_fasta, *trna_fasta]
        gtf_df = ExternalData.parse_gtf(dest_gtf, filters=None, strip_version=strip_version)
        if blocklist_biotypes:
            biotype_fa = extract_biotype_blocklist(
                gtf_df,
                transcripts_fasta,
                [b.strip() for b in blocklist_biotypes.split(",") if b.strip()],
                dataset_dir / "biotype_blocklist.fasta",
            )
            if biotype_fa:
                blocklist_files.append(biotype_fa)

        # 6. Build the dataset (bowtie2 + jellyfish + definitions).
        ds = Dataset.from_components(
            dataset_dir,
            transcripts_fasta,
            species=species,
            gtf_file=dest_gtf,
            fasta_key_regex=fasta_key_regex,
            strip_version=strip_version,
            blocklist_fasta=blocklist_files or None,
            genome_fasta=plain_genome if keep_genome else None,
            annotations=tables or None,
            gene_name_column=gene_name_column,
            interactive=False,
            overwrite=overwrite,
        )

        # 7. Round-trip validation.
        tx_ids = (
            gtf_df.filter(pl.col("feature") == "transcript")["transcript_id"].drop_nulls().to_list()
        )
        validate_dataset(ds, tx_ids)

        # 8. Provenance manifest.
        manifest_path = write_intake_manifest(
            dataset_dir,
            species=species,
            genome=genome,
            annotation=annotation,
            annotation_format=annotation_format,
            extract_mode=extract_mode,
            fasta_key_regex=fasta_key_regex,
            strip_version=strip_version,
            tool_versions=tool_versions,
            report=report,
            blocklist_files=[str(p) for p in blocklist_files],
            annotation_tables={k: str(v) for k, v in tables.items()},
            gene_name_column=gene_name_column,
            argv=sys.argv,
        )
        console.print(
            f"[bold green]Ingest complete.[/] Dataset at {dataset_dir}; "
            f"manifest at {manifest_path.name} (complete the provenance stubs). "
            f"Next: `mkprobes transcripts {dataset_dir} <genes.txt> --longest`."
        )
    finally:
        if plain_genome.name.startswith(".tmp_"):
            plain_genome.unlink(missing_ok=True)
            Path(str(plain_genome) + ".fai").unlink(missing_ok=True)
