import shutil
from functools import cache
from pathlib import Path
from typing import Literal, Sequence, cast

import click
import polars as pl
import pyfastx
from loguru import logger
from pydantic import BaseModel

from ..utils._alignment import jellyfish
from ..utils.sequtils import kmers
from .external_data import ExternalData, ExternalDataDefinition


def parse_jellyfish(path: Path | str) -> pl.DataFrame:
    """
    Parse a jellyfish output file.

    An empty (but existing) file is legitimate: jellyfish ran to completion and
    no k-mer cleared the count threshold (e.g. `-L 10` on a small
    transcriptome). Returns an empty frame with the right schema in that case.
    A *missing* file is still an error (raised upstream by `ExternalData.kmer`).
    """
    try:
        return pl.read_csv(path, separator=" ", has_header=False, new_columns=["kmer", "count"])
    except pl.exceptions.NoDataError:
        logger.warning(
            f"Jellyfish file {path} has no k-mers above threshold (normal for small transcriptomes)."
        )
        return pl.DataFrame({"kmer": [], "count": []}, schema={"kmer": pl.Utf8, "count": pl.Int64})


ANNOTATION_JOIN_COLUMNS = ("transcript_id", "gene_id")


def _read_annotation_table(path: Path | str) -> pl.DataFrame:
    """Reads an annotation table by extension: .parquet, .csv, or .tsv/.txt."""
    path = Path(path)
    match path.suffix:
        case ".parquet":
            return pl.read_parquet(path)
        case ".csv":
            return pl.read_csv(path)
        case ".tsv" | ".txt":
            return pl.read_csv(path, separator="\t")
        case _:
            raise ValueError(
                f"Unsupported annotation table format {path.suffix!r} for {path}. "
                "Use .parquet, .csv, or .tsv."
            )


def _validate_annotation_table(name: str, path: Path | str) -> None:
    """
    Validates that an annotation table is readable and joinable.

    Each table must carry at least one of the join columns
    (`transcript_id`/`gene_id`) so it can be joined against the dataset.
    """
    df = _read_annotation_table(path)
    if not any(c in df.columns for c in ANNOTATION_JOIN_COLUMNS):
        raise ValueError(
            f"Annotation table {name!r} ({path}) has no join column: expected at least one of "
            f"{ANNOTATION_JOIN_COLUMNS}, got columns {df.columns}."
        )


class DatasetDefinition(BaseModel):
    """
    Defines the structure for storing dataset configuration.

    This Pydantic model is used to serialize and deserialize dataset metadata,
    primarily for saving and loading dataset configurations from a JSON file.
    It ensures that the dataset information is consistent and correctly structured.

    - Ensure all required fields (species, external_data) are present when creating an instance.
    - `blocklist_kmer_name` is optional (added for custom-species datasets): the
      name of a Jellyfish k-mer file built from rRNA/tRNA/other unwanted
      sequences, used by `Dataset.check_kmers` to reject probe candidates.
    - `genome_fasta_name` is optional: the source genome FASTA, when kept inside
      the dataset directory (e.g. `mkprobes ingest --keep-genome`). The probe
      pipeline itself is transcriptome-based; the genome is used by downstream
      steps such as genome-mode probe simulation for non-model species.
    - `annotations` is an optional registry of named annotation tables
      (parquet/csv/tsv file names relative to the dataset directory). Each table
      must carry a `transcript_id` and/or `gene_id` column to join against the
      dataset. Typical uses: ortholog mappings, functional annotation
      (eggNOG/InterPro), expression tables (FPKM).
    - All optional fields have back-compat defaults: older `dataset.json` files
      load unchanged.
    """

    species: str
    external_data: dict[str, ExternalDataDefinition]
    blocklist_kmer_name: str | None = None
    genome_fasta_name: str | None = None
    annotations: dict[str, str] = {}


class Dataset:
    """
    Represents a dataset for probe design, including sequence data and k-mer information.

    This class manages the core data files (FASTA, k-mer counts) and provides
    methods to access and process them. It can be initialized either by building
    from component files (FASTA) or by loading an existing dataset from a folder.

    Intended Use:
        - Creating new datasets from FASTA files.
        - Loading pre-existing datasets for analysis or probe design.
        - Accessing k-mer sets and performing k-mer-based checks.

    Potential Pitfalls:
        - If `from_components` is not used or if files are manually moved/deleted,
          the dataset might become inconsistent, leading to `FileNotFoundError` or
          `ValueError` during operations.
        - Ensure `kmer18_path` and `trna_rna_kmers_path` point to valid Jellyfish output
          files if provided; otherwise, related functionalities might not work as expected.
        - The `overwrite` flag in `from_components` should be used with caution to avoid
          accidental data loss.
    """

    def __init__(
        self,
        path: str | Path,
        external_data: ExternalData,
        species: str | None = None,
        kmer18_path: str | Path | None = None,
        trna_rna_kmers_path: str | Path | None = None,
        genome_fasta_path: str | Path | None = None,
        annotation_paths: dict[str, Path] | None = None,
    ):
        """
        Initializes a Dataset object.

        Args:
            path: The base directory path for the dataset.
            external_data: An `ExternalData` object containing sequence and annotation data.
            species: The species name (e.g., "human", "mouse"). Defaults to None.
            kmer18_path: Path to the 18-mer Jellyfish count file. Defaults to None.
            trna_rna_kmers_path: Path to a Jellyfish count file for tRNA/rRNA kmers.
                Defaults to None.
            genome_fasta_path: Optional path to the source genome FASTA (used by
                downstream steps such as genome-mode simulation, not by the
                probe pipeline itself).
            annotation_paths: Optional registry of named annotation tables
                (name -> file path); load with `self.annotation(name)`.
        """
        self.species = species
        self.path = Path(path)
        self.data = external_data
        self.kmer18 = parse_jellyfish(kmer18_path) if kmer18_path else None
        self.kmerset = set(self.kmer18["kmer"] if self.kmer18 is not None else [])
        self.trna_rna_kmers = (
            set(parse_jellyfish(trna_rna_kmers_path)["kmer"]) if trna_rna_kmers_path else None
        )
        self.genome_fasta_path = Path(genome_fasta_path) if genome_fasta_path else None
        self.annotation_paths: dict[str, Path] = dict(annotation_paths or {})

        # For backwards compatibility
        self.gencode = self.data
        self.ensembl: ExternalData | None = None

    def annotation(self, name: str) -> pl.DataFrame:
        """
        Loads a named annotation table registered in `dataset.json`.

        Tables are parquet/csv/tsv files carrying a `transcript_id` and/or
        `gene_id` column for joining against the dataset (e.g. ortholog
        mappings, functional annotation, expression/FPKM tables).

        Raises:
            KeyError: If `name` is not a registered annotation table.
        """
        try:
            table_path = self.annotation_paths[name]
        except KeyError:
            raise KeyError(
                f"No annotation table named {name!r} in this dataset. "
                f"Available: {sorted(self.annotation_paths) or 'none'}"
            ) from None
        return _read_annotation_table(table_path)

    @classmethod
    def from_components(
        cls,
        path: str | Path,
        fasta_file: str | Path,
        *,
        species: str,
        gtf_file: str | Path | None = None,
        fasta_key_regex: str = r"^(\S+)",
        strip_version: bool = True,
        blocklist_fasta: Sequence[str | Path] | None = None,
        genome_fasta: str | Path | None = None,
        annotations: dict[str, str | Path] | None = None,
        interactive: bool = True,
        overwrite: bool = False,
    ):
        """
        Creates a new dataset from a FASTA file (and optionally a GTF).

        This method will:
        1. Create the dataset directory at `path`.
        2. Copy the `fasta_file` (and `gtf_file`, if given) into this directory.
        3. Initialize `ExternalData` (parsing the GTF into a cache when provided).
        4. Build Bowtie2 index for the FASTA file.
        5. Run Jellyfish to count k-mers from the FASTA file.
        6. Optionally build a 15-mer blocklist (`blocklist15.jf`) from
           rRNA/tRNA FASTA files, used to reject probe candidates.
        7. Create a `dataset.json` file with the dataset definition.

        Args:
            path: The directory where the dataset will be created.
            fasta_file: Path to the input (transcriptome) FASTA file.
            species: The species name for this dataset.
            gtf_file: Optional GTF annotation. When given, gene-name lookups and
                transcript selection work on this dataset instead of requiring
                raw FASTA record IDs.
            fasta_key_regex: Regex extracting the lookup key from FASTA headers.
            strip_version: If True (historical default), trailing `.N` version
                suffixes are stripped from FASTA keys and GTF IDs. Use False for
                de novo annotations whose IDs embed meaningful dots (StringTie
                `STRG.1.1`, AUGUSTUS `g1.t1`).
            blocklist_fasta: FASTA files of sequences to blocklist (typically
                rRNA/tRNA); their 15-mers are counted into `blocklist15.jf`.
            genome_fasta: Optional source genome FASTA to keep inside the
                dataset directory (copied in). Not used by the probe pipeline
                itself, but enables downstream genome-mode steps (e.g.
                simulation for non-model species).
            annotations: Optional named annotation tables (name -> file path,
                parquet/csv/tsv). Copied into the dataset directory and
                validated: each must carry a `transcript_id` and/or `gene_id`
                column. Load later with `dataset.annotation(name)`.
            interactive: If False, skip confirmation prompts before the
                Bowtie2/Jellyfish builds (for scripted flows like ingest).
            overwrite: If True, allows overwriting existing files (e.g., Jellyfish output).
                Defaults to False.

        Returns:
            An instance of the `Dataset` class.
        """
        path = Path(path)
        fasta_file = Path(fasta_file)

        logger.info(f"Creating dataset at {path}")
        path.mkdir(exist_ok=True, parents=True)
        new_path = path / fasta_file.name
        logger.info(f"Copying {fasta_file} to {new_path}")
        if fasta_file.resolve() != new_path.resolve():
            shutil.copy(fasta_file, new_path)
        del fasta_file

        gtf_name: str | None = None
        if gtf_file is not None:
            gtf_file = Path(gtf_file)
            new_gtf = path / gtf_file.name
            if gtf_file.resolve() != new_gtf.resolve():
                logger.info(f"Copying {gtf_file} to {new_gtf}")
                shutil.copy(gtf_file, new_gtf)
            gtf_name = new_gtf.name

        cache_path = path / new_path.with_suffix(".parquet").name
        external_data = ExternalData(
            cache=cache_path,
            fasta=new_path,
            gtf_path=path / gtf_name if gtf_name else None,
            regen_cache=overwrite,
            fasta_key_regex=fasta_key_regex,
            strip_version=strip_version,
        )
        external_data.bowtie_build(overwrite=overwrite, interactive=interactive)
        # Adaptive initial hash: jellyfish grows the hash as needed, so a
        # smaller initial allocation is correct for any input and avoids a
        # multi-second 10G-hash init on small transcriptomes.
        hash_size = f"{max(new_path.stat().st_size // 4, 10_000_000)}"
        external_data.run_jellyfish(overwrite=overwrite, interactive=interactive, hash_size=hash_size)

        blocklist_kmer_name: str | None = None
        if blocklist_fasta:
            blocklist_path = path / "blocklist15.jf"
            if not blocklist_path.exists() or overwrite:
                seqs: list[str] = []
                for f in blocklist_fasta:
                    try:
                        seqs.extend(record[1] for record in pyfastx.Fastx(str(f)))
                    except RuntimeError as e:
                        raise ValueError(
                            f"No sequences found in blocklist FASTA file {f} (empty or not FASTA)."
                        ) from e
                if not seqs:
                    raise ValueError(f"No sequences found in blocklist FASTA file(s): {blocklist_fasta}")
                logger.info(
                    f"Building 15-mer blocklist from {len(seqs)} sequences "
                    f"({len(blocklist_fasta)} file(s)) -> {blocklist_path.name}"
                )
                jellyfish(
                    seqs,
                    blocklist_path,
                    15,
                    hash_size=f"{max(sum(map(len, seqs)) * 2, 10_000_000)}",
                )
            blocklist_kmer_name = blocklist_path.name

        genome_fasta_name: str | None = None
        if genome_fasta is not None:
            genome_fasta = Path(genome_fasta)
            new_genome = path / genome_fasta.name
            if genome_fasta.resolve() != new_genome.resolve():
                logger.info(f"Copying genome {genome_fasta} to {new_genome}")
                shutil.copy(genome_fasta, new_genome)
            genome_fasta_name = new_genome.name

        annotation_names: dict[str, str] = {}
        if annotations:
            for name, table_file in annotations.items():
                table_file = Path(table_file)
                _validate_annotation_table(name, table_file)
                new_table = path / table_file.name
                if table_file.resolve() != new_table.resolve():
                    logger.info(f"Copying annotation table {name!r}: {table_file} -> {new_table}")
                    shutil.copy(table_file, new_table)
                annotation_names[name] = new_table.name

        Path(path / "dataset.json").write_text(
            DatasetDefinition(
                external_data={
                    "default": ExternalDataDefinition(
                        fasta_name=new_path.name,
                        cache_name=cache_path.name if gtf_name else None,
                        gtf_name=gtf_name,
                        bowtie2_index_name=external_data.bowtie2_index.name,
                        kmer18_name=external_data.kmer.name,
                        fasta_key_regex=fasta_key_regex,
                        strip_version=strip_version,
                    )
                },
                species=species,
                blocklist_kmer_name=blocklist_kmer_name,
                genome_fasta_name=genome_fasta_name,
                annotations=annotation_names,
            ).model_dump_json()
        )

        return cls(
            path=path,
            external_data=external_data,
            species=species,
            kmer18_path=external_data.kmer,
            trna_rna_kmers_path=path / blocklist_kmer_name if blocklist_kmer_name else None,
            genome_fasta_path=path / genome_fasta_name if genome_fasta_name else None,
            annotation_paths={k: path / v for k, v in annotation_names.items()},
        )

    @classmethod
    def from_folder(cls, path: Path):
        """
        Loads an existing dataset from a specified folder.

        The folder must contain a `dataset.json` file that defines the dataset structure
        and references the necessary data files (FASTA, k-mer counts, etc.).

        Args:
            path: The path to the dataset folder.

        Returns:
            An instance of the `Dataset` class.

        Raises:
            FileNotFoundError: If `dataset.json` is not found in the specified path.
        """
        path = Path(path)
        if not (path / "dataset.json").exists():
            raise FileNotFoundError(f"Path {path} does not exist. Please create a dataset first.")

        definition = DatasetDefinition.model_validate_json((path / "dataset.json").read_text())
        external_data = ExternalData.from_definition(path, definition.external_data["default"])
        return cls(
            path=path,
            external_data=external_data,
            species=definition.species,
            kmer18_path=path / external_data.kmer,
            trna_rna_kmers_path=(
                path / definition.blocklist_kmer_name if definition.blocklist_kmer_name else None
            ),
            genome_fasta_path=(
                path / definition.genome_fasta_name if definition.genome_fasta_name else None
            ),
            annotation_paths={k: path / v for k, v in definition.annotations.items()},
        )

    def check_kmers(self, seq: str):
        """
        Checks if any blocklist k-mers from the input sequence are present in the `trna_rna_kmers` set.

        This is typically used to filter out sequences that might originate from
        tRNAs or rRNAs, based on a pre-compiled set of common kmers from these RNA types.

        Args:
            seq: The nucleotide sequence to check.

        Returns:
            True if any k-mer from the sequence is found in `self.trna_rna_kmers`,
            False otherwise. Logs a warning if `self.trna_rna_kmers` is not set.
        """
        if not self.trna_rna_kmers:
            logger.warning("No tRNA-RNA kmers found. Skipping.")
            return False

        k = len(next(iter(self.trna_rna_kmers)))
        return any(x in self.trna_rna_kmers for x in kmers(seq, k))

    @property
    @cache
    def appris(self):
        raise NotImplementedError("appris not implemented")


class ReferenceDataset(Dataset):
    """
    Represents a pre-configured reference dataset, typically for human or mouse.

    This class inherits from `Dataset` and is specialized for known reference
    genomes/transcriptomes where specific file names and structures are expected
    (e.g., gencode.gtf.gz, ensembl.gtf.gz, appris_data.principal.txt).

    Intended Use:
        - Working with standardized human or mouse reference data provided by
          the fishtools package or a similar pre-packaged dataset.
        - Accessing specific annotation files like APPRIS data.

    Potential Pitfalls:
        - The constructor expects a specific directory structure and file names
          within the provided `path`. If these files are missing or named differently,
          `FileNotFoundError` will be raised.
        - Currently, only "human" and "mouse" species are explicitly supported,
          and using other species names will raise a `ValueError`.
        - If APPRIS data (`appris_data.principal.txt`) is not found, the `appris`
          property will return `None` and log a warning.
    """

    def __init__(self, path: Path | str):
        """
        Initializes a ReferenceDataset object.

        This constructor assumes a specific directory structure and file naming convention
        within the `path` for standard reference data (e.g., GENCODE GTF, Ensembl GTF,
        precomputed k-mer files).

        Args:
            path: The base directory path for the reference dataset.
                The name of this directory is used to infer the species (e.g., "human", "mouse").

        Raises:
            FileNotFoundError: If the `path` does not exist or if essential files
                (e.g., gencode.gtf.gz for human/mouse) are missing.
            ValueError: If the inferred species is not "human" or "mouse".
        """
        self.path = path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Path {path} does not exist.")

        super().__init__(
            path,
            external_data=ExternalData(
                cache=self.path / "gencode.parquet",
                gtf_path=self.path / "gencode.gtf.gz",
                fasta=self.path / "txome.fasta",
                fasta_key_regex=r"^([^.\s]+)",
            ),
            kmer18_path=self.path / "cdna18.jf",
            trna_rna_kmers_path=self.path / "r_t_snorna15.jf",
            species=cast(Literal["human", "mouse"], self.path.name),
        )
        if self.species not in ("human", "mouse"):
            raise ValueError(f"Species not human or mouse, got {self.species}")

        try:
            self.ensembl = ExternalData(
                cache=self.path / "ensembl.parquet",
                gtf_path=self.path / "ensembl.gtf.gz",
                fasta=self.path / "txome.fasta",
            )
        except FileNotFoundError:
            if self.species in ("human", "mouse"):
                raise FileNotFoundError(f"No GENCODE data found for {self.species}.")
            self.ensembl = None

    @property
    @cache
    def appris(self):
        """
        Loads APPRIS principal isoform data from 'appris_data.principal.txt'.

        APPRIS (Annotation of Principal Splice Isoforms) data helps identify the
        main functional transcript(s) for a gene. This method attempts to load
        this data from a standard file name within the dataset path.

        The expected file format is a tab-separated file with columns:
        gene_name, gene_id, transcript_id, ccds, annotation.

        Returns:
            A Polars DataFrame containing the APPRIS data if the file is found
            and successfully parsed. Returns `None` and logs a warning if the
            file is not found.
        """
        try:
            return pl.read_csv(
                self.path / "appris_data.principal.txt",
                separator="\t",
                has_header=False,
                new_columns=[
                    "gene_name",
                    "gene_id",
                    "transcript_id",
                    "ccds",
                    "annotation",
                ],
            )
        except FileNotFoundError:
            logger.warning("No APPRIS data found.")
            return None


def load_dataset(path: Path | str) -> "Dataset | ReferenceDataset":
    """
    Resolve a dataset path into either a generic dataset or a reference dataset.

    Resolution order:
    1. If `dataset.json` exists, load as a generic `Dataset` (custom species).
    2. If the folder name is `human` or `mouse`, load as a `ReferenceDataset`.
    3. Otherwise fail with guidance on how to create a dataset.
    """
    path = Path(path)
    if (path / "dataset.json").exists():
        return Dataset.from_folder(path)

    if path.name in ("human", "mouse"):
        return ReferenceDataset(path)

    raise FileNotFoundError(
        f"Path {path} is not a recognized probe dataset. Expected a `dataset.json` "
        "(custom species) or a `human`/`mouse` reference dataset folder. Create one with "
        "`mkprobes ingest <dir> --genome <genome.fa> --gtf <annotation.gtf> --species <name>` "
        "(genome + annotation) or "
        "`mkprobes create-dataset <dir> --fasta <transcriptome.fa> --species <name>` "
        "(transcriptome FASTA)."
    )


@click.command()
@click.argument("path", type=click.Path(dir_okay=True, file_okay=False, path_type=Path))
@click.option("--fasta", type=click.Path(dir_okay=False, file_okay=True, path_type=Path))
@click.option("--species", "-s", type=str, help="Species name for metadata")
@click.option(
    "--gtf",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional GTF annotation matching the FASTA. Enables gene-name lookups "
    "and transcript selection for this dataset.",
)
@click.option(
    "--blocklist-fasta",
    "blocklist_fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    help="FASTA of sequences to blocklist (rRNA/tRNA). Repeatable. "
    "15-mers from these sequences are used to reject probe candidates.",
)
@click.option(
    "--fasta-key-regex",
    default=r"^(\S+)",
    show_default=True,
    help="Regex extracting the lookup key from FASTA headers.",
)
@click.option(
    "--strip-version/--no-strip-version",
    default=True,
    show_default=True,
    help="Strip trailing .N version suffixes from IDs. Use --no-strip-version for "
    "de novo annotations whose IDs embed meaningful dots (StringTie STRG.1.1, AUGUSTUS g1.t1).",
)
@click.option(
    "--annotation",
    "annotation",
    multiple=True,
    metavar="NAME=PATH",
    help="Register a named annotation table (parquet/csv/tsv with a transcript_id "
    "and/or gene_id column), e.g. --annotation orthologs=orthologs.tsv. Repeatable.",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing dataset")
def create_dataset(
    path: Path,
    fasta: Path,
    species: str,
    gtf: Path | None,
    blocklist_fasta: tuple[Path, ...],
    fasta_key_regex: str,
    strip_version: bool,
    annotation: tuple[str, ...],
    overwrite: bool,
):
    if not species:
        raise ValueError("Species name is required.")
    annotations: dict[str, str | Path] = {}
    for spec in annotation:
        name, sep, table_path = spec.partition("=")
        if not sep or not name or not table_path:
            raise click.BadParameter(f"Expected NAME=PATH, got {spec!r}", param_hint="--annotation")
        annotations[name] = Path(table_path)
    Dataset.from_components(
        path,
        fasta,
        species=species,
        gtf_file=gtf,
        fasta_key_regex=fasta_key_regex,
        strip_version=strip_version,
        blocklist_fasta=list(blocklist_fasta) or None,
        annotations=annotations or None,
        overwrite=overwrite,
    )
