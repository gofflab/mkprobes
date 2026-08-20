# %%
import gzip
import json
import re
from functools import cache
from io import StringIO
from pathlib import Path
from typing import Any, Sequence, overload

import polars as pl
import pyfastx
import requests
from loguru import logger
from pydantic import BaseModel

from ..utils._alignment import bowtie_build, jellyfish


def get_ensembl(path: Path | str, id_: str, overwrite: bool = False):
    path = Path(path)
    path.mkdir(exist_ok=True, parents=True)
    if (p := (path / f"{id_}.json")).exists() and not overwrite:
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            logger.warning(f"Error decoding {p}. Deleting.")
            p.unlink(missing_ok=True)

    logger.info(f"Fetching {id_} on ensembl")
    res = requests.get(
        f"https://rest.ensembl.org/lookup/id/{id_}?content-type=application/json",
        timeout=30,
    )
    res.raise_for_status()
    p.write_text(json.dumps(j := res.json(), indent=2))
    return j


def _strip_version_suffix(token: str) -> str:
    head, dot, tail = token.rpartition(".")
    if dot and tail.isdigit():
        return head
    return token


# One GTF attribute token: `key "quoted value"` or `key unquoted_value`,
# terminated by `;` or end of string. Matched sequentially so semicolons
# inside quoted values are consumed as part of the value, never as a
# record separator.
_GTF_ATTR_TOKEN = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^;\s][^;]*?))\s*(?:;|$)')


def _parse_gtf_attributes(attributes: Sequence[str | None]) -> dict[str, list[str | None]]:
    """
    Tokenizes GTF attribute blobs into per-key columns.

    Scans each blob left-to-right with a quote-aware tokenizer, so values
    containing semicolons (e.g. `description "foo; bar"`) stay intact and
    text inside quotes can never be mistaken for a key. Keys are discovered
    across ALL rows (deterministic — no sampling). Repeated keys within a
    row (e.g. GENCODE multi-`tag` rows) resolve to the LAST occurrence.
    """
    columns: dict[str, list[str | None]] = {}
    for i, blob in enumerate(attributes):
        row: dict[str, str] = {}
        if blob:
            for m in _GTF_ATTR_TOKEN.finditer(blob):
                row[m.group(1)] = m.group(2) if m.group(2) is not None else m.group(3)
        for key in row:
            if key not in columns:
                columns[key] = [None] * i
        for key, values in columns.items():
            values.append(row.get(key))
    return columns


class MockGTF:
    """
    A placeholder for GTF data when a GTF file is not provided or cannot be loaded.

    This class is used internally by `ExternalData` when GTF information is unavailable.
    Any attempt to access attributes or items that would normally come from a parsed
    GTF file will result in a `NotImplementedError`.
    """

    def __getattr__(self, name: str) -> Any:
        raise NotImplementedError(f"ExternalData not created with GTF path. Cannot access {name}")

    def __getitem__(self, name: str | list[str]) -> Any:
        raise NotImplementedError("ExternalData not created with GTF path. Cannot index.")


class ExternalDataDefinition(BaseModel):
    """
    Defines the structure for storing external data source file names.

    This Pydantic model is used to serialize and deserialize the names of files
    associated with an `ExternalData` instance, such as the cache, GTF, FASTA,
    k-mer count file, and Bowtie2 index. It's primarily used when saving or
    loading dataset configurations that include `ExternalData`.

    Potential Pitfalls:
        - All file names are relative to a base path that is typically managed by
          the `Dataset` class. Ensure these names correctly point to existing files
          within that context when an `ExternalData` instance is created using
          this definition.
        - `cache_name` and `gtf_name` are optional. If `gtf_name` is not provided,
          the resulting `ExternalData` instance will use `MockGTF` for its GTF component,
          limiting its functionality.
        - `fasta_key_regex` and `strip_version` default to the historical behavior
          (first whitespace-delimited token, trailing `.N` version suffix stripped).
          Datasets ingested from de novo annotations set `strip_version=False` so
          GTF and FASTA identifiers match exactly (blanket `.N`-stripping would,
          e.g., collapse StringTie's `STRG.1.1`/`STRG.1.2` isoforms).
    """

    cache_name: str | None = None
    gtf_name: str | None = None
    fasta_name: str
    kmer18_name: str
    bowtie2_index_name: str
    fasta_key_regex: str = r"^(\S+)"
    strip_version: bool = True


class ExternalData:
    """
    Manages access to genomic data from GTF (Gene Transfer Format) and FASTA files.

    It provides methods to parse GTF files, retrieve gene and transcript information,
    and fetch sequences from FASTA files. The class uses caching for parsed GTF
    data to speed up subsequent initializations.

    Normal Use:
    -----------
    1. Initialize with paths to a cache location, a FASTA file, and a GTF file:
       ```python
       from pathlib import Path
       ext_data = ExternalData(
           cache="path/to/cache.parquet",  # Path for caching parsed GTF data
           fasta="path/to/genome.fasta",
           gtf_path="path/to/annotations.gtf"
       )
       ```
       The `cache` file does not need to exist beforehand; it will be created if a
       `gtf_path` is provided and the file is not already present or `regen_cache` is True.

    2. Access GTF data as a Polars DataFrame:
       ```python
       gtf_df = ext_data.gtf
       gene_info_df = ext_data.gene_info("MyGeneName")
       ```
    3. Retrieve sequences:
       ```python
       sequence = ext_data.get_seq("ENST00000123456")
       ```
    4. Convert between different ID types:
       ```python
       gene_id = ext_data.convert("MyGeneName", "gene_name", "gene_id")
       ```

    Potential Pitfalls:
    -------------------
    - **Caching**:
        - Parsed GTF data is cached to `cache` path. If the GTF file changes,
          `regen_cache=True` must be set during initialization to re-parse and
          update the cache. Otherwise, stale data might be used.
        - If `cache` path exists and `regen_cache=False` (default), `gtf_path` is
          ignored, and data is loaded directly from the cache.

    - **Missing GTF Path**:
        - If `gtf_path` is `None` AND no valid cache file exists at the `cache` path,
          `self.gtf` will be a `MockGTF` instance.
        - Most methods relying on GTF data (e.g., `gene_info`, `convert`, `get_transcripts`)
          will raise `NotImplementedError` or return default/input values.
        - This mode is intended for scenarios where only FASTA access is needed, or
          GTF data is managed externally.

    - **FASTA Key Function (`fasta_key_regex`)**:
        - This function normalizes sequence headers from the FASTA file to generate
          keys for sequence lookup (e.g., extracting transcript IDs).
        - The default function captures the first word.
        - If your FASTA headers have a different format, you MUST provide a custom
          `fasta_key_regex` to ensure IDs match those used/derived from the GTF.
          Mismatched keys will lead to `KeyError` or `ValueError` when fetching sequences.

    - **Required GTF Attributes**:
        - The `parse_gtf` method (and thus initialization with a `gtf_path`) requires
          `gene_id` and `transcript_id` to be present in the GTF attributes; a
          `ValueError` is raised otherwise. `gene_name`/`transcript_name` are
          optional — missing values are filled from a fallback chain
          (`gene_name` ← `Name` ← `gene` ← `gene_id`; `transcript_name` ← `transcript_id`).

    - **GTF Attribute Discovery**:
        - Attribute keys (e.g., `gene_name`, `transcript_type`) are discovered
          deterministically across ALL rows with a quote-aware tokenizer.
          Semicolons inside quoted values are safe; repeated keys within a row
          (e.g. GENCODE `tag`) resolve to the last occurrence.

    - **ID Versioning**:
        - When `strip_version=True` (default), a trailing `.N` version suffix
          (e.g., ".1", ".5") is stripped from FASTA keys, GTF
          `gene_id`/`transcript_id`, and lookup inputs (via `key_func`), so all
          three stay in agreement.
        - When `strip_version=False`, IDs pass through verbatim everywhere. Use
          this for de novo annotations whose IDs embed meaningful dots
          (StringTie `STRG.1.1` vs `STRG.1.2`, AUGUSTUS `g1.t1`).
    """

    def __init__(
        self,
        cache: Path | str | None = None,
        *,
        fasta: Path | str,
        gtf_path: Path | str | None = None,
        regen_cache: bool = False,
        fasta_key_regex: str = r"^(\S+)",
        bowtie2_index: str | None = None,
        kmer18: str | None = None,
        strip_version: bool = True,
    ) -> None:
        """
        Initializes an ExternalData object.

        See the class docstring for detailed information on parameters and behavior.

        Args:
            cache: Path to the Parquet file for caching parsed GTF data.
            fasta: Path to the FASTA file.
            gtf_path: Path to the GTF file. Optional if cache exists or only FASTA access is needed.
            regen_cache: If True, forces re-parsing of the GTF file and overwrites the cache.
            fasta_key_regex: Regex used to extract a lookup key from FASTA headers.
            bowtie2_index: Optional explicit name for the Bowtie2 index files (stem).
            kmer18: Optional explicit name for the 18-mer Jellyfish output file.
            strip_version: If True (historical default), trailing `.N` version
                suffixes are stripped from FASTA keys, GTF gene/transcript IDs,
                and lookup inputs alike. Set False for datasets whose IDs must
                match exactly (e.g. de novo annotations, where StringTie's
                `STRG.1.1`/`STRG.1.2` would otherwise collide).
        """
        self.fasta_path = Path(fasta).resolve()

        self.key_regex = fasta_key_regex
        self.strip_version = strip_version
        regex = re.compile(fasta_key_regex)

        def key_func(value: str) -> str:
            match = regex.match(value)
            token = (match.group(1) if match else value).strip()
            return _strip_version_suffix(token) if strip_version else token

        self.key_func = key_func

        try:
            self.fa = pyfastx.Fasta(Path(fasta).as_posix(), key_func=self.key_func)
        except Exception as e:
            raise Exception(
                f"Error reading FASTA file {fasta}. Ensure the file is valid and not empty."
            ) from e

        self._ts_gene_map: dict[str, str] | None = None

        if cache and Path(cache).exists() and not regen_cache:
            self.gtf: pl.DataFrame | MockGTF = pl.read_parquet(cache)
        else:
            if gtf_path is None:
                # logger.warning("GTF path not specified. Must be specified for reference species.")
                self.gtf = MockGTF()
            else:
                if not cache:
                    raise ValueError("Cache path must be specified if GTF path is provided.")
                self.gtf = self.parse_gtf(Path(gtf_path).resolve(), strip_version=strip_version)
                self.gtf.write_parquet(cache)

        self._override_bowtie2_index = bowtie2_index
        self._override_kmer18 = kmer18

    @classmethod
    def from_definition(cls, path: Path, definition: ExternalDataDefinition):
        """
        Creates an ExternalData instance from an ExternalDataDefinition.
        """
        return cls(
            cache=path / definition.cache_name if definition.cache_name else None,
            fasta=path / definition.fasta_name,
            gtf_path=path / definition.gtf_name if definition.gtf_name else None,
            fasta_key_regex=definition.fasta_key_regex,
            bowtie2_index=definition.bowtie2_index_name,
            kmer18=definition.kmer18_name,
            strip_version=definition.strip_version,
        )

    @property
    def bowtie2_index(self):
        """
        Gets the base path for the Bowtie2 index files.

        If `bowtie2_index` was provided during class initialization, that is used.
        Otherwise, it searches for index files (e.g., `{self.fasta_path}.1.bt2`) in the
        same directory as the FASTA file.

        Returns:
            Path: The base path (stem) of the Bowtie2 index files.

        Raises:
            FileNotFoundError: If the Bowtie2 index cannot be found.
        """
        if self._override_bowtie2_index:
            return self.fasta_path.parent / self._override_bowtie2_index

        bt = self.fasta_path.parent.glob(f"{self.fasta_path.stem}*.bt2")
        if not len(list(bt)):
            raise FileNotFoundError(
                f"Bowtie2 index not found for {self.fasta_path.stem}. Please build with `bowtie2-build {{fasta_path}} {{fasta file name}}` or call self.bowtie_build()."
            )

        return self.fasta_path.with_suffix("")

    def bowtie_build(self, overwrite: bool = False, interactive: bool = True):
        """
        Builds the Bowtie2 index for the FASTA file if it doesn't already exist.

        Checks for the existence of `fasta_stem.1.bt2`. If not found, it prompts
        the user before running `bowtie2-build` (unless `interactive=False`,
        for scripted flows like `mkprobes ingest`).

        Raises:
            FileNotFoundError: If `bowtie2-build` fails to create the index files.
        """
        if self.fasta_path.with_suffix(".1.bt2").exists() and not overwrite:
            return
        logger.info(f"Bowtie2 index not found for {self.fasta_path.stem}.")
        if interactive:
            input("\nPress Enter to start building...")
        bowtie_build(self.fasta_path, self.fasta_path.stem)
        if not self.fasta_path.with_suffix(".1.bt2").exists():
            raise FileNotFoundError(
                f"Bowtie2 index not found for {self.fasta_path.stem}. Bowtie2 build failed."
            )
        logger.info(f"Bowtie2 index successfully built for {self.fasta_path.stem}.")

    @property
    def kmer(self):
        """
        Gets the path to the k-mer count file (Jellyfish output).

        If `kmer18` was provided, that is used. Otherwise, it looks for
        a file named `{self.fasta_stem}.jf` in the same directory as the FASTA file.

        Returns:
            Path: The path to the k-mer file.

        Raises:
            FileNotFoundError: If the k-mer file cannot be found.
        """
        if self._override_kmer18:
            kmer18 = self.fasta_path.parent / self._override_kmer18
            if not kmer18.exists():
                raise FileNotFoundError(
                    f"Kmer file {kmer18} not found. Please run jellyfish or recreate the dataset."
                )
            return self.fasta_path.parent / kmer18

        if not self.fasta_path.with_suffix(".jf").exists():
            raise FileNotFoundError("Kmer file not found. Please run jellyfish.")
        return self.fasta_path.with_suffix(".jf")

    def run_jellyfish(
        self,
        kmer: int = 18,
        overwrite: bool = False,
        interactive: bool = True,
        hash_size: str = "10G",
    ):
        """
        Runs Jellyfish to count k-mers from the sequences in the FASTA file.

        Generates a `.jf` file (e.g., `fasta_stem.jf`). If the output file
        already exists and `overwrite` is False, the operation is skipped.
        Prompts the user before running if not overwriting (unless
        `interactive=False`, for scripted flows like `mkprobes ingest`).

        Args:
            kmer: The k-mer size to use for counting. Defaults to 18.
            overwrite: If True, run Jellyfish even if the output file exists.
                Defaults to False.
            interactive: If False, skip the confirmation prompt.

        Raises:
            FileNotFoundError: If Jellyfish fails to create the output file.
        """
        if self.fasta_path.with_suffix(".jf").exists() and not overwrite:
            logger.info(f"Jellyfish file {self.fasta_path.with_suffix('.jf')} already exists. Skipping.")
            return

        logger.info("Need to run jellyfish to get 18-mers in cDNA.")
        if not overwrite and interactive:
            input("Press Enter to start running...")

        jellyfish(
            [x.seq for x in self.fa],
            self.fasta_path.with_suffix(".jf"),
            kmer,
            hash_size=hash_size,
            minimum=10,
            counter=4,
        )
        if not self.fasta_path.with_suffix(".jf").exists():
            raise FileNotFoundError("cdna18.jf not found. Jellyfish run failed.")
        logger.info(f"Jellyfish file {self.fasta_path.with_suffix('.jf')} successfully created.")

    @cache
    def gene_info(self, gene: str) -> pl.DataFrame:
        """
        Retrieves all GTF entries for a given gene name.

        Args:
            gene: The gene name (e.g., "Actb").

        Returns:
            A Polars DataFrame containing rows from the GTF that match the gene name.
            Returns an empty DataFrame if the gene is not found or if GTF data is unavailable.
        """
        return self.gtf.filter(pl.col("gene_name") == gene)

    @cache
    def gene_to_eid(self, gene: str) -> pl.Series:
        """
        Converts a gene name to its corresponding gene ID(s) (e.g., Ensembl ID).

        Args:
            gene: The gene name.

        Returns:
            A Polars Series containing the gene ID(s) for the given gene name.

        Raises:
            ValueError: If the gene name is not found in the GTF data.
        """
        ret = self.gene_info(gene)
        if ret.is_empty():
            raise ValueError(f"Could not find {gene}")
        return ret[:, "gene_id"]

    @cache
    def ts_to_gene(self, ts: str) -> str:
        """
        Maps a transcript ID (version stripped) to its corresponding gene name.

        Uses an internal cached mapping. If the transcript ID is not found,
        it returns the input transcript ID.

        Args:
            ts: The transcript ID (e.g., "ENSMUST00000000001.4" or "ENSMUST00000000001").

        Returns:
            The gene name associated with the transcript ID, or the input `ts` if not found.
        """
        ts = self.key_func(ts)
        if self._ts_gene_map is None:
            self._ts_gene_map = {k: v for k, v in zip(self.gtf["transcript_id"], self.gtf["gene_name"])}
        return self._ts_gene_map.get(ts, ts)

    @cache
    def ts_to_tsname(self, eid: str) -> str | None:
        """
        Converts a transcript ID (version stripped) to its transcript name.

        If the transcript ID is found, returns the 'transcript_name' attribute.
        Otherwise, returns the input transcript ID.

        Args:
            eid: The transcript ID (e.g., "ENSMUST00000000001.4" or "ENSMUST00000000001").

        Returns:
            The transcript name, or the input `eid` if not found or if 'transcript_name'
            is missing.
        """
        eid = self.key_func(eid)
        try:
            res = self.gtf.filter(pl.col("transcript_id") == eid)["transcript_name"].first()
            return res if res else eid  # type: ignore
        except pl.exceptions.ComputeError:
            return eid

    @cache
    def eid_to_ts(self, eid: str) -> str:
        """
        Converts a gene ID (Ensembl ID, version stripped) to a transcript ID.

        Returns the first transcript ID associated with the given gene ID.

        Args:
            eid: The gene ID (e.g., "ENSMUSG00000000001.4" or "ENSMUSG00000000001").

        Returns:
            The first transcript ID found for that gene ID.

        Raises:
            Polars expression error or IndexError if the gene ID is not found or has no transcripts.
        """
        eid = self.key_func(eid)
        return self.gtf.filter(pl.col("gene_id") == eid)[0, "transcript_id"]

    def batch_convert(self, val: Sequence[str], src: str, dst: str) -> pl.DataFrame:
        """Batch convert attributes. See available attributes in `self.gtf.columns`.
        Will take the first value found for each attribute.
        !! Will skip non-existent values.

        Args:
            val: list of values to convert.
            src: Attribute to convert from (column name in GTF).
            dst: Attribute to convert to (column name in GTF).

        Returns:
            pl.DataFrame with two columns: `src` and `dst`, containing the mappings.

        Raises:
            ValueError: If none of the input values are found in the `src` column.
        """

        res = pl.DataFrame({src: val}).join(self.gtf.group_by(src).first(), on=src, how="inner")[[src, dst]]
        if not len(res):
            raise ValueError(f"Could not find {val} in {src}")
        if len(res) != len(val):
            logger.warning(
                f"Mapping not bijective. Some values are non-existent in the source column {len(res)} != {len(val)}"
            )
        return res

    def convert(self, val: str, src: str, dst: str) -> str:
        """
        Converts a single value from a source attribute to a destination attribute using GTF data.

        Args:
            val: The value to convert.
            src: The source attribute column name in the GTF data.
            dst: The destination attribute column name in the GTF data.

        Returns:
            The converted value from the `dst` attribute.

        Raises:
            ValueError: If the `val` is not found in the `src` column or if multiple
                matches are found (indicating non-uniqueness).
        """
        res = self.gtf.filter(pl.col(src) == val)[dst]
        if not len(res):
            raise ValueError(f"Could not find {val} in {src}")
        if len(res) > 1:
            raise ValueError(f"Found multiple {val} in {src}")
        return res[0]

    @cache
    def get_transcripts(self, gene: str | None = None, *, eid: str | None = None) -> pl.Series:
        """
        Retrieves transcript IDs for a given gene name or gene ID (Ensembl ID).

        Exactly one of `gene` or `eid` must be provided.

        Args:
            gene: The gene name.
            eid: The gene ID (Ensembl ID).

        Returns:
            A Polars Series containing all transcript IDs associated with the specified gene.
        """
        if gene is not None:
            return self.gtf.filter(pl.col("gene_name") == gene)["transcript_id"]
        return self.gtf.filter(pl.col("gene_id") == eid)["transcript_id"]

    @cache
    def get_seq(self, eid: str, convert: bool = True) -> str:
        """
        Retrieves a sequence from the FASTA file using an ID.

        The ID is typically a transcript ID. If `convert` is True and the ID contains
        a hyphen (suggesting it might be a transcript name), it first attempts to
        convert the transcript name to a transcript ID using `self.convert`.
        The ID (original or converted) is then version-stripped before FASTA lookup.
        If lookup with the version-stripped ID fails, it tries with the original ID.

        Args:
            eid: The identifier (transcript ID or transcript name) for the sequence.
            convert: If True and `eid` contains "-", attempt to convert it from
                'transcript_name' to 'transcript_id' first. Defaults to True.

        Returns:
            The sequence string.

        Raises:
            ValueError: If the ID cannot be found in the FASTA file after attempts,
                or if the sequence is empty.
            ValueError: If `convert` is True and the conversion from transcript name
                to ID fails.
        """
        if "-" in eid and convert:
            try:
                eid = self.convert(eid, "transcript_name", "transcript_id")
            except ValueError:
                logger.warning(
                    f"Could not convert {eid} from transcript_name to transcript_id. Trying to get the sequence directly."
                )

        eid_key = self.key_func(eid)
        keys = (eid_key,) if eid_key == eid else (eid_key, eid)
        for key in keys:
            try:
                res = self.fa[key].seq
                break
            except KeyError:
                continue
        else:
            raise ValueError(f"Could not find {eid} in fasta file.")

        if not res:
            raise ValueError(f"Could not find {eid}")
        return res

    def filter_gene(self, gene: str) -> pl.DataFrame:
        """
        Filters the GTF DataFrame for entries matching a specific gene name.

        This is a convenience method, equivalent to `self.gtf.filter(pl.col("gene_name") == gene)`.

        Args:
            gene: The gene name to filter by.

        Returns:
            A Polars DataFrame containing only rows related to the specified gene.
        """
        return self.gtf.filter(pl.col("gene_name") == gene)

    @overload
    def __getitem__(self, eid: str) -> pl.Series: ...

    @overload
    def __getitem__(self, eid: list[str]) -> pl.DataFrame: ...

    def __getitem__(self, eid: str | list[str]) -> pl.Series | pl.DataFrame:
        return self.gtf[eid]

    def filter(self, *args: Any, **kwargs: Any):
        return self.gtf.filter(*args, **kwargs)

    @staticmethod
    def parse_gtf(
        path: str | Path | StringIO,
        filters: Sequence[str] | None = ("transcript",),
        strip_version: bool = True,
    ) -> pl.DataFrame:
        """
        Parses a GTF file into a Polars DataFrame.

        Handles gzipped GTF files. Attribute keys are discovered deterministically
        across ALL rows and parsed with a quote-aware tokenizer (semicolons inside
        quoted values are safe; repeated keys such as GENCODE `tag` resolve to the
        last occurrence). 'gene_id' and 'transcript_id' are mandatory. If a
        `gene_name`/`transcript_name` attribute is absent, the columns are
        synthesized from the best available fallback so downstream name lookups
        never hit a missing column (`gene_name` ← `Name` ← `gene` ← `gene_id`;
        `transcript_name` ← `transcript_id`).

        Args:
            path: Path to the GTF file or an StringIO object containing GTF data.
            filters: A sequence of feature types (e.g., "transcript", "exon") to keep.
                If None, all features are kept. Defaults to ("transcript",).
            strip_version: If True (default), strips a trailing `.N` version suffix
                from `gene_id`/`transcript_id` — the same rule applied to FASTA
                keys, so the two stay in agreement. Set False to keep IDs verbatim
                (required for de novo annotations such as StringTie, where
                `STRG.1.1` and `STRG.1.2` are distinct isoforms).

        Returns:
            A Polars DataFrame representing the parsed GTF data.

        Raises:
            ValueError: If the file looks like GFF3 rather than GTF, if no rows
                match `filters`, or if 'gene_id'/'transcript_id' attributes are
                not found.
        """
        if not isinstance(path, StringIO) and Path(path).suffix == ".gz":
            path = StringIO(gzip.open(path, "rt").read())

        df = pl.read_csv(
            path,
            comment_prefix="#",
            separator="\t",
            has_header=False,
            new_columns=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
            schema_overrides=[
                pl.Utf8,
                pl.Utf8,
                pl.Utf8,
                pl.UInt32,
                pl.UInt32,
                pl.Utf8,
                pl.Utf8,
                pl.Utf8,
                pl.Utf8,
            ],
        )

        # GFF3 masquerading as GTF is the most common format mix-up for de novo
        # annotations. Detect it before feature-filtering, since GFF3 feature
        # names (mRNA, not transcript) would otherwise yield a confusing
        # "no rows" error.
        sample_attr = next((a for a in df["attribute"].head(50).to_list() if a), None)
        if (
            sample_attr
            and re.match(r"[\w.-]+=", sample_attr)
            and not re.search(r'\w+ +"', sample_attr)
        ):
            raise ValueError(
                "This looks like GFF3 (attributes use 'key=value'), not GTF. "
                "Convert it first with gffread: `gffread <in.gff3> -T -o <out.gtf>`. "
                "The `mkprobes ingest` command does this automatically."
            )

        available_features = df["feature"].unique().sort().to_list()
        df = df.filter(
            pl.col("feature").is_in(filters) if filters else pl.col("feature").is_not_null()
        )
        if df.is_empty():
            raise ValueError(
                f"No rows with feature in {tuple(filters or ())} found in GTF "
                f"(available features: {available_features}). "
                "Annotations lacking transcript rows can be normalized with "
                "`gffread <in> -T -o <out.gtf>`, or pass filters=None."
            )

        attr_columns = _parse_gtf_attributes(df["attribute"].to_list())
        logger.info(f"Found {len(attr_columns)} attributes in GTF file: {sorted(attr_columns)}")

        if "gene_id" not in attr_columns:
            raise ValueError("Gene ID not found in GTF file. Required attribute per GTF 2.0 spec.")
        if "transcript_id" not in attr_columns:
            raise ValueError("Transcript ID not found in GTF file. Required attribute per GTF 2.0 spec.")

        df = df.with_columns([
            pl.Series(name, values, dtype=pl.Utf8) for name, values in attr_columns.items()
        ])

        if strip_version:
            df = df.with_columns([
                pl.col("gene_id").str.replace(r"\.\d+$", "").alias("gene_id"),
                pl.col("transcript_id").str.replace(r"\.\d+$", "").alias("transcript_id"),
            ])

        gene_name_sources = [c for c in ("gene_name", "Name", "gene") if c in df.columns]
        if "gene_name" not in df.columns or df["gene_name"].null_count() > 0:
            logger.info(
                "gene_name incomplete or absent in GTF; filling from fallback chain: "
                f"{' -> '.join([*gene_name_sources, 'gene_id'])}"
            )
        df = df.with_columns(
            pl.coalesce([pl.col(c) for c in (*gene_name_sources, "gene_id")]).alias("gene_name")
        )
        if "transcript_name" in df.columns:
            df = df.with_columns(
                pl.coalesce([pl.col("transcript_name"), pl.col("transcript_id")]).alias(
                    "transcript_name"
                )
            )
        else:
            df = df.with_columns(pl.col("transcript_id").alias("transcript_name"))

        return pl.DataFrame(df)
