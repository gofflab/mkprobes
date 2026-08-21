import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from mkprobes.ext.dataset import (
    Dataset,
    parse_jellyfish,
)
from mkprobes.ext.external_data import (
    ExternalData,
    ExternalDataDefinition,
    MockGTF,
)


@pytest.fixture
def dummy_fasta_content() -> str:
    return ">seq1\nAGCTAGCT\n>seq2\nTCGATCGA\n"


@pytest.fixture
def dummy_jellyfish_content() -> str:
    return "AGCTAGCTAGCTAGCTAG 2\nTCGATCGATCGATCGATC 3\n"


@pytest.fixture
def dummy_appris_content() -> str:
    return "GENE_A\tENSG00A\tENST00A\tCCDS1\tPRINCIPAL:1\nGENE_B\tENSG00B\tENST00B\tCCDS2\tALTERNATIVE:2\n"


@pytest.fixture
def dummy_gtf_content() -> str:
    return """
# comment
chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1";
chr1\tHAVANA\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";
"""


@pytest.fixture
def mock_external_data(tmp_path: Path, dummy_fasta_content: str) -> ExternalData:
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(dummy_fasta_content)

    # Mock methods that might involve heavy computation or external calls
    mock_ed = MagicMock(spec=ExternalData)
    mock_ed.fasta_path = fasta_file
    mock_ed.fa = MagicMock()  # Simplified mock for pyfastx.Fasta
    mock_ed.gtf = MockGTF()
    mock_ed.bowtie2_index = tmp_path / "test_bowtie_index"
    mock_ed.kmer = tmp_path / "test_kmer.jf"
    mock_ed.bowtie_build = MagicMock()
    mock_ed.run_jellyfish = MagicMock()
    return mock_ed


class TestParseJellyfish:
    def test_parse_valid_file(self, tmp_path: Path, dummy_jellyfish_content: str):
        jf_file = tmp_path / "test.jf"
        jf_file.write_text(dummy_jellyfish_content)
        df = parse_jellyfish(jf_file)
        assert isinstance(df, pl.DataFrame)
        assert df.columns == ["kmer", "count"]
        assert len(df) == 2
        assert df["kmer"][0] == "AGCTAGCTAGCTAGCTAG"
        assert df["count"][1] == 3

    def test_parse_empty_file_returns_empty_frame(self, tmp_path: Path):
        # An existing-but-empty file means jellyfish ran and nothing cleared
        # the count threshold - legitimate for small transcriptomes.
        jf_file = tmp_path / "empty.jf"
        jf_file.write_text("")
        df = parse_jellyfish(jf_file)
        assert df.is_empty()
        assert df.columns == ["kmer", "count"]

    def test_parse_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_jellyfish("non_existent_file.jf")


class TestDataset:
    def test_init_minimal(self, tmp_path: Path, mock_external_data: MagicMock):
        dataset = Dataset(path=tmp_path, external_data=mock_external_data, species="test_species")
        assert dataset.path == tmp_path
        assert dataset.data == mock_external_data
        assert dataset.species == "test_species"
        assert dataset.kmer18 is None
        assert not dataset.kmerset
        assert dataset.trna_rna_kmers is None
        assert dataset.gencode == mock_external_data  # Backwards compatibility
        assert dataset.ensembl is None

    def test_init_with_kmers(
        self, tmp_path: Path, mock_external_data: MagicMock, dummy_jellyfish_content: str
    ):
        kmer18_file = tmp_path / "k18.jf"
        kmer18_file.write_text(dummy_jellyfish_content)
        trna_kmers_file = tmp_path / "trna.jf"
        trna_kmers_file.write_text("CCCC 1\n")

        dataset = Dataset(
            path=tmp_path,
            external_data=mock_external_data,
            kmer18_path=kmer18_file,
            trna_rna_kmers_path=trna_kmers_file,
        )
        assert dataset.kmer18 is not None
        assert len(dataset.kmer18) == 2
        assert "AGCTAGCTAGCTAGCTAG" in dataset.kmerset
        assert dataset.trna_rna_kmers is not None
        assert "CCCC" in dataset.trna_rna_kmers

    @patch("mkprobes.ext.dataset.ExternalData.from_definition")
    def test_from_folder(self, mock_ed_from_def: MagicMock, tmp_path: Path):
        dataset_path = tmp_path / "my_dataset_folder"
        dataset_path.mkdir()

        mock_ed_instance = MagicMock(spec=ExternalData)
        mock_ed_instance.kmer = "test.jf"  # kmer path for constructor
        mock_ed_from_def.return_value = mock_ed_instance

        definition_content = {
            "species": "folder_species",
            "external_data": {
                "default": {
                    "fasta_name": "test.fa",
                    "bowtie2_index_name": "test.bt2",
                    "kmer18_name": "test.jf",
                }
            },
        }
        (dataset_path / "dataset.json").write_text(json.dumps(definition_content))
        (dataset_path / "test.jf").write_text("AGCTAGCTAGCTAGCTAG 2\nTCGATCGATCGATCGATC 3\n")

        dataset = Dataset.from_folder(dataset_path)

        expected_def = ExternalDataDefinition(**definition_content["external_data"]["default"])
        mock_ed_from_def.assert_called_once_with(dataset_path, expected_def)
        assert dataset.species == "folder_species"
        assert dataset.data == mock_ed_instance
        assert dataset.path == dataset_path

    def test_from_folder_no_json(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="does not exist. Please create a dataset first."):
            Dataset.from_folder(tmp_path / "non_existent_dataset_dir")

    def test_check_kmers(self, tmp_path: Path, mock_external_data: MagicMock):
        dataset = Dataset(path=tmp_path, external_data=mock_external_data)

        # Case 1: trna_rna_kmers is None
        assert not dataset.check_kmers("AGCT")

        # Case 2: trna_rna_kmers is set
        dataset.trna_rna_kmers = {"AGCT", "GGGG"}  # Assuming kmer length 4 for simplicity here
        with patch("mkprobes.ext.dataset.kmers", return_value=["AGCT", "GCTA", "CTAG"]):
            assert dataset.check_kmers("AGCTAG")
        with patch("mkprobes.ext.dataset.kmers", return_value=["TTTT", "AAAA"]):
            assert not dataset.check_kmers("TTTTAA")

    def test_check_kmers_uses_blacklist_k(self, tmp_path: Path, mock_external_data: MagicMock):
        dataset = Dataset(path=tmp_path, external_data=mock_external_data)
        dataset.trna_rna_kmers = {"A" * 15}
        assert dataset.check_kmers("A" * 15)

    def test_appris_property_not_implemented(self, tmp_path: Path, mock_external_data: MagicMock):
        dataset = Dataset(path=tmp_path, external_data=mock_external_data)
        with pytest.raises(NotImplementedError):
            _ = dataset.appris


class TestDatasetFromComponents:
    @patch("mkprobes.ext.dataset.ExternalData")
    def test_from_components_success(
        self,
        mock_external_data_cls: MagicMock,
        tmp_path: Path,
        dummy_fasta_content: str,  # Make sure this fixture provides non-empty FASTA content
        dummy_jellyfish_content: str,  # Make sure this fixture provides non-empty k-mer content
    ):
        # 1. Setup paths and source FASTA file
        source_fasta_dir = tmp_path / "source_files"
        source_fasta_dir.mkdir()
        source_fasta_file = source_fasta_dir / "input.fasta"
        source_fasta_file.write_text(dummy_fasta_content)
        assert source_fasta_file.read_text() != "", "Dummy FASTA content should not be empty"

        dataset_path = tmp_path / "my_dataset"
        species = "test_organism"

        # 2. Configure the mock ExternalData instance
        mock_ed_instance = MagicMock(spec=ExternalData)

        # ExternalData is initialized with the original fasta_file path
        mock_ed_instance.fasta_path = source_fasta_file

        # Mock properties that return paths/names based on the original fasta_file
        # .kmer property will return source_fasta_file.with_suffix(".jf")
        expected_kmer_file_path = source_fasta_file.with_suffix(".jf")
        mock_ed_instance.kmer = expected_kmer_file_path

        # .bowtie2_index property will return source_fasta_file.with_suffix("") (the stem)
        expected_bowtie_index_stem_path = source_fasta_file.with_suffix("")
        mock_ed_instance.bowtie2_index = expected_bowtie_index_stem_path

        # Mock methods
        mock_ed_instance.bowtie_build = MagicMock()
        mock_ed_instance.run_jellyfish = MagicMock()

        mock_external_data_cls.return_value = mock_ed_instance

        # 3. Create the dummy k-mer file that Dataset.__init__ will parse.
        # This file is located relative to the *source* FASTA, as per external_data.kmer.
        expected_kmer_file_path.write_text(dummy_jellyfish_content)
        assert expected_kmer_file_path.read_text() != "", "Dummy k-mer content should not be empty"

        # 4. Call Dataset.from_components
        dataset = Dataset.from_components(
            path=dataset_path, fasta_file=source_fasta_file, species=species, overwrite=True
        )

        # 5. Assertions
        # Assertions for ExternalData method calls
        mock_external_data_cls.assert_called_once_with(
            cache=dataset_path / source_fasta_file.with_suffix(".parquet").name,  # Cache is in dataset_path
            fasta=dataset_path / source_fasta_file.name,  # Initialized with original fasta_file
            gtf_path=None,
            regen_cache=True,
            fasta_key_regex=r"^(\S+)",
            strip_version=True,
        )
        mock_ed_instance.bowtie_build.assert_called_once()
        mock_ed_instance.run_jellyfish.assert_called_once()

        # Assertion for file copying
        # shutil.copy is called if source and target are different.
        # Target is dataset_path / source_fasta_file.name
        expected_copied_fasta_path = dataset_path / source_fasta_file.name

        # Assertions for dataset.json
        dataset_json_path = dataset_path / "dataset.json"
        assert dataset_json_path.exists(), "dataset.json was not created"

        json_data = json.loads(dataset_json_path.read_text())

        assert json_data["species"] == species
        assert "default" in json_data["external_data"]
        default_ext_data_def = json_data["external_data"]["default"]

        assert default_ext_data_def["fasta_name"] == source_fasta_file.name
        assert default_ext_data_def["bowtie2_index_name"] == expected_bowtie_index_stem_path.name
        assert default_ext_data_def["kmer18_name"] == expected_kmer_file_path.name

        # Assertions for the returned Dataset object
        assert isinstance(dataset, Dataset)
        assert dataset.path == dataset_path
        assert dataset.species == species
        assert dataset.data == mock_ed_instance
        assert dataset.kmer18 is not None, "Dataset kmer18 data should be loaded"
        assert len(dataset.kmer18) > 0, "Dataset kmer18 data should not be empty"
        assert "AGCTAGCTAGCTAGCTAG" in dataset.kmerset  # Based on dummy_jellyfish_content

        # Check that the FASTA file was indeed copied (even though shutil.copy is mocked,
        # the path used by ExternalData for its cache, etc., implies this structure)
        assert expected_copied_fasta_path.exists(), "FASTA file was not copied to dataset directory"
        assert expected_copied_fasta_path.read_text() == dummy_fasta_content


class TestLoadDataset:
    """Resolution matrix for load_dataset()."""

    def test_resolves_dataset_json_to_generic(self, tmp_path: Path):
        from mkprobes.ext.dataset import load_dataset

        ds_dir = tmp_path / "axolotl"
        ds_dir.mkdir()
        (ds_dir / "dataset.json").write_text(
            json.dumps({
                "species": "axolotl",
                "external_data": {
                    "default": {
                        "fasta_name": "t.fa",
                        "bowtie2_index_name": "t",
                        "kmer18_name": "t.jf",
                    }
                },
            })
        )
        (ds_dir / "t.jf").write_text("AGCTAGCTAGCTAGCTAG 2\n")
        with patch("mkprobes.ext.dataset.ExternalData.from_definition") as mock_from_def:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.kmer = "t.jf"
            mock_from_def.return_value = mock_ed
            ds = load_dataset(ds_dir)
        assert type(ds) is Dataset
        assert ds.species == "axolotl"

    def test_resolves_human_mouse_to_reference(self, tmp_path: Path):
        from mkprobes.ext import dataset as dataset_mod

        for name in ("human", "mouse"):
            ref_dir = tmp_path / name
            ref_dir.mkdir()
            with patch.object(dataset_mod, "ReferenceDataset") as mock_ref:
                dataset_mod.load_dataset(ref_dir)
                mock_ref.assert_called_once_with(ref_dir)

    def test_dataset_json_wins_over_folder_name(self, tmp_path: Path):
        # A folder literally named "mouse" that carries dataset.json is a
        # custom dataset, not a reference dataset.
        from mkprobes.ext.dataset import load_dataset

        ds_dir = tmp_path / "mouse"
        ds_dir.mkdir()
        (ds_dir / "dataset.json").write_text(
            json.dumps({
                "species": "mouse_custom",
                "external_data": {
                    "default": {
                        "fasta_name": "t.fa",
                        "bowtie2_index_name": "t",
                        "kmer18_name": "t.jf",
                    }
                },
            })
        )
        (ds_dir / "t.jf").write_text("AGCTAGCTAGCTAGCTAG 2\n")
        with patch("mkprobes.ext.dataset.ExternalData.from_definition") as mock_from_def:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.kmer = "t.jf"
            mock_from_def.return_value = mock_ed
            ds = load_dataset(ds_dir)
        assert type(ds) is Dataset

    def test_unrecognized_path_fails_with_guidance(self, tmp_path: Path):
        from mkprobes.ext.dataset import load_dataset

        with pytest.raises(FileNotFoundError, match="mkprobes ingest.*create-dataset"):
            load_dataset(tmp_path / "squid")

    def test_dataset_json_over_reference_build_is_ambiguous(self, tmp_path: Path):
        # The silent-downgrade case: a dataset.json dropped into a downloaded
        # reference build used to win, quietly losing pseudogene screening.
        from mkprobes.ext.dataset import load_dataset

        ref_dir = tmp_path / "mouse"
        ref_dir.mkdir()
        (ref_dir / "gencode.gtf.gz").write_bytes(b"")
        (ref_dir / "dataset.json").write_text(json.dumps({"species": "mouse", "external_data": {}}))

        with pytest.raises(ValueError, match="ambiguous"):
            load_dataset(ref_dir)

    def test_custom_reference_species_dataset_warns(self, tmp_path: Path):
        # A custom mouse dataset is still allowed, but must not be silent about
        # using the reduced path.
        from loguru import logger

        from mkprobes.ext.dataset import load_dataset

        ds_dir = tmp_path / "my_mouse"
        ds_dir.mkdir()
        (ds_dir / "dataset.json").write_text(
            json.dumps({
                "species": "mouse",
                "external_data": {
                    "default": {
                        "fasta_name": "t.fa",
                        "bowtie2_index_name": "t",
                        "kmer18_name": "t.jf",
                    }
                },
            })
        )
        (ds_dir / "t.jf").write_text("AGCTAGCTAGCTAGCTAG 2\n")

        messages: list[str] = []
        sink_id = logger.add(messages.append, level="WARNING")
        try:
            with patch("mkprobes.ext.dataset.ExternalData.from_definition") as mock_from_def:
                mock_ed = MagicMock(spec=ExternalData)
                mock_ed.kmer = "t.jf"
                mock_from_def.return_value = mock_ed
                ds = load_dataset(ds_dir)
        finally:
            logger.remove(sink_id)

        assert type(ds) is Dataset
        assert any("no pseudogene screening" in m for m in messages), messages

    def test_reference_build_in_misnamed_folder_explains_rename(self, tmp_path: Path):
        from mkprobes.ext.dataset import load_dataset

        ref_dir = tmp_path / "mm39"
        ref_dir.mkdir()
        (ref_dir / "gencode.gtf.gz").write_bytes(b"")

        with pytest.raises(ValueError, match="must be exactly `human` or `mouse`"):
            load_dataset(ref_dir)


class TestFromComponentsGtfAndBlocklist:
    def _run(
        self,
        tmp_path: Path,
        dummy_fasta_content: str,
        dummy_jellyfish_content: str,
        dummy_gtf_content: str,
        **kwargs,
    ):
        src = tmp_path / "src"
        src.mkdir(exist_ok=True)
        fasta = src / "txome.fasta"
        fasta.write_text(dummy_fasta_content)
        gtf = src / "anno.gtf"
        gtf.write_text(dummy_gtf_content)

        ds_path = tmp_path / "ds"

        with patch("mkprobes.ext.dataset.ExternalData") as mock_cls:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.fasta_path = ds_path / "txome.fasta"
            kmer_path = ds_path / "txome.jf"
            mock_ed.kmer = kmer_path
            mock_ed.bowtie2_index = ds_path / "txome"
            mock_cls.return_value = mock_ed

            with patch("mkprobes.ext.dataset.jellyfish") as mock_jf:

                def fake_jellyfish(seqs, out, k, **kw):
                    Path(out).write_text("CCCCCCCCCCCCCCC 5\n")

                mock_jf.side_effect = fake_jellyfish
                ds_path.mkdir(exist_ok=True, parents=True)
                kmer_path.write_text(dummy_jellyfish_content)
                ds = Dataset.from_components(
                    ds_path, fasta, species="axolotl", gtf_file=gtf, **kwargs
                )
        return ds, ds_path, mock_cls, mock_jf

    def test_gtf_wired_into_definition(
        self,
        tmp_path: Path,
        dummy_fasta_content: str,
        dummy_jellyfish_content: str,
        dummy_gtf_content: str,
    ):
        ds, ds_path, mock_cls, _ = self._run(
            tmp_path,
            dummy_fasta_content,
            dummy_jellyfish_content,
            dummy_gtf_content,
            strip_version=False,
        )
        assert (ds_path / "anno.gtf").exists()
        d = json.loads((ds_path / "dataset.json").read_text())
        default = d["external_data"]["default"]
        assert default["gtf_name"] == "anno.gtf"
        assert default["cache_name"] == "txome.parquet"
        assert default["strip_version"] is False
        # ExternalData constructed with the copied GTF and strip flag
        _, called_kwargs = mock_cls.call_args
        assert called_kwargs["gtf_path"] == ds_path / "anno.gtf"
        assert called_kwargs["strip_version"] is False

    def test_blocklist_built_and_threaded(
        self,
        tmp_path: Path,
        dummy_fasta_content: str,
        dummy_jellyfish_content: str,
        dummy_gtf_content: str,
    ):
        rrna = tmp_path / "rrna.fa"
        rrna.write_text(">rrna1\nCCCCCCCCCCCCCCCCCC\n")
        ds, ds_path, _, mock_jf = self._run(
            tmp_path,
            dummy_fasta_content,
            dummy_jellyfish_content,
            dummy_gtf_content,
            blocklist_fasta=[rrna],
        )
        mock_jf.assert_called_once()
        (seqs, out, k) = mock_jf.call_args[0]
        assert k == 15
        assert Path(out).name == "blocklist15.jf"
        assert seqs == ["CCCCCCCCCCCCCCCCCC"]

        d = json.loads((ds_path / "dataset.json").read_text())
        assert d["blocklist_kmer_name"] == "blocklist15.jf"
        # threaded into the constructed Dataset: check_kmers is live
        assert ds.trna_rna_kmers == {"CCCCCCCCCCCCCCC"}
        assert ds.check_kmers("AAACCCCCCCCCCCCCCCAAA")
        assert not ds.check_kmers("ATATATATATATATATATATAT")

    def test_empty_blocklist_fasta_raises(
        self,
        tmp_path: Path,
        dummy_fasta_content: str,
        dummy_jellyfish_content: str,
        dummy_gtf_content: str,
    ):
        empty = tmp_path / "empty.fa"
        empty.write_text("")
        with pytest.raises(ValueError, match="No sequences found in blocklist"):
            self._run(
                tmp_path,
                dummy_fasta_content,
                dummy_jellyfish_content,
                dummy_gtf_content,
                blocklist_fasta=[empty],
            )


class TestFromFolderBlocklist:
    def test_blocklist_threaded_from_definition(self, tmp_path: Path):
        ds_dir = tmp_path / "ds"
        ds_dir.mkdir()
        (ds_dir / "dataset.json").write_text(
            json.dumps({
                "species": "axolotl",
                "external_data": {
                    "default": {
                        "fasta_name": "t.fa",
                        "bowtie2_index_name": "t",
                        "kmer18_name": "t.jf",
                    }
                },
                "blocklist_kmer_name": "blocklist15.jf",
            })
        )
        (ds_dir / "t.jf").write_text("AGCTAGCTAGCTAGCTAG 2\n")
        (ds_dir / "blocklist15.jf").write_text("GGGGGGGGGGGGGGG 7\n")
        with patch("mkprobes.ext.dataset.ExternalData.from_definition") as mock_from_def:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.kmer = "t.jf"
            mock_from_def.return_value = mock_ed
            ds = Dataset.from_folder(ds_dir)
        assert ds.trna_rna_kmers == {"GGGGGGGGGGGGGGG"}
        assert ds.check_kmers("AAGGGGGGGGGGGGGGGAA")

    def test_old_format_dataset_json_loads_without_blocklist(self, tmp_path: Path):
        ds_dir = tmp_path / "ds"
        ds_dir.mkdir()
        (ds_dir / "dataset.json").write_text(
            json.dumps({
                "species": "squid",
                "external_data": {
                    "default": {
                        "fasta_name": "t.fa",
                        "bowtie2_index_name": "t",
                        "kmer18_name": "t.jf",
                    }
                },
            })
        )
        (ds_dir / "t.jf").write_text("AGCTAGCTAGCTAGCTAG 2\n")
        with patch("mkprobes.ext.dataset.ExternalData.from_definition") as mock_from_def:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.kmer = "t.jf"
            mock_from_def.return_value = mock_ed
            ds = Dataset.from_folder(ds_dir)
        assert ds.trna_rna_kmers is None


class TestGenomeAndAnnotations:
    def test_annotation_tables_copied_validated_loadable(
        self,
        tmp_path: Path,
        dummy_fasta_content: str,
        dummy_jellyfish_content: str,
    ):
        src = tmp_path / "src"
        src.mkdir()
        fasta = src / "txome.fasta"
        fasta.write_text(dummy_fasta_content)
        ortho = src / "orthologs.tsv"
        ortho.write_text("transcript_id\thuman_symbol\nseq1\tSHANK3\nseq2\tGRIN1\n")
        genome = src / "genome.fa"
        genome.write_text(">scaffold_1\nACGTACGTACGT\n")

        ds_path = tmp_path / "ds"
        with patch("mkprobes.ext.dataset.ExternalData") as mock_cls:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.fasta_path = ds_path / "txome.fasta"
            kmer_path = ds_path / "txome.jf"
            mock_ed.kmer = kmer_path
            mock_ed.bowtie2_index = ds_path / "txome"
            mock_cls.return_value = mock_ed
            ds_path.mkdir(parents=True)
            kmer_path.write_text(dummy_jellyfish_content)
            ds = Dataset.from_components(
                ds_path,
                fasta,
                species="axolotl",
                genome_fasta=genome,
                annotations={"orthologs": ortho},
            )

        d = json.loads((ds_path / "dataset.json").read_text())
        assert d["genome_fasta_name"] == "genome.fa"
        assert d["annotations"] == {"orthologs": "orthologs.tsv"}
        assert (ds_path / "genome.fa").exists()
        assert (ds_path / "orthologs.tsv").exists()

        assert ds.genome_fasta_path == ds_path / "genome.fa"
        table = ds.annotation("orthologs")
        assert table.filter(pl.col("transcript_id") == "seq1")[0, "human_symbol"] == "SHANK3"
        with pytest.raises(KeyError, match="No annotation table named 'fpkm'"):
            ds.annotation("fpkm")

    def test_annotation_without_join_column_rejected(
        self,
        tmp_path: Path,
        dummy_fasta_content: str,
        dummy_jellyfish_content: str,
    ):
        src = tmp_path / "src"
        src.mkdir()
        fasta = src / "txome.fasta"
        fasta.write_text(dummy_fasta_content)
        bad = src / "bad.csv"
        bad.write_text("some_column,other\na,b\n")

        ds_path = tmp_path / "ds"
        with patch("mkprobes.ext.dataset.ExternalData") as mock_cls:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.fasta_path = ds_path / "txome.fasta"
            kmer_path = ds_path / "txome.jf"
            mock_ed.kmer = kmer_path
            mock_ed.bowtie2_index = ds_path / "txome"
            mock_cls.return_value = mock_ed
            ds_path.mkdir(parents=True)
            kmer_path.write_text(dummy_jellyfish_content)
            with pytest.raises(ValueError, match="no join column"):
                Dataset.from_components(
                    ds_path, fasta, species="axolotl", annotations={"bad": bad}
                )

    def test_from_folder_threads_genome_and_annotations(self, tmp_path: Path):
        ds_dir = tmp_path / "ds"
        ds_dir.mkdir()
        (ds_dir / "dataset.json").write_text(
            json.dumps({
                "species": "axolotl",
                "external_data": {
                    "default": {
                        "fasta_name": "t.fa",
                        "bowtie2_index_name": "t",
                        "kmer18_name": "t.jf",
                    }
                },
                "genome_fasta_name": "genome.fa",
                "annotations": {"fpkm": "fpkm.csv"},
            })
        )
        (ds_dir / "t.jf").write_text("AGCTAGCTAGCTAGCTAG 2\n")
        (ds_dir / "fpkm.csv").write_text("gene_id,fpkm\ng1,12.5\n")
        with patch("mkprobes.ext.dataset.ExternalData.from_definition") as mock_from_def:
            mock_ed = MagicMock(spec=ExternalData)
            mock_ed.kmer = "t.jf"
            mock_from_def.return_value = mock_ed
            ds = Dataset.from_folder(ds_dir)
        assert ds.genome_fasta_path == ds_dir / "genome.fa"
        assert ds.annotation("fpkm")[0, "fpkm"] == 12.5
