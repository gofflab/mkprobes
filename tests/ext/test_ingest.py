import gzip
import json
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest
from conftest import flatten_cli_output
import yaml
from click.testing import CliRunner

from mkprobes.ext.ingest import (
    IngestValidationReport,
    check_external_tools,
    detect_annotation_format,
    extract_biotype_blocklist,
    ingest,
    validate_gtf,
    write_intake_manifest,
)

DATA = Path(__file__).parent.parent / "data"
MINI_GENOME = DATA / "mini_genome.fa"

EXTERNAL_TOOLS_PRESENT = all(shutil.which(t) for t in ("gffread", "bowtie2-build", "jellyfish"))


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    # pyfastx writes .fxi index files next to the FASTA; keep fixtures pristine.
    dest = tmp_path / "mini_genome.fa"
    shutil.copy(MINI_GENOME, dest)
    return dest


class TestDetectFormat:
    def test_gtf(self):
        assert detect_annotation_format(DATA / "augustus.gtf") == "gtf"

    def test_gff3(self, tmp_path: Path):
        p = tmp_path / "a.gff3"
        p.write_text("scaffold_1\tmaker\tmRNA\t100\t500\t.\t+\t.\tID=rna1;Parent=gene1\n")
        assert detect_annotation_format(p) == "gff3"

    def test_gz(self, tmp_path: Path):
        p = tmp_path / "a.gtf.gz"
        with gzip.open(p, "wt") as fh:
            fh.write('scaffold_1\tX\ttranscript\t1\t99\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
        assert detect_annotation_format(p) == "gtf"


class TestCheckExternalTools:
    def test_aggregates_missing(self):
        with patch("mkprobes.ext.ingest.shutil.which", return_value=None):
            with pytest.raises(RuntimeError) as exc:
                check_external_tools()
        msg = str(exc.value)
        assert "gffread" in msg and "bowtie2-build" in msg and "jellyfish" in msg
        assert "bioconda" in msg

    def test_versions_captured_when_present(self):
        if not shutil.which("jellyfish"):
            pytest.skip("jellyfish not on PATH")
        versions = check_external_tools(required=("jellyfish",))
        assert "jellyfish" in versions


class TestValidateGtf:
    @pytest.mark.parametrize("style", ["augustus", "stringtie", "maker", "ncbi_style"])
    def test_mini_species_all_styles_pass(self, style: str, genome: Path):
        report = validate_gtf(DATA / f"{style}.gtf", genome)
        assert report.ok, [i.message for i in report.errors]
        assert report.n_genes == 6
        assert report.n_transcripts == 7
        assert report.feature_counts["exon"] == 15
        assert report.transcripts_per_gene_max == 2

    def test_gene_name_sources(self, genome: Path):
        assert validate_gtf(DATA / "augustus.gtf", genome).gene_name_source == "gene_id (fallback)"
        assert validate_gtf(DATA / "ncbi_style.gtf", genome).gene_name_source == "gene"

    def test_seqname_mismatch_all_missing_is_error(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "bad.gtf"
        gtf.write_text(
            'chr1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tX\texon\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        report = validate_gtf(gtf, genome)
        assert not report.ok
        assert any(i.code == "SEQNAME_MISMATCH" and i.severity == "error" for i in report.issues)

    def test_seqname_mismatch_partial_is_warning(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "bad.gtf"
        gtf.write_text(
            'scaffold_1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\texon\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chrUn\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n'
            'chrUn\tX\texon\t10\t500\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n'
        )
        report = validate_gtf(gtf, genome)
        assert any(i.code == "SEQNAME_MISMATCH" and i.severity == "warning" for i in report.issues)

    def test_duplicate_transcript_ids(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "dup.gtf"
        gtf.write_text(
            'scaffold_1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\texon\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\ttranscript\t600\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\texon\t600\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        report = validate_gtf(gtf, genome)
        assert any(i.code == "DUPLICATE_TRANSCRIPT_ID" for i in report.errors)

    def test_no_exon_rows(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "noexon.gtf"
        gtf.write_text('scaffold_1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
        report = validate_gtf(gtf, genome)
        assert any(i.code == "NO_EXON_ROWS" for i in report.errors)

    def test_transcript_without_exons_is_warning(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "orphan.gtf"
        gtf.write_text(
            'scaffold_1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\texon\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\ttranscript\t600\t900\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n'
        )
        report = validate_gtf(gtf, genome)
        assert report.ok
        assert any(i.code == "TRANSCRIPT_WITHOUT_EXONS" for i in report.warnings)

    def test_forbidden_id_chars(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "pipe.gtf"
        gtf.write_text(
            'scaffold_1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t|1";\n'
            'scaffold_1\tX\texon\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t|1";\n'
        )
        report = validate_gtf(gtf, genome)
        assert any(i.code == "ID_FORBIDDEN_CHARS" for i in report.errors)

    def test_inverted_coordinates(self, tmp_path: Path, genome: Path):
        gtf = tmp_path / "coords.gtf"
        gtf.write_text(
            'scaffold_1\tX\ttranscript\t500\t10\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'scaffold_1\tX\texon\t500\t10\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        report = validate_gtf(gtf, genome)
        assert any(i.code == "COORDINATES_INVERTED" for i in report.errors)

    def test_gff3_input_reports_parse_failed(self, tmp_path: Path, genome: Path):
        gff3 = tmp_path / "a.gff3"
        gff3.write_text("scaffold_1\tmaker\tmRNA\t100\t500\t.\t+\t.\tID=rna1;Parent=gene1\n")
        report = validate_gtf(gff3, genome)
        assert any(i.code == "PARSE_FAILED" for i in report.errors)


class TestBiotypeBlocklist:
    def test_extracts_rrna(self, tmp_path: Path):
        from mkprobes.ext.external_data import ExternalData

        gtf_df = ExternalData.parse_gtf(DATA / "ncbi_style.gtf", filters=None, strip_version=False)
        tx_fasta = tmp_path / "transcripts.fasta"
        ids = gtf_df.filter(gtf_df["gene_biotype"] == "rRNA")["transcript_id"].drop_nulls().unique().to_list()
        assert len(ids) == 1
        rrna_id = ids[0]
        tx_fasta.write_text(f">{rrna_id}\n" + "ACGT" * 30 + "\n>rna-XM_001001.1\n" + "GGCC" * 30 + "\n")

        out = extract_biotype_blocklist(gtf_df, tx_fasta, ["rRNA"], tmp_path / "biotype_blocklist.fasta")
        assert out is not None
        content = out.read_text()
        assert rrna_id in content
        assert "rna-XM_001001.1" not in content

    def test_no_biotype_column_returns_none(self, tmp_path: Path):
        from mkprobes.ext.external_data import ExternalData

        gtf_df = ExternalData.parse_gtf(DATA / "augustus.gtf", filters=None, strip_version=False)
        tx_fasta = tmp_path / "t.fasta"
        tx_fasta.write_text(">g1.t1\nACGTACGT\n")
        assert (
            extract_biotype_blocklist(gtf_df, tx_fasta, ["rRNA"], tmp_path / "out.fasta") is None
        )


class TestManifest:
    def test_manifest_written_and_loadable(self, tmp_path: Path, genome: Path):
        report = validate_gtf(DATA / "augustus.gtf", genome)
        gtf_copy = tmp_path / "annotation.gtf"
        shutil.copy(DATA / "augustus.gtf", gtf_copy)
        out = write_intake_manifest(
            tmp_path,
            species="mini",
            genome=genome,
            annotation=gtf_copy,
            annotation_format="gtf",
            extract_mode="transcripts",
            fasta_key_regex=r"^(\S+)",
            strip_version=False,
            tool_versions={"gffread": "0.12.7"},
            report=report,
            blocklist_files=[],
            annotation_tables={},
            argv=["mkprobes", "ingest", "x"],
        )
        m = yaml.safe_load(out.read_text())
        assert m["manifest_version"] == 2
        assert m["species"]["display_name"] == "mini"
        assert len(m["inputs"]["genome_sha256"]) == 64
        assert m["quality_control"]["n_transcripts"] == 7
        assert m["quality_control"]["passed_for_probe_generation"] is True
        assert m["processing"]["strip_version"] is False


class TestIngestCli:
    def test_reserved_species_rejected(self, tmp_path: Path, genome: Path):
        runner = CliRunner()
        res = runner.invoke(
            ingest,
            [
                str(tmp_path / "ds"),
                "--genome",
                str(genome),
                "--gtf",
                str(DATA / "augustus.gtf"),
                "--species",
                "mouse",
            ],
        )
        assert res.exit_code != 0
        assert "reserved" in res.output

    def test_validate_only_writes_report(self, tmp_path: Path, genome: Path):
        runner = CliRunner()
        with patch("mkprobes.ext.ingest.check_external_tools", return_value={}):
            res = runner.invoke(
                ingest,
                [
                    str(tmp_path / "ds"),
                    "--genome",
                    str(genome),
                    "--gtf",
                    str(DATA / "augustus.gtf"),
                    "--species",
                    "mini",
                    "--validate-only",
                ],
            )
        assert res.exit_code == 0, res.output
        report_file = tmp_path / "ds" / "validation_report.json"
        assert report_file.exists()
        report = IngestValidationReport.model_validate_json(report_file.read_text())
        assert report.ok
        assert (tmp_path / "ds" / "annotation.gtf").exists()

    def test_validate_only_fails_on_bad_gtf(self, tmp_path: Path, genome: Path):
        bad = tmp_path / "bad.gtf"
        bad.write_text(
            'chr1\tX\ttranscript\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tX\texon\t10\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        runner = CliRunner()
        with patch("mkprobes.ext.ingest.check_external_tools", return_value={}):
            res = runner.invoke(
                ingest,
                [
                    str(tmp_path / "ds"),
                    "--genome",
                    str(genome),
                    "--gtf",
                    str(bad),
                    "--species",
                    "mini",
                    "--validate-only",
                ],
            )
        assert res.exit_code != 0
        assert "Validation failed" in res.output


@pytest.mark.external
@pytest.mark.skipif(not EXTERNAL_TOOLS_PRESENT, reason="gffread/bowtie2/jellyfish not on PATH")
class TestIngestEndToEnd:
    def test_full_ingest_mini_species(self, tmp_path: Path, genome: Path):
        from mkprobes.ext.dataset import Dataset, load_dataset

        ds_dir = tmp_path / "mini_ds"
        runner = CliRunner()
        res = runner.invoke(
            ingest,
            [
                str(ds_dir),
                "--genome",
                str(genome),
                "--gtf",
                str(DATA / "ncbi_style.gtf"),
                "--species",
                "mini",
                "--blocklist-biotypes",
                "rRNA",
            ],
            catch_exceptions=False,
        )
        assert res.exit_code == 0, res.output

        for f in (
            "annotation.gtf",
            "transcripts.fasta",
            "transcripts.parquet",
            "transcripts.jf",
            "blocklist15.jf",
            "dataset.json",
            "solar_intake.yaml",
            "validation_report.json",
        ):
            assert (ds_dir / f).exists(), f"missing {f}"

        d = json.loads((ds_dir / "dataset.json").read_text())
        assert d["species"] == "mini"
        assert d["external_data"]["default"]["gtf_name"] == "annotation.gtf"
        assert d["external_data"]["default"]["strip_version"] is False
        assert d["blocklist_kmer_name"] == "blocklist15.jf"

        ds = load_dataset(ds_dir)
        assert type(ds) is Dataset
        # 2700 nt spliced g1.t1 (ncbi id), exact length check
        seq = ds.data.get_seq("rna-XM_001001.1", convert=False)
        assert len(seq) == 800 + 900 + 1000
        # rRNA 15-mers are blocklisted
        rrna_seq = ds.data.get_seq(
            ds.data.gtf.filter(ds.data.gtf["gene_biotype"] == "rRNA")[0, "transcript_id"],
            convert=False,
        )
        assert ds.check_kmers(rrna_seq)


class TestUnstrandedTranscripts:
    """
    StringTie emits `.` for single-exon transcripts it cannot orient. gffread
    keeps those but extracts the plus strand for them, so probes for any that
    are really minus-strand are antisense and never bind - and nothing
    downstream notices, because the sequence itself is valid.
    """

    def _gtf(self, tmp_path: Path, strands: list[str]) -> Path:
        rows = []
        for i, strand in enumerate(strands):
            attrs = f'gene_id "g{i}"; transcript_id "t{i}";'
            start, end = 100 + i * 200, 200 + i * 200
            for feature in ("transcript", "exon"):
                rows.append(f"chr1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
        path = tmp_path / "annotation.gtf"
        path.write_text("\n".join(rows) + "\n")
        return path

    def _genome(self, tmp_path: Path) -> Path:
        path = tmp_path / "genome.fa"
        path.write_text(">chr1\n" + "ACGT" * 500 + "\n")
        return path

    def test_counts_transcripts_not_rows(self, tmp_path: Path):
        from mkprobes.ext.ingest import validate_gtf

        # Two unstranded transcripts, each with an exon: four rows, two transcripts.
        report = validate_gtf(self._gtf(tmp_path, [".", ".", "+"]), self._genome(tmp_path))

        assert report.n_unstranded_transcripts == 2
        strand_issue = next(i for i in report.issues if i.code == "STRAND_MISSING")
        assert "2 transcript(s)" in strand_issue.message
        assert "4 rows" not in strand_issue.message

    def test_reports_the_share_of_the_annotation(self, tmp_path: Path):
        from mkprobes.ext.ingest import validate_gtf

        report = validate_gtf(self._gtf(tmp_path, [".", "+", "+", "-"]), self._genome(tmp_path))

        strand_issue = next(i for i in report.issues if i.code == "STRAND_MISSING")
        assert "25.0%" in strand_issue.message

    def test_says_what_actually_happens(self, tmp_path: Path):
        from mkprobes.ext.ingest import validate_gtf

        report = validate_gtf(self._gtf(tmp_path, ["."]), self._genome(tmp_path))

        message = next(i for i in report.issues if i.code == "STRAND_MISSING").message
        # gffread keeps them and silently uses the plus strand; it does not skip.
        assert "kept" in message
        assert "skips" not in message
        assert "antisense" in message

    def test_collects_the_ids(self, tmp_path: Path):
        from mkprobes.ext.ingest import validate_gtf

        report = validate_gtf(self._gtf(tmp_path, ["+", ".", "-", "."]), self._genome(tmp_path))

        assert report.unstranded_transcript_ids == ["t1", "t3"]

    def test_ids_stay_out_of_the_json_report(self, tmp_path: Path):
        import json

        from mkprobes.ext.ingest import validate_gtf

        report = validate_gtf(self._gtf(tmp_path, ["."] * 3), self._genome(tmp_path))

        # There can be tens of thousands; they belong in their own file.
        serialized = json.loads(report.model_dump_json())
        assert "unstranded_transcript_ids" not in serialized
        assert serialized["n_unstranded_transcripts"] == 3

    def test_silent_when_every_transcript_is_stranded(self, tmp_path: Path):
        from mkprobes.ext.ingest import validate_gtf

        report = validate_gtf(self._gtf(tmp_path, ["+", "-", "+"]), self._genome(tmp_path))

        assert report.n_unstranded_transcripts == 0
        assert not [i for i in report.issues if i.code == "STRAND_MISSING"]


class TestOverwriteGuard:
    """
    `--validate-only` writes annotation.gtf into the dataset directory. Guarding
    the real ingest on that file made the documented flow - validate, then
    ingest - fail on its own leftovers, so the guard is on dataset.json, the
    marker that a build actually completed.
    """

    def _inputs(self, tmp_path: Path) -> tuple[Path, Path]:
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\n" + "ACGT" * 500 + "\n")
        gtf = tmp_path / "in.gtf"
        gtf.write_text(
            'chr1\tt\ttranscript\t100\t200\t.\t+\t.\tgene_id "g0"; transcript_id "t0";\n'
            'chr1\tt\texon\t100\t200\t.\t+\t.\tgene_id "g0"; transcript_id "t0";\n'
        )
        return genome, gtf

    def _run(self, runner: CliRunner, dataset_dir: Path, genome: Path, gtf: Path, *extra: str):
        from mkprobes import cli

        return runner.invoke(
            cli.main,
            ["ingest", str(dataset_dir), "--genome", str(genome), "--gtf", str(gtf),
             "--species", "test", *extra],
        )

    def test_validate_only_leftovers_do_not_block_a_real_ingest(self, tmp_path: Path):
        runner = CliRunner()
        genome, gtf = self._inputs(tmp_path)
        dataset_dir = tmp_path / "ds"

        first = self._run(runner, dataset_dir, genome, gtf, "--validate-only")
        assert first.exit_code == 0, first.output
        assert (dataset_dir / "annotation.gtf").exists()  # the leftover in question

        second = self._run(runner, dataset_dir, genome, gtf)

        # It may still fail later for want of external tools; what it must not
        # do is refuse up front because validation ran.
        assert "already holds a built dataset" not in flatten_cli_output(second.output)

    def test_a_built_dataset_is_still_protected(self, tmp_path: Path):
        runner = CliRunner()
        genome, gtf = self._inputs(tmp_path)
        dataset_dir = tmp_path / "ds"
        dataset_dir.mkdir()
        (dataset_dir / "dataset.json").write_text("{}")

        result = self._run(runner, dataset_dir, genome, gtf)

        assert result.exit_code != 0
        assert "already holds a built dataset" in flatten_cli_output(result.output)
        assert "--overwrite" in flatten_cli_output(result.output)

    def test_overwrite_bypasses_the_guard(self, tmp_path: Path):
        runner = CliRunner()
        genome, gtf = self._inputs(tmp_path)
        dataset_dir = tmp_path / "ds"
        dataset_dir.mkdir()
        (dataset_dir / "dataset.json").write_text("{}")

        result = self._run(runner, dataset_dir, genome, gtf, "--overwrite", "--validate-only")

        assert "already holds a built dataset" not in flatten_cli_output(result.output)
