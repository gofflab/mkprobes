"""
Design settings: where they come from, how they are checked, and that they
reach the stages that use them.

The crawler's Tm and length windows and the split-arm Tm decide how many
probes an AT-rich transcript yields, and none of them used to be reachable.
These pin the resolution order (flag over manifest over default), that an
unset value designs exactly what it did before, and that `construct` builds
from the overlap the screen settled on rather than always the no-overlap file.
"""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import click
import polars as pl
import pytest
from click.testing import CliRunner
from conftest import flatten_cli_output

from mkprobes import cli
from mkprobes.candidates import candidates
from mkprobes.codebook.codebook import ProbeSet
from mkprobes.codebook.finalconstruct import construct, pick_screened_overlap, screened_overlaps
from mkprobes.design import (
    DEFAULT_LENGTH_RANGE_CUSTOM,
    DEFAULT_LENGTH_RANGE_REFERENCE,
    DEFAULT_SPLIT_TM,
    DEFAULT_TM_RANGE,
    LENGTH_RANGE,
    TM_RANGE,
    DesignParameters,
    matches_recorded,
)
from mkprobes.init_project import check_manifest, manifest_stub
from mkprobes.run_panel import design_for_codebook, run_panel_cli
from mkprobes.utils.provenance import read_provenance


class TestDesignParameters:
    def test_unset_resolves_to_the_historical_defaults(self):
        reference = DesignParameters().resolve(reference=True)
        custom = DesignParameters().resolve(reference=False)

        assert reference.tm_range == custom.tm_range == DEFAULT_TM_RANGE
        assert reference.split_tm == custom.split_tm == DEFAULT_SPLIT_TM
        # The two paths inherited different length caps; both are preserved.
        assert reference.length_range == DEFAULT_LENGTH_RANGE_REFERENCE == (43, 55)
        assert custom.length_range == DEFAULT_LENGTH_RANGE_CUSTOM == (43, 54)
        assert custom.min_probes == 60
        assert custom.max_overlap == 0

    def test_set_fields_override_defaults(self):
        resolved = DesignParameters(split_tm=55, length_range=(43, 60)).resolve(reference=False)
        assert resolved.split_tm == 55
        assert resolved.length_range == (43, 60)
        assert resolved.tm_range == DEFAULT_TM_RANGE

    def test_merged_lets_the_override_win_field_by_field(self):
        manifest = DesignParameters(split_tm=55, tm_range=(50, 68))
        flags = DesignParameters(split_tm=58)

        merged = manifest.merged(flags)

        assert merged.split_tm == 58
        assert merged.tm_range == (50, 68)
        assert manifest.merged(None) is manifest

    def test_length_above_primer3_ceiling_is_refused(self):
        with pytest.raises(ValueError, match="60"):
            DesignParameters(length_range=(43, 61))

    def test_inverted_tm_range_is_refused(self):
        with pytest.raises(ValueError, match="low < high"):
            DesignParameters(tm_range=(68, 54))

    def test_overlap_must_step_by_five(self):
        with pytest.raises(ValueError, match="multiple of 5"):
            DesignParameters(max_overlap=7)
        assert DesignParameters(max_overlap=20).max_overlap == 20

    def test_unknown_field_is_refused(self):
        # A typo in the manifest must not silently design under the defaults.
        with pytest.raises(ValueError):
            DesignParameters(split_tmp=55)  # type: ignore[call-arg]

    def test_describe_names_only_what_was_set(self):
        assert DesignParameters().describe() == "defaults"
        assert (
            DesignParameters(tm_range=(50, 68), split_tm=55).describe() == "tm_range 50.0-68.0, split_tm 55.0"
        )

    def test_recorded_output_matches_either_dataset_kind(self):
        design = DesignParameters(split_tm=55)
        recorded_reference = design.resolve(reference=True).model_dump(mode="json")
        recorded_custom = design.resolve(reference=False).model_dump(mode="json")

        assert matches_recorded(design, recorded_reference)
        assert matches_recorded(design, recorded_custom)
        assert not matches_recorded(DesignParameters(), recorded_custom)
        # Outputs written before settings were recorded are trusted.
        assert matches_recorded(design, None)


class TestRangeFlag:
    def test_comma_and_hyphen_both_separate(self):
        assert TM_RANGE.convert("50,68", None, None) == (50.0, 68.0)
        assert LENGTH_RANGE.convert("43-60", None, None) == (43, 60)

    def test_anything_else_is_an_error(self):
        with pytest.raises(click.BadParameter, match="LOW,HIGH"):
            TM_RANGE.convert("50", None, None)
        with pytest.raises(click.BadParameter, match="LOW,HIGH"):
            LENGTH_RANGE.convert("43,sixty", None, None)


class TestManifestDesignBlock:
    def test_stub_writes_the_defaults_for_the_species(self):
        mouse = manifest_stub("p", "mouse")[0]["design"]
        squid = manifest_stub("p", "squid")[0]["design"]

        assert mouse["length_range"] == [43, 55]
        assert squid["length_range"] == [43, 54]
        assert squid["split_tm"] == DEFAULT_SPLIT_TM

    def test_stub_validates_and_designs_as_the_defaults(self, tmp_path: Path):
        path = tmp_path / "manifest.json"
        path.write_text(json.dumps(manifest_stub("p", "squid")))
        (tmp_path / "codebook.json").write_text("{}")

        [probeset] = check_manifest(path)

        assert probeset.design.resolve(reference=False) == DesignParameters().resolve(reference=False)

    def test_missing_block_means_defaults(self):
        probeset = ProbeSet(name="p", species="squid", codebook="codebook.json", bcidx=0)
        assert probeset.design.is_default()

    def test_bad_value_is_reported_by_check_manifest(self, tmp_path: Path):
        entry = manifest_stub("p", "squid")[0]
        entry["design"]["length_range"] = [43, 80]
        path = tmp_path / "manifest.json"
        path.write_text(json.dumps([entry]))
        (tmp_path / "codebook.json").write_text("{}")

        with pytest.raises(ValueError, match="not a valid manifest"):
            check_manifest(path)

    def test_check_manifest_cli_reports_the_design(self, tmp_path: Path):
        entry = manifest_stub("p", "squid")[0]
        entry["design"] = {"split_tm": 55}
        path = tmp_path / "manifest.json"
        path.write_text(json.dumps([entry]))
        (tmp_path / "codebook.json").write_text("{}")

        result = CliRunner().invoke(cli.main, ["check-manifest", str(path)])

        assert result.exit_code == 0, result.output
        assert "design split_tm 55" in flatten_cli_output(result.output)


def _write_manifest(directory: Path, codebook: str = "codebook.json", **design) -> Path:
    entry = {"name": "p", "species": "squid", "codebook": codebook, "bcidx": 0, "design": design}
    path = directory / "manifest.json"
    path.write_text(json.dumps([entry]))
    return path


class TestDesignForCodebook:
    def test_manifest_beside_the_codebook_is_used(self, tmp_path: Path):
        codebook = tmp_path / "codebook.json"
        codebook.write_text("{}")
        manifest = _write_manifest(tmp_path, split_tm=55)

        design, source = design_for_codebook(codebook, None)

        assert design.split_tm == 55
        assert source == manifest

    def test_no_manifest_means_defaults(self, tmp_path: Path):
        codebook = tmp_path / "codebook.json"
        codebook.write_text("{}")

        design, source = design_for_codebook(codebook, None)

        assert design.is_default()
        assert source is None

    def test_manifest_beside_the_codebook_that_does_not_list_it_is_ignored(self, tmp_path: Path):
        codebook = tmp_path / "other.json"
        codebook.write_text("{}")
        _write_manifest(tmp_path, codebook="codebook.json", split_tm=55)

        design, source = design_for_codebook(codebook, None)

        assert design.is_default()
        assert source is None

    def test_explicit_manifest_must_list_the_codebook(self, tmp_path: Path):
        codebook = tmp_path / "other.json"
        codebook.write_text("{}")
        manifest = _write_manifest(tmp_path, codebook="codebook.json", split_tm=55)

        with pytest.raises(ValueError, match="no probe set whose codebook"):
            design_for_codebook(codebook, manifest)

    def test_each_panel_in_a_manifest_keeps_its_own_settings(self, tmp_path: Path):
        # One manifest can describe several panels, designed differently: a
        # mouse panel under the defaults next to a squid panel with relaxed
        # thermodynamics. Each codebook resolves to its own block.
        for name in ("mouse.json", "squid.json"):
            (tmp_path / name).write_text("{}")
        entries = [
            {"name": "mouse", "species": "mouse", "codebook": "mouse.json", "bcidx": 0},
            {
                "name": "squid",
                "species": "squid",
                "codebook": "squid.json",
                "bcidx": 1,
                "design": {"split_tm": 55, "length_range": [43, 60], "max_overlap": 10},
            },
        ]
        (tmp_path / "manifest.json").write_text(json.dumps(entries))

        mouse, _ = design_for_codebook(tmp_path / "mouse.json", None)
        squid, _ = design_for_codebook(tmp_path / "squid.json", None)

        assert mouse.is_default()
        assert squid.split_tm == 55
        assert squid.max_overlap == 10
        assert [ps.design.describe() for ps in check_manifest(tmp_path / "manifest.json")] == [
            "defaults",
            "length_range 43-60, split_tm 55.0, max_overlap 10",
        ]

    def test_conflicting_entries_are_refused(self, tmp_path: Path):
        codebook = tmp_path / "codebook.json"
        codebook.write_text("{}")
        entries = [
            {
                "name": "a",
                "species": "squid",
                "codebook": "codebook.json",
                "bcidx": 0,
                "design": {"split_tm": 55},
            },
            {
                "name": "b",
                "species": "squid",
                "codebook": "codebook.json",
                "bcidx": 1,
                "design": {"split_tm": 58},
            },
        ]
        (tmp_path / "manifest.json").write_text(json.dumps(entries))

        with pytest.raises(ValueError, match="different design settings"):
            design_for_codebook(codebook, None)


class TestRunPanelCli:
    @pytest.fixture
    def project(self, tmp_path: Path) -> Path:
        (tmp_path / "codebook.json").write_text(json.dumps({"G": [1, 2, 3]}))
        return tmp_path

    def _invoke(self, project: Path, *args: str):
        with (
            patch(
                "mkprobes.run_panel.run_panel", return_value={"done": [], "skipped": [], "failed": []}
            ) as run,
            patch("mkprobes.ext.ingest.check_external_tools"),
        ):
            result = CliRunner().invoke(run_panel_cli, [str(project), str(project / "codebook.json"), *args])
        return result, run

    def test_flags_reach_the_driver(self, project: Path):
        result, run = self._invoke(
            project, "--tm-range", "50,68", "--length-range", "43,60", "--split-tm", "55"
        )

        assert result.exit_code == 0, result.output
        design = run.call_args.kwargs["design"]
        assert design.tm_range == (50, 68)
        assert design.length_range == (43, 60)
        assert design.split_tm == 55

    def test_nothing_given_designs_under_the_defaults(self, project: Path):
        result, run = self._invoke(project)

        assert result.exit_code == 0, result.output
        assert run.call_args.kwargs["design"].is_default()

    def test_manifest_supplies_settings_and_a_flag_overrides_one(self, project: Path):
        _write_manifest(project, split_tm=55, max_overlap=20, min_probes=40)

        result, run = self._invoke(project, "--split-tm", "58")

        assert result.exit_code == 0, result.output
        design = run.call_args.kwargs["design"]
        assert design.split_tm == 58
        assert design.max_overlap == 20
        assert design.min_probes == 40

    def test_bad_range_is_reported_as_usage(self, project: Path):
        result, _ = self._invoke(project, "--length-range", "43,80")

        assert result.exit_code == 2
        assert "60" in flatten_cli_output(result.output)

    def test_explicit_manifest_that_does_not_list_the_codebook_fails(self, project: Path, tmp_path: Path):
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        manifest = _write_manifest(elsewhere, codebook="another.json", split_tm=55)

        result, _ = self._invoke(project, "--manifest", str(manifest))

        assert result.exit_code != 0
        assert "no probe set whose codebook" in flatten_cli_output(result.output)


class TestCandidatesCli:
    def test_flags_reach_candidate_generation(self, tmp_path: Path):
        with (
            patch("mkprobes.candidates.get_candidates") as get_candidates,
            patch("mkprobes.candidates.load_dataset", return_value=object()),
            patch("mkprobes.ext.ingest.check_external_tools"),
        ):
            result = CliRunner().invoke(
                candidates, [str(tmp_path), "--gene", "G", "--split-tm", "55", "--length-range", "43-60"]
            )

        assert result.exit_code == 0, result.output
        design = get_candidates.call_args.kwargs["design"]
        assert design.split_tm == 55
        assert design.length_range == (43, 60)
        assert design.tm_range is None


def _screened(names: list[str]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "name": names,
            "padlock": ["ACGTACGTACGTACGTACGT"] * len(names),
            "pad_start": [20] * len(names),
            "seq": ["ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"] * len(names),
        }
    )


class TestConstructOverlap:
    """
    `screen --minimum` writes one file per overlap it tried. `construct` used
    to read the no-overlap file regardless, so the search could never reach the
    pool.
    """

    def test_finds_every_overlap_for_the_gene_only(self, tmp_path: Path):
        for name in (
            "G_screened_ol-2_BamHIKpnI.parquet",
            "G_screened_ol10_BamHIKpnI.parquet",
            "G_screened_ol-2.parquet",
            "G2_screened_ol20_BamHIKpnI.parquet",
        ):
            (tmp_path / name).touch()

        assert sorted(screened_overlaps(tmp_path, "G", "_BamHIKpnI")) == [-2, 10]
        assert sorted(screened_overlaps(tmp_path, "G", "")) == [-2]
        assert pick_screened_overlap(tmp_path, "G", "_BamHIKpnI") == 10

    def test_no_screened_file_is_named_clearly(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="mkprobes screen"):
            pick_screened_overlap(tmp_path, "G", "_BamHIKpnI")

    def test_builds_from_the_overlap_the_screen_settled_on(self, tmp_path: Path):
        _screened(["G_G:1-40"]).write_parquet(tmp_path / "G_screened_ol-2_BamHIKpnI.parquet")
        _screened(["G_G:1-40", "G_G:30-70"]).write_parquet(tmp_path / "G_screened_ol10_BamHIKpnI.parquet")
        dataset = SimpleNamespace(path=tmp_path)
        design = DesignParameters(max_overlap=10).resolve(reference=False)

        construct(
            dataset,
            tmp_path,
            transcript="G",
            codebook={"G": [1, 2, 3]},
            restriction=["BamHI", "KpnI"],
            design=design,
        )  # type: ignore[arg-type]

        final = tmp_path / "G_final_BamHIKpnI_1,2,3.parquet"
        assert len(pl.read_parquet(final)) == 2
        record = read_provenance(final)
        assert record is not None
        assert record["overlap"] == 10
        assert record["design"]["max_overlap"] == 10

    def test_an_explicit_overlap_wins(self, tmp_path: Path):
        _screened(["G_G:1-40"]).write_parquet(tmp_path / "G_screened_ol-2_BamHIKpnI.parquet")
        _screened(["G_G:1-40", "G_G:30-70"]).write_parquet(tmp_path / "G_screened_ol10_BamHIKpnI.parquet")

        construct(
            SimpleNamespace(path=tmp_path),
            tmp_path,
            transcript="G",
            codebook={"G": [1, 2, 3]},
            restriction=["BamHI", "KpnI"],
            overlap=-2,
        )  # type: ignore[arg-type]

        final = tmp_path / "G_final_BamHIKpnI_1,2,3.parquet"
        assert len(pl.read_parquet(final)) == 1
        assert read_provenance(final)["overlap"] == -2  # type: ignore[index]

    def test_a_missing_explicit_overlap_names_what_exists(self, tmp_path: Path):
        _screened(["G_G:1-40"]).write_parquet(tmp_path / "G_screened_ol-2_BamHIKpnI.parquet")

        with pytest.raises(FileNotFoundError, match=r"overlap 20.*\[-2\]"):
            construct(
                SimpleNamespace(path=tmp_path),
                tmp_path,
                transcript="G",
                codebook={"G": [1, 2, 3]},
                restriction=["BamHI", "KpnI"],
                overlap=20,
            )  # type: ignore[arg-type]


class TestNothingPassesTheFilter:
    def test_is_an_error_naming_the_gene_not_a_process_exit(self):
        from mkprobes.utils._filtration import the_filter

        # Every arm fails the priority floors, so no tier admits anything.
        arms = pl.DataFrame(
            {
                "name": ["G_G:1-40_splint", "G_G:1-40_padlock"],
                "gene": ["G", "G"],
                "pos_start": [1, 1],
                "pos_end": [40, 40],
                "tm": [60.0, 60.0],
                "oks": [0, 0],
                "hp": [10.0, 10.0],
                "max_tm_offtarget": [0.0, 0.0],
                "maps_to_pseudo": ["", ""],
            }
        )

        with pytest.raises(ValueError, match="No probes passed the filters for G"):
            the_filter(arms, overlap=-2)
