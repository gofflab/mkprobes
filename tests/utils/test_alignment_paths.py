"""
Subprocess invocations must survive paths containing spaces.

Commands were built by formatting a string and splitting it on whitespace, so
any dataset under a path like `My Drive/Goff Lab/` was passed to bowtie2 as
several arguments and every alignment failed with exit 255.
"""

import ast
import importlib
from pathlib import Path

import pytest

MODULES_INVOKING_TOOLS = [
    "mkprobes.utils._alignment",
    "mkprobes.assembly",
    "mkprobes.ext.prepare",
    "mkprobes.ext.ingest",
]


def _calls(module: str) -> list[ast.Call]:
    source = importlib.import_module(module).__loader__.get_source(module)
    return [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Call)]


@pytest.mark.parametrize("module", MODULES_INVOKING_TOOLS)
def test_no_command_is_built_by_splitting_a_string(module: str):
    """`shlex.split` on an interpolated path silently breaks it at spaces."""
    offenders = [
        node.lineno
        for node in _calls(module)
        if isinstance(node.func, ast.Attribute)
        and node.func.attr == "split"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "shlex"
        # A literal with no interpolation carries no user path.
        and node.args
        and not isinstance(node.args[0], ast.Constant)
    ]
    assert not offenders, f"{module}: command built by splitting a formatted string at {offenders}"


@pytest.mark.parametrize("module", MODULES_INVOKING_TOOLS)
def test_no_subprocess_uses_a_shell(module: str):
    """shell=True re-splits the command and needs every path escaped."""
    offenders = [
        node.lineno
        for node in _calls(module)
        for kw in node.keywords
        if kw.arg == "shell" and getattr(kw.value, "value", False) is True
    ]
    assert not offenders, f"{module}: shell=True at {offenders}"


def test_bowtie2_receives_the_index_path_as_one_argument(tmp_path: Path):
    """The regression itself: a spaced path must arrive intact."""
    import subprocess
    from unittest.mock import patch

    from mkprobes.utils._alignment import run_bowtie

    reference = "/tmp/My Drive/Goff Lab/transcripts"
    with patch("subprocess.run") as run:
        run.return_value = subprocess.CompletedProcess(args=[], returncode=0, stdout="")
        run_bowtie(">p\nACGT\n", reference, fasta=True)

    command = run.call_args[0][0]
    assert isinstance(command, list)
    assert reference in command, f"index path was mangled: {command}"
