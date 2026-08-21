"""Shared test helpers."""

import re

import pytest

_ANSI = re.compile(r"\x1b\[[0-9;]*m")
# Box-drawing block: rich_click frames errors in a panel, so a wrapped phrase
# comes back with border characters sitting between its words.
_BOX = re.compile(r"[─-╿]")


def flatten_cli_output(output: str) -> str:
    """
    CLI output with colour codes and panel borders stripped, and wrapping undone.

    rich_click renders errors into a bordered panel wrapped to the terminal
    width, which differs between a developer's terminal and CI. A phrase that
    fits on one line locally arrives in CI split across lines with `│` between
    its halves, so asserting on the raw output fails there for no real reason.
    Tests match against this instead.
    """
    return " ".join(_BOX.sub(" ", _ANSI.sub("", output)).split())


@pytest.fixture
def flat():
    """Returns :func:`flatten_cli_output`, for asserting on wrapped CLI text."""
    return flatten_cli_output
