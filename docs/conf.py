from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

project = "mkprobes"
author = "Goff Lab"
copyright = f"{datetime.now():%Y}, {author}"
release = "0.1.0"

# `sphinx_click.ext` renders reference/cli.md straight from the Click tree, so
# the CLI reference cannot drift from the code. It imports `mkprobes.cli` for
# real, which means the docs build must install the package and its
# dependencies - see docs/requirements.txt and .github/workflows/docs.yml.
#
# There is deliberately no autodoc/autosummary here: no page documents the
# Python API, and an autodoc stack configured with mocked-out imports fails in
# confusing ways the moment someone adds one. Add them back together with the
# page that needs them.
extensions = [
    "myst_parser",
    "sphinx_click.ext",
]

# Off because these pages are full of command-line flags: smartquotes rewrites
# a `--flag` mentioned in running text (including inside generated CLI help)
# into an en-dash, which is wrong and uncopyable.
smartquotes = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]
myst_heading_anchors = 4

html_theme = "sphinx_rtd_theme"
html_static_path = []
html_title = "mkprobes documentation"

os.environ.setdefault("PYTHONUNBUFFERED", "1")
