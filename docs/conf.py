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

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"

# Keep docs buildable in lean environments and CI where scientific deps are absent.
autodoc_mock_imports = [
    "anndata",
    "Bio",
    "matplotlib",
    "mygene",
    "numpy",
    "pandas",
    "polars",
    "primer3",
    "pyarrow",
    "pyfastx",
    "pydantic",
    "questionary",
    "requests",
    "rich",
    "rich_click",
    "click",
    "loguru",
    "scipy",
    "seaborn",
    "sklearn",
    "tqdm",
    "yaml",
]

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
