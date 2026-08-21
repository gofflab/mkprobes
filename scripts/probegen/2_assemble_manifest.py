# %%
"""
DEPRECATED shim: probeset assembly now lives in the package as
`mkprobes.assembly` (CLI: `mkprobes assemble <manifest> {short|gen}`).
"""

from mkprobes.assembly import (  # noqa: F401  (re-exports)
    backfill,
    cli,
    handle_checks,
    manual_accept,
    run,
)

if __name__ == "__main__":
    cli()
