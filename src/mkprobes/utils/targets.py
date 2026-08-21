"""
Reading target lists.

A target list is the one file every user writes by hand, so it has to tolerate
the things people put in hand-written files: blank lines, trailing newlines, and
comments explaining why a gene is on the list. Reading it with `splitlines()`
turned a trailing newline into an empty target and a `#` note into a gene name
nobody could look up.
"""

from pathlib import Path


def read_target_list(path: Path | str) -> list[str]:
    """
    Reads a target list: one gene or transcript per line.

    Blank lines are skipped, `#` starts a comment, and inline comments are
    stripped. Order is preserved, because it decides bit assignment downstream.

    Raises `ValueError` when the file holds no targets, or names any target
    twice - a duplicate would take a second set of bits and quietly break the
    codebook.
    """
    path = Path(path)
    targets: list[str] = []
    for line in path.read_text().splitlines():
        entry = line.split("#", 1)[0].strip()
        if entry:
            targets.append(entry)

    if not targets:
        raise ValueError(f"{path} lists no targets. Add one gene or transcript per line.")

    seen: set[str] = set()
    duplicated = sorted({t for t in targets if t in seen or seen.add(t)})
    if duplicated:
        raise ValueError(
            f"{path} lists {len(duplicated)} target(s) more than once: "
            f"{', '.join(duplicated)}. Remove the repeats before designing."
        )
    return targets
