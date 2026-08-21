"""Regression tests for the two probe-tiling selectors.

`find_overlap` (greedy, priority tier 1) and `OverlapWeighted.q` (weighted,
tiers 2+) currently disagree on the sign of `overlap`: the weighted selector
treats probe `i` as compatible with `j` when `end[i] < start[j] + overlap`,
while the greedy selector uses `end[i] < start[j] - overlap`. They agree only
at `overlap == 0`.

`overlap` means *permitted overlap in nt* -- higher is more permissive. That is
what `--maxoverlap` is documented as, what `OverlapWeighted.q`'s own comment
says, and what `run_screen(minimum=...)` relies on when it escalates `overlap`
through `(-2, 5, 10, 15, 20)` to reach `--minimum`. The greedy selector responds
to it backwards.

The sign fix is approved pending confirmation against a real crawl, so the tests
pinning the corrected convention are `xfail(strict=True)`. They turn into XPASS
the moment `find_overlap` is fixed, which is the signal to drop the marks and to
delete `test_greedy_currently_permits_overlap_at_negative_overlap` below.

See docs/notes/probe_selection_defects.md.
"""

from itertools import pairwise

import pytest

from mkprobes.utils._algorithms import find_overlap, find_overlap_weighted

PENDING_SIGN_FIX = pytest.mark.xfail(
    strict=True,
    reason="find_overlap's sign convention is inverted; fix pending real-crawl confirmation",
)


def tiling(n: int, step: int, length: int = 50) -> tuple[list[int], list[int]]:
    """`n` probes of `length` nt starting every `step` nt. Ends are sorted."""
    start = list(range(0, n * step, step))
    return start, [s + length - 1 for s in start]


def assert_respects_overlap(sel: list[int], start: list[int], end: list[int], overlap: int) -> None:
    """Every adjacent pair in a selection must satisfy `end[prev] < start[cur] + overlap`.

    This is the assay constraint: at `overlap=0` two selected probes share no
    base, at `overlap=-2` they are separated by a 2 nt gap, and at `overlap=5`
    they may share up to 5 nt.
    """
    assert sel == sorted(sel), "selection must be returned in position order"
    for prev, cur in pairwise(sel):
        assert end[prev] < start[cur] + overlap, (
            f"probes {prev} ({start[prev]}-{end[prev]}) and {cur} ({start[cur]}-{end[cur]}) "
            f"violate overlap={overlap}"
        )


class TestWeightedSelector:
    """`OverlapWeighted` -- the convention the assay requires."""

    @pytest.mark.parametrize("overlap", [-10, -5, -2, 0, 5, 10])
    def test_selection_respects_overlap(self, overlap: int):
        start, end = tiling(80, step=5)
        sel = find_overlap_weighted(start, end, [1.0] * len(start), overlap=overlap)
        assert_respects_overlap(sel, start, end, overlap)

    def test_overlap_zero_means_no_shared_base(self):
        # B starts exactly where A ends (shares 1 nt); C starts 1 nt after.
        start, end = [0, 49, 50], [49, 98, 99]
        assert find_overlap_weighted(start, end, [1.0] * 3, overlap=0) == [0, 2]

    def test_negative_overlap_forces_a_gap(self):
        # B overlaps A by 2 nt; C starts 3 nt after A ends.
        start, end = [0, 48, 52], [49, 97, 101]
        assert find_overlap_weighted(start, end, [1.0] * 3, overlap=-2) == [0, 2]

    def test_positive_overlap_permits_sharing(self):
        start, end = [0, 48, 52], [49, 97, 101]
        assert find_overlap_weighted(start, end, [1.0] * 3, overlap=5) == [0, 1]

    def test_probe_count_is_monotonic_in_overlap(self):
        """Raising `overlap` must never reduce the probe count -- `run_screen`
        escalates it precisely to gain probes."""
        start, end = tiling(80, step=5)
        counts = [
            len(find_overlap_weighted(start, end, [1.0] * len(start), overlap=o))
            for o in (-10, -5, -2, 0, 5, 10, 20)
        ]
        assert counts == sorted(counts), counts
        assert counts[-1] > counts[0], "escalating overlap should eventually gain probes"

    def test_prefers_higher_priority_probe(self):
        """With two mutually exclusive candidates the higher weight wins."""
        start, end = [0, 10], [49, 59]
        assert find_overlap_weighted(start, end, [1.0, 5.0], overlap=0) == [1]
        assert find_overlap_weighted(start, end, [5.0, 1.0], overlap=0) == [0]

    def test_rejects_unsorted_ends(self):
        with pytest.raises(ValueError, match="Ends not sorted"):
            find_overlap_weighted([0, 10], [99, 49], [1.0, 1.0])

    def test_rejects_mismatched_lengths(self):
        with pytest.raises(ValueError, match="Lengths not equal"):
            find_overlap_weighted([0, 10], [49, 59], [1.0])


class TestGreedySelector:
    def test_agrees_with_weighted_at_overlap_zero(self):
        """The one value at which the two conventions coincide."""
        start, end = tiling(60, step=7)
        assert find_overlap(start, end, overlap=0) == find_overlap_weighted(
            start, end, [1.0] * len(start), overlap=0
        )

    def test_selection_respects_overlap_zero(self):
        start, end = tiling(60, step=7)
        assert_respects_overlap(find_overlap(start, end, overlap=0), start, end, 0)

    def test_rejects_unsorted_ends(self):
        with pytest.raises(ValueError, match="Ends not sorted"):
            find_overlap([0, 10], [99, 49])

    def test_greedy_currently_permits_overlap_at_negative_overlap(self):
        """Pins the defect so it cannot change unnoticed. Delete when fixed.

        With `overlap=-2` the greedy selector picks A+B, which share 2 nt of
        target sequence, where the weighted selector picks the gapped A+C.
        """
        start, end = [0, 48, 52], [49, 97, 101]
        assert find_overlap(start, end, overlap=-2) == [0, 1]
        assert find_overlap_weighted(start, end, [1.0] * 3, overlap=-2) == [0, 2]

    @pytest.mark.parametrize(
        "overlap",
        [
            # Negative overlap: the inverted sign makes the greedy selector too
            # loose, so it emits probes that share target sequence.
            pytest.param(-10, marks=PENDING_SIGN_FIX),
            pytest.param(-5, marks=PENDING_SIGN_FIX),
            pytest.param(-2, marks=PENDING_SIGN_FIX),
            # Positive overlap: the inverted sign makes it too *strict*, so the
            # selection is still feasible -- it just leaves probes on the table.
            # That loss is caught by the monotonicity test below, not here.
            5,
            10,
        ],
    )
    def test_selection_respects_overlap(self, overlap: int):
        start, end = tiling(80, step=5)
        assert_respects_overlap(find_overlap(start, end, overlap=overlap), start, end, overlap)

    @PENDING_SIGN_FIX
    def test_probe_count_is_monotonic_in_overlap(self):
        start, end = tiling(80, step=5)
        counts = [len(find_overlap(start, end, overlap=o)) for o in (-10, -5, -2, 0, 5, 10, 20)]
        assert counts == sorted(counts), counts

    @PENDING_SIGN_FIX
    @pytest.mark.parametrize("overlap", [-10, -5, -2, 5, 10])
    def test_matches_weighted_feasibility_for_nonzero_overlap(self, overlap: int):
        """Both selectors must accept the same set of adjacent-probe pairs."""
        start, end = [0, 48, 52], [49, 97, 101]
        greedy = find_overlap(start, end, overlap=overlap)
        weighted = find_overlap_weighted(start, end, [1.0] * 3, overlap=overlap)
        assert_respects_overlap(greedy, start, end, overlap)
        assert len(greedy) == len(weighted)
