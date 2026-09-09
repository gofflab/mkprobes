"""
primer3 must never be entered from two threads at once.

libprimer3 keeps its hairpin working buffers in globals and primer3-py releases
the GIL around the call, while polars spreads a `map_elements` UDF over its
thread pool once a frame has several chunks. A hairpin column over a real
candidate table put sixteen threads inside primer3 together and the worker
died with a bus error, taking a 500-target panel down with it. The call is now
serialised; this runs the crashing shape in a subprocess so a regression shows
as a failed test rather than a dead pytest.
"""

import subprocess
import sys

CRASHING_SHAPE = """
import random
import polars as pl
from mkprobes.utils.seqcalc import hp

random.seed(0)
chunks = [
    pl.DataFrame({"seq": ["".join(random.choices("ACGT", k=random.randint(18, 60))) for _ in range(200)]})
    for _ in range(32)
]
df = pl.concat(chunks, rechunk=False)
assert df.n_chunks() > 1
out = df.with_columns(hp=pl.col("seq").map_elements(lambda s: hp(s, "hybrid", formamide=40), return_dtype=pl.Float64))
assert out["hp"].null_count() == 0
print("ok")
"""


def test_hairpin_over_a_chunked_frame_does_not_crash_the_process():
    result = subprocess.run(
        [sys.executable, "-c", CRASHING_SHAPE], capture_output=True, text=True, timeout=300, check=False
    )
    assert result.returncode == 0, f"exit {result.returncode} (138 is SIGBUS)\n{result.stderr[-2000:]}"
    assert result.stdout.strip().endswith("ok")
