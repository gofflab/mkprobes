"""
Design parameters: the thermodynamic and tiling settings a panel is designed under.

The crawler, the split into two arms, and the screen each carry thresholds
that were calibrated for mammalian transcripts at about 50% GC. On an AT-rich
transcriptome (cephalopods sit near 36% GC) those thresholds discard most of a
transcript before any off-target check runs, and none of them used to be
reachable without editing code.

They are reachable from two places, resolved in this order:

1. flags on `mkprobes run-panel` or `mkprobes candidates`, for one run;
2. the `"design"` block of the probe set in `manifest.json`, which is where a
   panel's values belong so that design and assembly agree on them;
3. the built-in defaults, which are unchanged from before these were exposed.

A field left out (or `null`) means "use the default", so an existing manifest
or command line designs exactly what it did before.
"""

from __future__ import annotations

from typing import Any

import click
from pydantic import BaseModel, ConfigDict, field_validator

#: primer3 refuses hairpin calculations on anything longer than this, so the
#: crawler cannot emit a longer probe whatever the length range says.
MAX_PROBE_LENGTH = 60

#: Crawler Tm window in degrees C at the design formamide concentration. A window
#: grows until it clears the low end, and is dropped if it exceeds the high end.
DEFAULT_TM_RANGE: tuple[float, float] = (54.0, 68.0)

#: Probe length window in nt. The reference (human/mouse) path and the custom
#: path inherited different upper bounds from the original drivers; both are
#: kept so that neither designs differently than it did.
DEFAULT_LENGTH_RANGE_REFERENCE: tuple[int, int] = (43, 55)
DEFAULT_LENGTH_RANGE_CUSTOM: tuple[int, int] = (43, 54)

#: Tm each arm of the split probe has to reach (within 27 nt) before the
#: candidate is kept, in degrees C without formamide.
DEFAULT_SPLIT_TM: float = 60.0

#: Screen-stage target: probes per gene the overlap search tries to reach.
DEFAULT_MIN_PROBES = 60

#: How far the overlap search may go (nt of overlap between neighbouring
#: probes, in steps of 5) to reach `min_probes`. 0 disables the search.
DEFAULT_MAX_OVERLAP = 0


class DesignParameters(BaseModel):
    """
    Settings that shape which probes a target yields.

    Every field is optional; `None` means the built-in default. `resolve` turns
    a partial set into concrete values, and `merged` layers one set over
    another, which is how a command-line flag overrides the manifest.
    """

    model_config = ConfigDict(extra="forbid")

    #: Crawler Tm window (low, high), degrees C at the design formamide
    #: concentration. Lowering the low end admits weaker-binding probes.
    tm_range: tuple[float, float] | None = None
    #: Probe length window (min, max) in nt. Raising the maximum lets an
    #: AT-rich window grow long enough to reach the Tm floor; 60 is the ceiling.
    length_range: tuple[int, int] | None = None
    #: Tm each arm of the split probe must reach, degrees C without formamide.
    split_tm: float | None = None
    #: Screen-stage target number of probes per gene.
    min_probes: int | None = None
    #: Maximum overlap between neighbouring probes the screen may use to reach
    #: `min_probes`; a multiple of 5, 0 to disable.
    max_overlap: int | None = None

    @field_validator("tm_range")
    @classmethod
    def _tm_range_ordered(cls, value: tuple[float, float] | None):
        if value is not None and not value[0] < value[1]:
            raise ValueError(f"tm_range must be (low, high) with low < high, got {value}.")
        return value

    @field_validator("length_range")
    @classmethod
    def _length_range_sane(cls, value: tuple[int, int] | None):
        if value is None:
            return value
        low, high = value
        if not 1 <= low <= high:
            raise ValueError(f"length_range must be (min, max) with 1 <= min <= max, got {value}.")
        if high > MAX_PROBE_LENGTH:
            raise ValueError(
                f"length_range maximum is {high}, but primer3 cannot evaluate hairpins above "
                f"{MAX_PROBE_LENGTH} nt, so probes longer than that can never be designed."
            )
        return value

    @field_validator("split_tm")
    @classmethod
    def _split_tm_positive(cls, value: float | None):
        if value is not None and value <= 0:
            raise ValueError(f"split_tm must be positive, got {value}.")
        return value

    @field_validator("min_probes")
    @classmethod
    def _min_probes_positive(cls, value: int | None):
        if value is not None and value < 1:
            raise ValueError(f"min_probes must be at least 1, got {value}.")
        return value

    @field_validator("max_overlap")
    @classmethod
    def _max_overlap_step(cls, value: int | None):
        if value is not None and (value < 0 or value % 5):
            raise ValueError(f"max_overlap must be 0 or a positive multiple of 5, got {value}.")
        return value

    def merged(self, override: DesignParameters | None) -> DesignParameters:
        """These parameters with every field set in `override` replacing its own."""
        if override is None:
            return self
        return DesignParameters(**{**self.explicit(), **override.explicit()})

    def explicit(self) -> dict[str, Any]:
        """Only the fields that were actually set, as plain JSON-able values."""
        return self.model_dump(exclude_none=True)

    def is_default(self) -> bool:
        return not self.explicit()

    def resolve(self, *, reference: bool) -> ResolvedDesign:
        """
        Concrete values for one dataset kind.

        `reference` selects the length window: the human/mouse path and the
        custom path have different upper bounds by inheritance.
        """
        return ResolvedDesign(
            tm_range=self.tm_range or DEFAULT_TM_RANGE,
            length_range=self.length_range
            or (DEFAULT_LENGTH_RANGE_REFERENCE if reference else DEFAULT_LENGTH_RANGE_CUSTOM),
            split_tm=self.split_tm if self.split_tm is not None else DEFAULT_SPLIT_TM,
            min_probes=self.min_probes if self.min_probes is not None else DEFAULT_MIN_PROBES,
            max_overlap=self.max_overlap if self.max_overlap is not None else DEFAULT_MAX_OVERLAP,
        )

    def describe(self) -> str:
        """One line for a human, naming only what was set."""
        if self.is_default():
            return "defaults"
        parts = []
        for key, value in self.explicit().items():
            shown = f"{value[0]}-{value[1]}" if isinstance(value, (tuple, list)) else value
            parts.append(f"{key} {shown}")
        return ", ".join(parts)


class ResolvedDesign(BaseModel):
    """`DesignParameters` with every field filled in."""

    model_config = ConfigDict(frozen=True)

    tm_range: tuple[float, float]
    length_range: tuple[int, int]
    split_tm: float
    min_probes: int
    max_overlap: int


def matches_recorded(design: DesignParameters, recorded: dict[str, Any] | None) -> bool:
    """
    Whether an output recorded as designed under `recorded` came from `design`.

    Outputs record the resolved values. The two dataset kinds resolve the
    length window differently and the reader here (the panel driver, assembly)
    does not always know which kind produced the file, so either resolution
    counts. Files with nothing recorded predate this and are trusted.
    """
    if recorded is None:
        return True
    return any(recorded == design.resolve(reference=kind).model_dump(mode="json") for kind in (True, False))


class RangeParam(click.ParamType):
    """
    A `LOW,HIGH` pair on the command line, e.g. `--tm-range 50,68`.

    A hyphen also separates the two, so `43-60` reads the way people write a
    range, but a comma is what the help text shows because `50-68` next to a
    negative number would be ambiguous.
    """

    name = "LOW,HIGH"

    def __init__(self, kind: type[int | float]):
        self.kind = kind

    def convert(self, value: Any, param: click.Parameter | None, ctx: click.Context | None):
        if isinstance(value, tuple):
            return value
        text = str(value).strip()
        parts = text.split(",") if "," in text else text.split("-")
        if len(parts) != 2:
            self.fail(f"expected two numbers as LOW,HIGH, got {value!r}.", param, ctx)
        try:
            low, high = (self.kind(p.strip()) for p in parts)
        except ValueError:
            self.fail(f"expected two numbers as LOW,HIGH, got {value!r}.", param, ctx)
        return (low, high)


TM_RANGE = RangeParam(float)
LENGTH_RANGE = RangeParam(int)


def design_options(command):
    """
    The thermodynamic flags shared by `candidates` and `run-panel`.

    Applied as a decorator. Every flag defaults to `None` so a command can tell
    "not given" from "given the default value", which is what lets the
    manifest supply a value that a flag can still override.
    """
    options = [
        click.option(
            "--tm-range",
            type=TM_RANGE,
            default=None,
            metavar="LOW,HIGH",
            help=f"Crawler Tm window in °C at the design formamide concentration "
            f"[default: {DEFAULT_TM_RANGE[0]:g},{DEFAULT_TM_RANGE[1]:g}]. Lower LOW on AT-rich "
            "transcripts to admit probes that bind less tightly.",
        ),
        click.option(
            "--length-range",
            type=LENGTH_RANGE,
            default=None,
            metavar="MIN,MAX",
            help=f"Probe length window in nt [default: {DEFAULT_LENGTH_RANGE_REFERENCE[0]},"
            f"{DEFAULT_LENGTH_RANGE_REFERENCE[1]} for human/mouse, "
            f"{DEFAULT_LENGTH_RANGE_CUSTOM[0]},{DEFAULT_LENGTH_RANGE_CUSTOM[1]} otherwise]. "
            f"MAX cannot exceed {MAX_PROBE_LENGTH}. Raise it so AT-rich windows can reach the Tm floor.",
        ),
        click.option(
            "--split-tm",
            type=float,
            default=None,
            metavar="TM",
            help=f"Tm in °C each arm of the split probe must reach [default: {DEFAULT_SPLIT_TM:g}]. "
            "The single largest lever on AT-rich transcripts; lower it with care, both arms "
            "must bind for ligation.",
        ),
    ]
    for option in reversed(options):
        command = option(command)
    return command


def design_from_flags(
    tm_range: tuple[float, float] | None,
    length_range: tuple[int, int] | None,
    split_tm: float | None,
    min_probes: int | None = None,
    max_overlap: int | None = None,
) -> DesignParameters:
    """The parameters a command line asked for, reported as a validation error if unusable."""
    try:
        return DesignParameters(
            tm_range=tm_range,
            length_range=length_range,
            split_tm=split_tm,
            min_probes=min_probes,
            max_overlap=max_overlap,
        )
    except ValueError as e:
        raise click.UsageError(_first_message(e)) from e


def _first_message(error: Exception) -> str:
    """Pydantic's report, cut down to the sentence a person needs."""
    try:
        errors = error.errors()  # type: ignore[attr-defined]
    except AttributeError:
        return str(error)
    if not errors:
        return str(error)
    message = str(errors[0].get("msg", error))
    return message.removeprefix("Value error, ")
