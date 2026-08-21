DT = "TTACACTCCATCCACTCAA"
SP6 = "ATTTAGGTGACACTATAG"
GOOD_SPECIES = {"human", "mouse"}

#: The enzyme pair the SOLAR construct is built around.
#:
#: This is assay chemistry, not a tunable parameter. The header/footer table
#: carries these two sites so the probe can be excised, and final assembly
#: releases it with a KpnI/BamHI double digest and asserts on the resulting
#: geometry. Designing against a different pair would produce probes that no
#: downstream step can cut out, so the CLI accepts only this pair - see
#: `validate_restriction`.
SOLAR_RESTRICTION: tuple[str, str] = ("BamHI", "KpnI")

#: How the enzyme pair appears in output filenames.
RESTRICTION_TOKEN = "".join(SOLAR_RESTRICTION)


def validate_restriction(enzymes: "tuple[str, ...] | list[str] | None") -> tuple[str, ...]:
    """
    Checks a requested enzyme set against the chemistry, returning it unchanged.

    Raises `ValueError` naming the mismatch. Order is not significant, but the
    set is: anything other than the SOLAR pair fails here rather than after the
    whole panel has been computed.
    """
    if not enzymes:
        return ()
    requested = tuple(enzymes)
    if set(requested) != set(SOLAR_RESTRICTION):
        raise ValueError(
            f"SOLAR probes are built for {' + '.join(SOLAR_RESTRICTION)}, but "
            f"{' + '.join(requested)} was requested. The header/footer sequences carry "
            "these two sites and final assembly excises the probe with a KpnI/BamHI "
            "double digest, so a different pair yields probes that cannot be cut out. "
            "Drop the option to use the default."
        )
    return requested
