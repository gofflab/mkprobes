import sys

from loguru import logger


def setup_logging(level: str = "INFO") -> None:
    """
    Configure the concise console logger used by the command-line interface.

    Defaults to INFO: DEBUG records carry tracebacks and internal detail that
    are noise on a bench scientist's terminal. `mkprobes --debug` raises it,
    and per-gene file sinks keep DEBUG regardless.
    """
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
        colorize=True,
        format=(
            "<green>{time:HH:mm:ss}</green> | "
            "<level>{level: >8}</level> | "
            "<cyan>{name}</cyan> | <level>{message}</level>"
        ),
    )

