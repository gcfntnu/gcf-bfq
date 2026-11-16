"""BFQ's known-instrument requirement, using configmaker's shared registry."""


class SequencerError(ValueError):
    """A run cannot be associated with a registered sequencing instrument."""


def require_sequencer(run_id):
    """Resolve the instrument or fail before BFQ starts processing."""
    from configmaker.configmaker import SEQUENCERS  # noqa: PLC0415

    parts = run_id.split("_") if isinstance(run_id, str) else []
    if len(parts) < 2 or not parts[0] or not parts[1].strip():
        raise SequencerError(
            f"Cannot determine sequencer for run {run_id!r}: malformed run ID; "
            "expected an instrument ID in the second underscore-separated field. "
            "Check the run identity before retrying."
        )
    instrument = parts[1]
    machine = SEQUENCERS.get(instrument, "")
    if not machine:
        raise SequencerError(
            f"Cannot determine sequencer for run {run_id!r}: instrument ID {instrument!r} "
            "has no machine mapping in configmaker's SEQUENCERS (gcf-tools). "
            "Check the run identity and install gcf-tools with the correct instrument "
            "mapping before retrying. BFQ requires a known sequencer."
        )
    return machine
