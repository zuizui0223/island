"""Classify Chapter 1 when/where support without conflating no test with a null result."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def classify_when_where_support(
    within: pd.DataFrame,
    contexts: list[str] | tuple[str, ...],
) -> pd.DataFrame:
    required = {"stratum", "support_tier", "context", "where_supported"}
    missing = required - set(within.columns)
    if missing:
        raise typer.BadParameter(f"within table missing columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for context in contexts:
        entry: dict[str, object] = {"context": context}
        for tier in ("confirmatory", "pilot"):
            for stratum in ("all_native", "native_nonendemic", "endemic"):
                subset = within.loc[
                    within["support_tier"].eq(tier)
                    & within["stratum"].eq(stratum)
                    & within["context"].eq(context)
                ]
                testable = not subset.empty
                supported = bool(testable and subset.iloc[0]["where_supported"])
                entry[f"{tier}_testable_{stratum}"] = testable
                entry[f"{tier}_supported_{stratum}"] = supported

        conf_testable = any(
            bool(entry[f"confirmatory_testable_{s}"])
            for s in ("all_native", "native_nonendemic", "endemic")
        )
        conf_supported = any(
            bool(entry[f"confirmatory_supported_{s}"])
            for s in ("all_native", "native_nonendemic", "endemic")
        )
        pilot_testable = any(
            bool(entry[f"pilot_testable_{s}"])
            for s in ("all_native", "native_nonendemic", "endemic")
        )
        pilot_supported = any(
            bool(entry[f"pilot_supported_{s}"])
            for s in ("all_native", "native_nonendemic", "endemic")
        )

        if (
            bool(entry["confirmatory_supported_all_native"])
            and bool(entry["confirmatory_supported_native_nonendemic"])
        ):
            cls = "confirmatory_persists_in_native_nonendemic"
        elif (
            bool(entry["confirmatory_supported_endemic"])
            and not bool(entry["confirmatory_supported_native_nonendemic"])
        ):
            cls = "confirmatory_endemic_concentrated"
        elif conf_supported:
            cls = "confirmatory_status_limited"
        elif conf_testable:
            cls = "confirmatory_tested_null"
        elif pilot_supported:
            cls = "pilot_signal_confirmatory_not_testable"
        elif pilot_testable:
            cls = "pilot_tested_null_confirmatory_not_testable"
        else:
            cls = "not_testable_current_support"

        entry["when_class"] = cls
        rows.append(entry)
    return pd.DataFrame(rows)


@app.command("run")
def run(
    within_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    contexts: str = typer.Option(
        "northern_midlatitude,northern_high_latitude,tropical,southern_extratropical"
    ),
) -> None:
    result = classify_when_where_support(
        pd.read_csv(within_csv),
        [x.strip() for x in contexts.split(",") if x.strip()],
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)


if __name__ == "__main__":
    app()
