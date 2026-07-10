from pathlib import Path
import re

path = Path("src/island_v2/global_trait_campaign.py")
text = path.read_text(encoding="utf-8")

text = text.replace(
    "import json\nimport re\nimport time\n",
    "import json\nimport os\nimport re\nimport time\nfrom datetime import datetime, timezone\nfrom email.utils import parsedate_to_datetime\n",
    1,
)

fetch_pattern = re.compile(
    r"def fetch_openalex_sources\(.*?(?=def _error_species)",
    re.DOTALL,
)
fetch_replacement = '''class OpenAlexTransientError(RuntimeError):
    """A retryable OpenAlex infrastructure, quota, or transport failure."""


def _retry_after_seconds(value: str, now: datetime | None = None) -> float | None:
    """Parse Retry-After seconds or an HTTP date into a non-negative delay."""
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return max(0.0, float(text))
    except ValueError:
        pass
    try:
        parsed = parsedate_to_datetime(text)
    except (TypeError, ValueError, OverflowError):
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    reference = now or datetime.now(timezone.utc)
    return max(0.0, (parsed - reference).total_seconds())


def _openalex_request(
    client: Any,
    params: dict[str, Any],
    *,
    max_retries: int,
    backoff_seconds: float,
    max_backoff_seconds: float,
    sleeper: Any = time.sleep,
) -> Any:
    """Issue one OpenAlex request with bounded transient retry handling."""
    import httpx

    retryable_statuses = {429, 500, 502, 503, 504}
    attempt = 0
    while True:
        try:
            response = client.get(ecology.OPENALEX_WORKS_URL, params=params)
        except httpx.TransportError as exc:
            if attempt >= max_retries:
                raise OpenAlexTransientError(f"transport:{exc}") from exc
            delay = min(max_backoff_seconds, backoff_seconds * (2**attempt))
            sleeper(delay)
            attempt += 1
            continue

        if response.status_code not in retryable_statuses:
            response.raise_for_status()
            return response

        retry_after = _retry_after_seconds(response.headers.get("Retry-After", ""))
        # A large Retry-After generally means the daily budget is exhausted. Do not
        # sleep through the workflow timeout or spend more search calls immediately.
        if response.status_code == 429 and (
            retry_after is None or retry_after > max_backoff_seconds
        ):
            raise OpenAlexTransientError(
                f"status=429 retry_after={retry_after if retry_after is not None else 'missing'}"
            )
        if attempt >= max_retries:
            raise OpenAlexTransientError(
                f"status={response.status_code} retries_exhausted={max_retries}"
            )
        delay = (
            retry_after
            if retry_after is not None
            else min(max_backoff_seconds, backoff_seconds * (2**attempt))
        )
        sleeper(delay)
        attempt += 1


def fetch_openalex_sources(
    batch: pd.DataFrame,
    task: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    import httpx

    terms = _openalex_terms(config, task)
    task_config = config["tasks"][task]
    per_page = int(task_config.get("per_page", 5))
    pause_seconds = float(task_config.get("pause_seconds", 0.25))
    max_retries = int(task_config.get("max_retries", 3))
    backoff_seconds = float(task_config.get("backoff_seconds", 5.0))
    max_backoff_seconds = float(task_config.get("max_backoff_seconds", 60.0))
    require_api_key = bool(task_config.get("require_api_key", True))
    api_key = os.environ.get("OPENALEX_API_KEY", "").strip()
    rows: list[dict[str, str]] = []
    errors: list[str] = []

    species_names = [str(value) for value in batch["accepted_species"].astype(str)]
    if require_api_key and not api_key:
        return (
            pd.DataFrame(rows, columns=ecology.TEXT_SOURCE_COLUMNS),
            [
                f"openalex:{species}:transient:missing_OPENALEX_API_KEY"
                for species in species_names
            ],
        )

    headers = {"User-Agent": "island-floral-v2 global-trait-campaign"}
    with httpx.Client(timeout=45.0, follow_redirects=True, headers=headers) as client:
        for species in species_names:
            params: dict[str, Any] = {
                "search": f'"{species}" {terms}',
                "per_page": per_page,
                "select": "id,doi,display_name,publication_year,abstract_inverted_index",
            }
            if api_key:
                params["api_key"] = api_key
            try:
                response = _openalex_request(
                    client,
                    params,
                    max_retries=max_retries,
                    backoff_seconds=backoff_seconds,
                    max_backoff_seconds=max_backoff_seconds,
                )
                works = response.json().get("results") or []
                for work in works:
                    title = _text(re.sub(r"<[^>]+>", " ", work.get("display_name") or ""))
                    abstract = _text(openalex_abstract(work))
                    source_text = _text(f"{title}. {abstract}")
                    if not source_text:
                        continue
                    scope = (
                        "species_direct"
                        if species.casefold() in source_text.casefold()
                        else "species_indirect"
                    )
                    rows.append(
                        {
                            "accepted_species": species,
                            "source_text": source_text,
                            "source_url": _text(work.get("doi") or work.get("id")),
                            "source_citation": _text(
                                f"{title} {work.get('publication_year') or ''}"
                            ),
                            "source_type": "openalex_title_abstract",
                            "evidence_scope": scope,
                        }
                    )
            except OpenAlexTransientError as exc:
                errors.append(f"openalex:{species}:transient:{exc}")
            except httpx.HTTPStatusError as exc:
                errors.append(f"openalex:{species}:permanent:{exc}")
            except Exception as exc:  # noqa: BLE001
                errors.append(f"openalex:{species}:transient:unexpected:{exc}")
            if pause_seconds > 0:
                time.sleep(pause_seconds)
    return pd.DataFrame(rows, columns=ecology.TEXT_SOURCE_COLUMNS), errors


'''
text, count = fetch_pattern.subn(fetch_replacement, text, count=1)
if count != 1:
    raise SystemExit(f"fetch_openalex_sources replacement count={count}")

apply_pattern = re.compile(
    r"def apply_task_result\(.*?return prepare_dependent_statuses\(result, config\)\n",
    re.DOTALL,
)
apply_replacement = '''def _is_transient_source_error(error: str) -> bool:
    return ":transient:" in str(error)


def apply_task_result(
    ledger: pd.DataFrame,
    batch: pd.DataFrame,
    candidates: pd.DataFrame,
    errors: list[str],
    source_species: set[str],
    task: str,
    wave_id: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    result = ledger.copy()
    max_attempts = int(config.get("max_attempts", 3))
    error_by_species: dict[str, list[str]] = {}
    for error in errors:
        parts = str(error).split(":", 2)
        if len(parts) >= 2:
            error_by_species.setdefault(parts[1], []).append(error)

    target = candidates.loc[candidates["target_for_task"].astype(bool)].copy()
    target_counts = target.groupby("accepted_species").size().to_dict() if not target.empty else {}
    direct_biotic = set(
        candidates.loc[
            candidates["trait_name"].eq("pollen_vector_mode")
            & candidates["candidate_value"].eq("biotic")
            & candidates["evidence_scope"].eq("species_direct"),
            "accepted_species",
        ].astype(str)
    )

    for species in batch["accepted_species"].astype(str):
        mask = result["accepted_species"].eq(species)
        attempts_column = f"{task}_attempts"
        previous_attempts = int(result.loc[mask, attempts_column].iloc[0])
        species_errors = error_by_species.get(species, [])
        transient_failure = bool(species_errors) and all(
            _is_transient_source_error(error) for error in species_errors
        )
        result.loc[mask, f"{task}_wave_id"] = wave_id
        result.loc[mask, f"{task}_candidate_count"] = int(target_counts.get(species, 0))

        if species in source_species or not species_errors:
            result.loc[mask, attempts_column] = previous_attempts + 1
            result.loc[mask, f"{task}_status"] = "processed"
            result.loc[mask, f"{task}_last_error"] = ""
        elif transient_failure:
            # API quota, missing credentials, 429, 5xx, and transport failures are
            # infrastructure states. They must not consume the biological lookup
            # attempt budget or become a terminal evidence absence.
            result.loc[mask, attempts_column] = previous_attempts
            result.loc[mask, f"{task}_status"] = "retry"
            result.loc[mask, f"{task}_last_error"] = " | ".join(species_errors)[:1000]
        else:
            attempts = previous_attempts + 1
            result.loc[mask, attempts_column] = attempts
            status = "exhausted" if attempts >= max_attempts else "retry"
            result.loc[mask, f"{task}_status"] = status
            result.loc[mask, f"{task}_last_error"] = " | ".join(species_errors)[:1000]
        if species in direct_biotic:
            result.loc[mask, "machine_biotic_candidate"] = True
    return prepare_dependent_statuses(result, config)
'''
text, count = apply_pattern.subn(apply_replacement, text, count=1)
if count != 1:
    raise SystemExit(f"apply_task_result replacement count={count}")

path.write_text(text, encoding="utf-8")

config_path = Path("config/global_trait_campaign.yml")
config = config_path.read_text(encoding="utf-8")
config = config.replace(
    "    pause_seconds: 0.05\n    per_page: 5\n",
    "    pause_seconds: 0.25\n    per_page: 5\n    require_api_key: true\n    max_retries: 3\n    backoff_seconds: 5.0\n    max_backoff_seconds: 60.0\n",
    1,
)
config = config.replace(
    "    pause_seconds: 0.05\n    per_page: 8\n",
    "    pause_seconds: 0.25\n    per_page: 8\n    require_api_key: true\n    max_retries: 3\n    backoff_seconds: 5.0\n    max_backoff_seconds: 60.0\n",
    1,
)
config_path.write_text(config, encoding="utf-8")
