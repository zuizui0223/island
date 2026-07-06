# Global trait batch failure recovery

The 100-taxon extraction workflow treats a failed API request as a recoverable
acquisition failure, not as a reason to discard the worklist or completed
sub-batches.

Each run now:

1. submits the first taxon as a real preflight request;
2. uses a bounded web-search output budget;
3. retries only transient errors;
4. snapshots raw responses, candidate rows, and a manifest after every completed
   sub-batch; and
5. uploads the artifact even when a later sub-batch fails.

The manifest records the failed batch, exception class, status code when
available, and a redacted error message. The workflow then reports failure after
artifact upload so the operator can inspect the diagnostics and resume from the
same offset without losing successful candidate rows.
