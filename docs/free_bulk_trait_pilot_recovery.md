# Free bulk trait pilot recovery

This PR exists to run the AusTraits all-dataset pilot as a normal pull-request status check after PR #119 was merged before its experiment could be observed on a live PR check.

The check downloads the current AusTraits release, audits all datasets against the 106,295-species island master, and uploads the coverage report as an artifact.
