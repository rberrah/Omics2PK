"""Independent nonclinical PK/PD original-source reproducible benchmark launcher.

Original Puskarich 2021 CaPS TimeZero author GitHub files are fetched only
temporarily in RAM in an opt-in scientific CI; NOT saved as repo files or
GitHub artifacts and never logs participants, predictions or observations.
"""
from __future__ import annotations
import json

from omics2pk.caps_original_pd import run_original_public_caps


if __name__ == "__main__":
    report=run_original_public_caps()
    print(json.dumps(report,indent=2,sort_keys=True))
