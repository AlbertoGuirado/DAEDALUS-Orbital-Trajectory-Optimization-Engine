"""Provenance of experimental results: code version and runtime environment (RR-03).

Every experiment report embeds ``environment_record()`` and
``code_version()`` next to its configuration, so a result can be traced to the
code and the dependency versions that produced it. Nothing here affects a
computation.
"""

from __future__ import annotations

import platform
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np
import scipy

DISTRIBUTION: str = "autonomous-orbital-transfer-optimizer"


def environment_record() -> dict[str, Any]:
    """Python, NumPy, SciPy and platform identifiers of the running interpreter."""
    return {
        "python": sys.version.split()[0],
        "python_implementation": platform.python_implementation(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
    }


def _git(arguments: list[str], cwd: Path) -> str | None:
    try:
        completed = subprocess.run(["git", *arguments], cwd=cwd, capture_output=True, text=True,
                                   timeout=30, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def code_version(repository: Path | str | None = None) -> dict[str, Any]:
    """Package version and, when the source is a Git checkout, its commit.

    Keys: ``package_version`` (installed distribution metadata, or None),
    ``git_commit`` (full hash, or None when Git or the repository is not
    available), ``git_branch`` and ``git_dirty`` (True if tracked files differ
    from the commit; None when unknown). ``repository`` defaults to the
    checkout that contains this file. A dirty tree means the commit alone does
    not identify the code that ran.
    """
    try:
        version: str | None = metadata.version(DISTRIBUTION)
    except metadata.PackageNotFoundError:
        version = None
    root = Path(repository) if repository is not None else Path(__file__).resolve().parents[2]
    commit = _git(["rev-parse", "HEAD"], root)
    branch = _git(["rev-parse", "--abbrev-ref", "HEAD"], root) if commit else None
    status = _git(["status", "--porcelain", "--untracked-files=no"], root) if commit else None
    return {
        "package_version": version,
        "git_commit": commit,
        "git_branch": branch,
        "git_dirty": None if status is None else bool(status),
    }
