#!/usr/bin/env python3
"""Verify the complete SHA256SUMS inventory of this repository."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HEX64 = re.compile(r"[0-9a-f]{64}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def ignored(path: Path) -> bool:
    relative = path.relative_to(ROOT)
    return (
        ".git" in relative.parts
        or "replay-results" in relative.parts
        or "__pycache__" in relative.parts
        or path.name == ".DS_Store"
        or path.suffix == ".pyc"
    )


def main() -> int:
    inventory_path = ROOT / "SHA256SUMS"
    if not inventory_path.is_file():
        raise SystemExit("SHA256SUMS is missing")

    expected: dict[str, str] = {}
    for number, line in enumerate(
        inventory_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        digest, separator, relative = line.partition("  ")
        if separator != "  " or HEX64.fullmatch(digest) is None or not relative:
            raise SystemExit(f"malformed SHA256SUMS line {number}")
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts or relative in expected:
            raise SystemExit(f"unsafe or duplicate checksum path on line {number}")
        expected[relative] = digest

    observed_paths: set[str] = set()
    for path in ROOT.rglob("*"):
        if ignored(path):
            continue
        if path.is_symlink():
            raise SystemExit(f"symlink is not allowed: {path.relative_to(ROOT)}")
        if path.is_file() and path.name != "SHA256SUMS":
            observed_paths.add(path.relative_to(ROOT).as_posix())

    if observed_paths != set(expected):
        missing = sorted(set(expected) - observed_paths)
        extra = sorted(observed_paths - set(expected))
        raise SystemExit(f"checksum inventory mismatch: missing={missing}, extra={extra}")

    for relative, digest in sorted(expected.items()):
        observed = sha256_file(ROOT / relative)
        if observed != digest:
            raise SystemExit(
                f"checksum mismatch for {relative}: expected {digest}, observed {observed}"
            )

    print(f"PASS: SHA256SUMS covers {len(expected)} repository files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
