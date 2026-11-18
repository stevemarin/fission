"""Shared helpers for building small ENDF-6 sections in tests."""


def endf_line(fields, mat=125, mf=3, mt=1):
    """Build an 80-column ENDF line from up to six already-formatted fields."""
    data = "".join(f"{f:>11}" for f in fields).ljust(66)
    return f"{data}{mat:4d}{mf:2d}{mt:3d}{0:5d}"
