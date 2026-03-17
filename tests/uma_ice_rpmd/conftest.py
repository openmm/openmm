"""Pytest configuration for uma_ice_rpmd tests."""
import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_lammps: mark test as requiring LAMMPS + UMA (fairchem-lammps)",
    )
