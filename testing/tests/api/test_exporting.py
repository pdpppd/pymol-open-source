import os
import pymol
import sys
from pymol import cmd, testing
from pymol import test_utils
import tempfile

import pytest


@test_utils.requires_version("3.2")
def test_bcif_export():
    """Test BCIF export and round-trip"""
    # Create a simple structure
    cmd.fragment("ala")
    orig_count = cmd.count_atoms("ala")
    assert orig_count == 10
    # Export to BCIF
    with tempfile.NamedTemporaryFile(suffix='.bcif', delete=False) as f:
        bcif_file = f.name

    try:
        cmd.save(bcif_file, "ala")
        assert os.path.exists(bcif_file)
        assert os.path.getsize(bcif_file) > 0

        # Load back and verify
        cmd.delete("all")
        cmd.load(bcif_file, "test_loaded")
        loaded_count = cmd.count_atoms("test_loaded")
        assert loaded_count == orig_count, f"Atom count mismatch: {loaded_count} != {orig_count}"
    finally:
        if os.path.exists(bcif_file):
            os.unlink(bcif_file)
