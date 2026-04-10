#!/usr/bin/env python3
"""Check if LAMMPS atom->sortme() reorders atoms during run 0."""

import numpy as np
from pathlib import Path
from lammps import lammps

data_path = Path("ipi/lammps_data_ipi_client.data").resolve()
lmp = lammps(cmdargs=["-nocite", "-log", "none", "-echo", "none"])

lmp.commands_list([
    "units metal",
    "atom_style atomic",
    "boundary p p p",
    "pair_style zero 1.0",
    f"read_data {data_path}",
    "pair_coeff * *",
])

# Before run 0
tags_before = np.array(lmp.numpy.extract_atom("id")[:192], dtype=int, copy=True)
types_before = np.array(lmp.numpy.extract_atom("type")[:192], dtype=int, copy=True)
pos_before = np.array(lmp.numpy.extract_atom("x")[:192], dtype=np.float64, copy=True)

print("Before run 0:")
print(f"  tags[0:15]: {tags_before[:15].tolist()}")
print(f"  types[0:15]: {types_before[:15].tolist()}")

# Do a run 0 (triggers setup → sortme)
lmp.command("run 0")

# After run 0
tags_after = np.array(lmp.numpy.extract_atom("id")[:192], dtype=int, copy=True)
types_after = np.array(lmp.numpy.extract_atom("type")[:192], dtype=int, copy=True)
pos_after = np.array(lmp.numpy.extract_atom("x")[:192], dtype=np.float64, copy=True)

print("\nAfter run 0:")
print(f"  tags[0:15]: {tags_after[:15].tolist()}")
print(f"  types[0:15]: {types_after[:15].tolist()}")

# Compare
tags_match = np.array_equal(tags_before, tags_after)
print(f"\nTags match before/after: {tags_match}")

if not tags_match:
    n_changed = np.sum(tags_before != tags_after)
    print(f"  {n_changed} of {len(tags_before)} atoms changed position in internal array")
    first_changes = np.where(tags_before != tags_after)[0][:10]
    for idx in first_changes:
        print(f"    idx={idx}: tag {tags_before[idx]} (type {types_before[idx]}) "
              f"→ tag {tags_after[idx]} (type {types_after[idx]})")

    # Check: are the same atoms present?
    same_set = set(tags_before.tolist()) == set(tags_after.tolist())
    print(f"  Same atom set: {same_set}")

    # Position diff: comparing by index (wrong order) vs by tag (right order)
    pos_diff_idx = np.abs(pos_before - pos_after).max()
    print(f"  Max pos diff by index: {pos_diff_idx:.4f} Å (expect large if sorted)")

    # Compare by tag
    before_map = {t: pos_before[i] for i, t in enumerate(tags_before)}
    after_map = {t: pos_after[i] for i, t in enumerate(tags_after)}
    pos_diff_tag = max(np.abs(before_map[t] - after_map[t]).max() for t in before_map)
    print(f"  Max pos diff by tag: {pos_diff_tag:.2e} Å (should be ~0 if no wrapping)")

# Now test: disable sorting
print("\n=== Test with atom_modify sort 0 ===")
lmp2 = lammps(cmdargs=["-nocite", "-log", "none", "-echo", "none"])
lmp2.commands_list([
    "units metal",
    "atom_style atomic",
    "atom_modify sort 0 0.0",
    "boundary p p p",
    "pair_style zero 1.0",
    f"read_data {data_path}",
    "pair_coeff * *",
])

tags2_before = np.array(lmp2.numpy.extract_atom("id")[:192], dtype=int, copy=True)
lmp2.command("run 0")
tags2_after = np.array(lmp2.numpy.extract_atom("id")[:192], dtype=int, copy=True)
print(f"  Tags match (sort disabled): {np.array_equal(tags2_before, tags2_after)}")

lmp.close()
lmp2.close()
