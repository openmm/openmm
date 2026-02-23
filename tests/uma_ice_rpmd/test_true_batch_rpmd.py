# test_true_batch_rpmd.py - Testing TRUE RPMD-STYLE BATCHING
"""
This tests if UMA can handle multiple INDEPENDENT systems batched together,
where each system gets its own separate energy and forces.
"""
import sys
import time
import torch
import numpy as np
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/fairchem/src')

from fairchem.core import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch
from fairchem.core.graph.compute import generate_graph
from ase import Atoms

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
torch.backends.cudnn.benchmark = False
#torch.use_deterministic_algorithms(True)

def create_batched_data_rpmd_style(positions_list, symbols, cell, charge=0, spin=0, task_name='omol'):
    """
    Create a properly batched AtomicData for RPMD where each system is independent.
    
    Args:
        positions_list: List of position arrays, one per bead
        symbols: Atomic symbols (same for all beads)
        cell: Unit cell (same for all beads)
        
    Returns:
        AtomicData object with proper batch indices for independent systems
    """
    num_beads = len(positions_list)
    num_atoms_per_bead = len(symbols)
    
    data_list = []
    for i, positions in enumerate(positions_list):
        atoms = Atoms(symbols=symbols, positions=positions, pbc=True, cell=cell)
        atoms.info['charge'] = charge
        atoms.info['spin'] = spin
        data = AtomicData.from_ase(
            atoms,
            task_name=task_name,
            r_edges=False,
            r_data_keys=['spin', 'charge'],
        )
        # Keep dataset as a string so batching yields list[str], not list[list[str]]
        data.sid = [f"bead-{i}"]
        data_list.append(data)

    return atomicdata_list_to_batch(data_list)


print("="*80)
print("TESTING: TRUE RPMD-STYLE BATCHING with UMA")
print("="*80)
print("\nSetup: 8 independent 30-atom systems (RPMD beads)")
print("Goal: Get 8 separate energies and force arrays\n")

# Create test data (30 atoms)
symbols = ['O', 'H', 'H'] * 10
base_positions = np.random.rand(30, 3) * 10
cell = np.eye(3) * 10

# Initialize model
model = pretrained_mlip.get_predict_unit('uma-s-1p1', device='cuda')

# Test 1: Sequential predictions (baseline)
print("[1] SEQUENTIAL: 8 predictions one-by-one")
print("-" * 80)

positions_list = []
sequential_results = []
for i in range(8):
    positions = base_positions + np.random.randn(30, 3) * 0.1
    positions_list.append(positions)
    
    atoms = Atoms(symbols=symbols, positions=positions, pbc=True, cell=cell)
    atoms.info['charge'] = 0
    atoms.info['spin'] = 0
    data = AtomicData.from_ase(atoms, task_name='omol', r_edges=False, r_data_keys=['spin', 'charge']).to('cuda')
    # For single-graph inference, dataset must be a list to avoid iterating chars
    data.dataset = ['omol']
    sequential_results.append(data)

# Warmup
_ = model.predict(sequential_results[0])
torch.cuda.synchronize()

# Time sequential
seq_times = []
seq_energies = []
seq_forces_list = []
for data in sequential_results:
    t0 = time.time()
    result = model.predict(data)
    torch.cuda.synchronize()
    elapsed = (time.time() - t0) * 1000
    seq_times.append(elapsed)
    seq_energies.append(float(result['energy'].cpu().item()))
    seq_forces_list.append(result['forces'].detach().cpu().numpy())

seq_total = sum(seq_times)
print(f"  Individual times: {[f'{t:.1f}ms' for t in seq_times]}")
print(f"  Average per prediction: {np.mean(seq_times):.1f} ms")
print(f"  Total time: {seq_total:.1f} ms")
print(f"  Energies: {[f'{e:.2f}' for e in seq_energies[:3]]}... eV")

# Test 2: RPMD-style batching (all 8 as independent systems)
print("\n[2] RPMD-STYLE BATCHING: 8 independent systems in one call")
print("-" * 80)

# Create properly batched data
batch_data = create_batched_data_rpmd_style(positions_list, symbols, cell)
batch_data = batch_data.to('cuda')

print(f"  Batch structure:")
print(f"    - Total atoms: {batch_data.pos.shape[0]} (should be 240)")
print(f"    - Batch indices shape: {batch_data.batch.shape}")
print(f"    - Num graphs: {batch_data.num_graphs} (should be 8)")
print(f"    - natoms: {batch_data.natoms.tolist()} (should be [30]*8)")
print(f"    - Edge indices: {batch_data.edge_index.shape} (empty, will be generated)")

# Warmup
_ = model.predict(batch_data)
torch.cuda.synchronize()

# Time batched prediction
batch_times = []
for _ in range(3):
    t0 = time.time()
    batch_result = model.predict(batch_data)
    torch.cuda.synchronize()
    elapsed = (time.time() - t0) * 1000
    batch_times.append(elapsed)

batch_avg = np.mean(batch_times)
print(f"\n  Batch prediction times: {[f'{t:.1f}ms' for t in batch_times]}")
print(f"  Average: {batch_avg:.1f} ms")

# Extract results
batch_energies = batch_result['energy'].detach().cpu().numpy()  # Should be [8] separate energies
batch_forces = batch_result['forces'].detach().cpu().numpy()    # Should be [240, 3]

# Split forces by system using batch indices to preserve bead order
batch_indices = batch_data.batch.detach().cpu().numpy()
batch_forces_split = [batch_forces[batch_indices == i] for i in range(8)]

print(f"\n  Energies shape: {batch_energies.shape} (should be [8])")
print(f"  Energies: {[f'{e:.2f}' for e in batch_energies[:3]]}... eV")
print(f"  Forces shape: {batch_forces.shape} (should be [240, 3])")

# Analysis
print("\n" + "="*80)
print("RESULTS")
print("="*80)
print(f"Sequential (8 separate predictions): {seq_total:.1f} ms")
print(f"Batched (8 independent systems):      {batch_avg:.1f} ms")
print(f"Speedup: {seq_total / batch_avg:.2f}x")
print(f"Efficiency: {(seq_total / batch_avg) / 8 * 100:.1f}% of ideal 8x speedup")

if batch_avg < seq_total * 0.3:
    print("\nEXCELLENT batching efficiency! >70% of ideal speedup")
elif batch_avg < seq_total * 0.5:
    print("\nGOOD batching efficiency! ~50% of ideal speedup")
elif batch_avg < seq_total * 0.7:
    print("\nMODERATE batching efficiency. ~30-50% of ideal speedup")
else:
    print("\nPOOR batching efficiency. Batching not helping much")

# Test 3: CORRECTNESS CHECK
print("\n" + "="*80)
print("CORRECTNESS CHECK")
print("="*80)

# Check if we got 8 separate energies
if batch_energies.shape[0] == 8:
    print("Got 8 separate energies!")
else:
    print(f"Expected 8 energies, got {batch_energies.shape[0]}")

# Compare first system's results
print("\nPer-bead energies (sequential vs batched):")
for i in range(8):
    energy_diff = abs(seq_energies[i] - batch_energies[i])
    print(
        f"  bead {i}: seq={seq_energies[i]: .6f} eV | "
        f"batch={batch_energies[i]: .6f} eV | diff={energy_diff:.2e}"
    )

print("\nPer-bead force max diff (sequential vs batched):")
force_diffs = []
for i in range(8):
    diff = np.abs(seq_forces_list[i] - batch_forces_split[i]).max()
    force_diffs.append(diff)
    print(f"  bead {i}: max |ΔF| = {diff:.2e} eV/Å")

print("\nGraph diagnostics (edge counts per bead):")
backbone = model.model.module.backbone
cutoff = backbone.cutoff
max_neighbors = backbone.max_neighbors
enforce_strict = backbone.enforce_max_neighbors_strictly
radius_pbc_version = backbone.radius_pbc_version
for i in range(8):
    seq_data = sequential_results[i].cpu()
    batch_data_i = batch_data.get_example(i).cpu()
    seq_graph = generate_graph(
        seq_data,
        cutoff=cutoff,
        max_neighbors=max_neighbors,
        enforce_max_neighbors_strictly=enforce_strict,
        radius_pbc_version=radius_pbc_version,
        pbc=seq_data.pbc,
    )
    batch_graph = generate_graph(
        batch_data_i,
        cutoff=cutoff,
        max_neighbors=max_neighbors,
        enforce_max_neighbors_strictly=enforce_strict,
        radius_pbc_version=radius_pbc_version,
        pbc=batch_data_i.pbc,
    )
    seq_edges = seq_graph["edge_index"].shape[1]
    batch_edges = batch_graph["edge_index"].shape[1]
    print(f"  bead {i}: seq_edges={seq_edges} | batch_edges={batch_edges}")

    if seq_edges != batch_edges:
        print(f"    WARNING: edge count mismatch for bead {i}")

def _edge_set(edge_index):
    edges = edge_index.cpu().numpy().T
    edges = edges[np.lexsort((edges[:, 1], edges[:, 0]))]
    return edges

print("\nGraph diagnostics (edge identity check for bead 0):")
seq_data0 = sequential_results[0].cpu()
batch_data0 = batch_data.get_example(0).cpu()
seq_graph0 = generate_graph(
    seq_data0,
    cutoff=cutoff,
    max_neighbors=max_neighbors,
    enforce_max_neighbors_strictly=enforce_strict,
    radius_pbc_version=radius_pbc_version,
    pbc=seq_data0.pbc,
)
batch_graph0 = generate_graph(
    batch_data0,
    cutoff=cutoff,
    max_neighbors=max_neighbors,
    enforce_max_neighbors_strictly=enforce_strict,
    radius_pbc_version=radius_pbc_version,
    pbc=batch_data0.pbc,
)
seq_edges0 = _edge_set(seq_graph0["edge_index"])
batch_edges0 = _edge_set(batch_graph0["edge_index"])
same_edges0 = seq_edges0.shape == batch_edges0.shape and np.array_equal(seq_edges0, batch_edges0)
print(f"  bead 0 edge sets identical: {same_edges0}")
if not same_edges0:
    print(f"    seq_edges={seq_edges0.shape[0]} batch_edges={batch_edges0.shape[0]}")

print("\nBead matching matrix (max |ΔF|) between batched i and sequential j):")
match_matrix = np.zeros((8, 8))
for i in range(8):
    for j in range(8):
        match_matrix[i, j] = np.abs(batch_forces_split[i] - seq_forces_list[j]).max()

for i in range(8):
    row = "  " + " ".join([f"{match_matrix[i, j]:.2e}" for j in range(8)])
    print(f"  batched {i} -> {row}")

best_matches = np.argmin(match_matrix, axis=1)
print("\nBest sequential match for each batched bead:")
for i in range(8):
    print(
        f"  batched {i} best matches sequential {best_matches[i]} "
        f"(max |ΔF| = {match_matrix[i, best_matches[i]]:.2e} eV/Å)"
    )

energy_diff0 = abs(seq_energies[0] - batch_energies[0])
force_diff0 = np.abs(seq_forces_list[0] - batch_forces_split[0]).max()
if energy_diff0 < 1e-4 and force_diff0 < 1e-4:
    print("\nPERFECT! Batched results match sequential!")
    print("TRUE RPMD-STYLE BATCHING IS WORKING!")
elif energy_diff0 < 1e-2 and force_diff0 < 1e-2:
    print("\nForces match within tolerance")
else:
    print("\nResults differ significantly - something is wrong")

# Memory
print("\n" + "="*80)
print("MEMORY ANALYSIS")
print("="*80)
print(f"  GPU memory allocated: {torch.cuda.memory_allocated() / 1024**2:.0f} MB")
print(f"  GPU memory reserved:  {torch.cuda.memory_reserved() / 1024**2:.0f} MB")
