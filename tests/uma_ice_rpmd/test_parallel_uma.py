# test_parallel_uma.py - Testing TRUE TENSOR BATCHING
import sys
import time
import torch
import numpy as np
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/fairchem/src')

from fairchem.core import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch
from ase import Atoms

# Create test data (30 atoms, fixed topology)
symbols = ['O', 'H', 'H'] * 10
base_positions = np.random.rand(30, 3) * 10
cell = np.eye(3) * 10

print("="*80)
print("TESTING: TRUE TENSOR BATCHING with UMA")
print("="*80)
print("\nSetup: 30-atom water system, same topology, different positions")
print("Goal: Compare sequential predictions vs batched graph input\n")

# Initialize model
model = pretrained_mlip.get_predict_unit('uma-s-1p1', device='cuda')

# Test 1: Sequential predictions (baseline)
print("[1] SEQUENTIAL: 8 predictions one-by-one")
print("-" * 80)

sequential_data_list = []
for i in range(8):
    # Create slightly perturbed positions for each "bead"
    positions = base_positions + np.random.randn(30, 3) * 0.1
    atoms = Atoms(symbols=symbols, positions=positions, pbc=True, cell=cell)
    atoms.info['charge'] = 0
    atoms.info['spin'] = 0
    # r_edges=True to compute edges before batching
    # max_neigh needs to be specified for CPU graph construction
    data = AtomicData.from_ase(atoms, task_name='omol', r_edges=True, radius=6.0, max_neigh=50, r_data_keys=['spin', 'charge']).to('cuda')
    data.dataset = ['omol']
    sequential_data_list.append(data)

# Warmup
_ = model.predict(sequential_data_list[0])
torch.cuda.synchronize()

# Time sequential predictions
seq_times = []
seq_results = []
for i, data in enumerate(sequential_data_list):
    t0 = time.time()
    result = model.predict(data)
    torch.cuda.synchronize()
    elapsed = (time.time() - t0) * 1000
    seq_times.append(elapsed)
    seq_results.append(result)

seq_total = sum(seq_times)
print(f"  Individual times: {[f'{t:.1f}ms' for t in seq_times]}")
print(f"  Average per prediction: {np.mean(seq_times):.1f} ms")
print(f"  Total time: {seq_total:.1f} ms")

# Test 2: Graph batching (all 8 systems in one batch)
print("\n[2] GRAPH BATCHING: All 8 systems in one batch")
print("-" * 80)

# Create batched input using FAIRChem's custom batching function
batch_data = atomicdata_list_to_batch(sequential_data_list)

# Fix dataset attribute - it should be a list of strings, not a list of lists
# After batching, dataset becomes [['omol'], ['omol'], ...], but we need ['omol', 'omol', ...]
if isinstance(batch_data.dataset, list) and len(batch_data.dataset) > 0:
    if isinstance(batch_data.dataset[0], list):
        batch_data.dataset = [d[0] if isinstance(d, list) else d for d in batch_data.dataset]
print(f"  Batch structure:")
print(f"    - Total atoms: {batch_data.pos.shape[0]} (should be 240)")
print(f"    - Batch dimension: {batch_data.batch.shape}")
print(f"    - Edge indices shape: {batch_data.edge_index.shape}")

# Warmup
_ = model.predict(batch_data)
torch.cuda.synchronize()

# Time batched prediction
batch_times = []
for _ in range(3):  # Multiple runs for average
    t0 = time.time()
    batch_result = model.predict(batch_data)
    torch.cuda.synchronize()
    elapsed = (time.time() - t0) * 1000
    batch_times.append(elapsed)

batch_avg = np.mean(batch_times)
print(f"  Batch prediction times: {[f'{t:.1f}ms' for t in batch_times]}")
print(f"  Average: {batch_avg:.1f} ms")

# Analysis
print("\n" + "="*80)
print("RESULTS")
print("="*80)
print(f"Sequential (8 separate predictions): {seq_total:.1f} ms")
print(f"Batched (1 graph with 8 systems):     {batch_avg:.1f} ms")
print(f"Speedup: {seq_total / batch_avg:.2f}x")
print(f"Efficiency: {(seq_total / batch_avg) / 8 * 100:.1f}% of ideal 8x speedup")

if batch_avg < seq_total * 0.3:
    print("\nEXCELLENT batching efficiency! >70% of ideal speedup")
elif batch_avg < seq_total * 0.5:
    print("\nGOOD batching efficiency! ~50% of ideal speedup")
elif batch_avg < seq_total * 0.7:
    print("\nMODERATE batching efficiency. ~30-50% of ideal speedup")
else:
    print("\nPOOR batching efficiency. Graph batching not helping much")

# Test 3: Verify correctness
print("\n[3] CORRECTNESS CHECK")
print("-" * 80)
print("Comparing batched vs sequential forces for first system...")

# Extract forces for first system from batched result
batch_forces = batch_result['forces'][:30].detach().numpy()#cpu().numpy()  # First 30 atoms
seq_forces = seq_results[0]['forces'].detach().numpy()# cpu().numpy()

force_diff = np.abs(batch_forces - seq_forces).max()
print(f"  Max force difference: {force_diff:.2e} eV/Å")

if force_diff < 1e-5:
    print("  Forces match perfectly!")
elif force_diff < 1e-3:
    print("  Forces match within tolerance")
else:
    print(f"  WARNING: Forces differ by {force_diff:.2e} - potential issue!")

# Additional analysis: Check memory usage
print("\n[4] MEMORY ANALYSIS")
print("-" * 80)
print(f"  GPU memory allocated: {torch.cuda.memory_allocated() / 1024**2:.0f} MB")
print(f"  GPU memory reserved:  {torch.cuda.memory_reserved() / 1024**2:.0f} MB")
