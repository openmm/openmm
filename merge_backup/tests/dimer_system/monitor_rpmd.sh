#!/bin/bash
# Monitor RPMD simulation progress

echo "=== RPMD Simulation Monitor ==="
echo ""

# Check if process is running
if ps aux | grep -q "[r]un_simulation_rpmd.py"; then
    echo "✓ Simulation is RUNNING"
    PID=$(ps aux | grep "[r]un_simulation_rpmd.py" | awk '{print $2}')
    echo "  PID: $PID"
    
    # Get CPU and memory usage
    CPU=$(ps aux | grep "[r]un_simulation_rpmd.py" | awk '{print $3}')
    MEM=$(ps aux | grep "[r]un_simulation_rpmd.py" | awk '{print $4}')
    echo "  CPU: ${CPU}%"
    echo "  MEM: ${MEM}%"
else
    echo "✗ Simulation is NOT running"
fi

echo ""
echo "=== Latest Log Output ==="
if [ -f rpmd_simulation.log ]; then
    tail -30 rpmd_simulation.log
else
    echo "No log file found"
fi

echo ""
echo "=== Output File Status ==="
if [ -f cavity_diamer_lambda0.0700.npz ]; then
    ls -lh cavity_diamer_lambda0.0700.npz
    python3 -c "
import numpy as np
try:
    data = np.load('cavity_diamer_lambda0.0700.npz', allow_pickle=True)
    metadata = data['metadata'].item()
    print(f\"  Current step: {metadata.get('step', 'unknown')}/{metadata.get('total_steps', 'unknown')}\")
    print(f\"  Status: {metadata.get('status', 'unknown')}\")
    print(f\"  Data points: {len(data['time_ps'])}\")
except Exception as e:
    print(f\"  Error reading file: {e}\")
" 2>&1
else
    echo "No output file yet"
fi
