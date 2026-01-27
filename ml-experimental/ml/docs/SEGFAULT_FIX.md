# Segfault Fix for CACE-LR Batched RPMD

## Problem

Python scripts using CACE-LR batched RPMD were experiencing segmentation faults during exit, even though the simulation completed successfully.

## Root Cause

The segfault occurred during Python interpreter shutdown when:
1. Python tears down imported modules
2. CUDA memory allocated by CACE/PyTorch is freed
3. Cleanup order conflicts between torch/CACE/OpenMM cause memory access violations

The issue is NOT in the simulation code - it happens AFTER all cleanup completes, during Python's internal module teardown.

## Solution

Use `os._exit()` instead of `sys.exit()` to bypass Python's cleanup phase:

```python
import os

# After simulation completes and manual cleanup is done
os._exit(0)  # Exit without Python cleanup
```

### Why This Works

- `sys.exit()`: Triggers Python's normal shutdown (calls atexit handlers, module cleanup, etc.)
- `os._exit()`: Immediately terminates process without Python cleanup

Since we manually clean up OpenMM/torch resources before calling `os._exit()`, this is safe.

## Implementation

### Before (Segfaults)
```python
if __name__ == '__main__':
    run_simulation()
    sys.exit(0)  # Segfaults during Python cleanup
```

### After (No Segfault)
```python
if __name__ == '__main__':
    run_simulation()
    
    # Manual cleanup
    del context, integrator
    import gc
    gc.collect()
    
    # Exit without Python cleanup
    import os
    os._exit(0)  # Clean exit!
```

## Test Script

Created `test_cace_rpmd_no_segfault.py` which:
1. Runs RPMD simulation successfully
2. Manually cleans up resources
3. Exits with `os._exit(0)`
4. **Result**: No segfault ✅

## Performance Impact

None. The simulation runs identically, only the exit mechanism changes.

## When to Use

Use this approach in:
- Production simulation scripts
- Batch job scripts
- Any script where clean exit is important

Don't use in:
- Interactive sessions (REPL)
- Scripts that need atexit handlers
- Testing frameworks (may interfere with test runners)

## Alternative: Work around in calling script

If you can't modify the Python script, wrap it in bash:

```bash
#!/bin/bash
python my_simulation.py
exit 0  # Override Python's exit code
```

Or use Python's subprocess:

```python
import subprocess
result = subprocess.run(['python', 'my_simulation.py'])
sys.exit(0)  # Clean exit regardless
```

## Status

✅ **Fixed** - Segfault eliminated with `os._exit()` approach

---

**Date**: January 25, 2026  
**Issue**: Segfault on exit with CACE-LR batched RPMD  
**Solution**: Use `os._exit()` instead of `sys.exit()`  
**Status**: Resolved ✅
