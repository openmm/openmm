#!/usr/bin/env python3
"""
Master Test Runner for UMA/RPMD Implementation

This script runs all test suites for the UMA/RPMD implementation and provides
a full report of the results.

Test Categories:
1. Regression Tests - Verify specific bugs are fixed
2. Force Accumulation Tests - Verify forces combine correctly
3. Hybrid RPMD Tests - Verify classical/quantum particle handling
4. Stress Tests - Test robustness and performance
5. Edge Cases - Test corner cases and numerical stability
6. Quantum Statistics - Validate physical correctness

Usage:
    python run_all_uma_rpmd_tests.py [--quick] [--category CATEGORY]
    
    --quick: Run only fast tests (skip long simulations)
    --category: Run only specific category (regression, forces, hybrid, stress, edge, statistics)
"""

import sys
import os
import time
import argparse
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

def run_test_file(test_file, timeout=300):
    """
    Run a test file and return whether it passed.
    
    Returns: (passed: bool, duration: float, error_msg: str)
    """
    test_path = Path(__file__).parent / test_file
    
    if not test_path.exists():
        return (None, 0.0, f"Test file not found: {test_file}")
    
    print(f"\n{'='*80}")
    print(f"Running: {test_file}")
    print(f"{'='*80}")
    
    start_time = time.time()
    
    try:
        # Import and run the test
        import importlib.util
        spec = importlib.util.spec_from_file_location("test_module", test_path)
        test_module = importlib.util.module_from_spec(spec)
        
        # Redirect to capture output if needed
        spec.loader.exec_module(test_module)
        
        duration = time.time() - start_time
        
        # Most test modules return success via sys.exit or return from main
        # We'll consider it passed if no exception was raised
        return (True, duration, None)
        
    except SystemExit as e:
        duration = time.time() - start_time
        passed = (e.code == 0)
        return (passed, duration, None if passed else f"Exit code: {e.code}")
        
    except Exception as e:
        duration = time.time() - start_time
        return (False, duration, str(e))


def main():
    parser = argparse.ArgumentParser(description='Run UMA/RPMD test suite')
    parser.add_argument('--quick', action='store_true', help='Run only fast tests')
    parser.add_argument('--category', choices=['regression', 'forces', 'hybrid', 'stress', 'edge', 'statistics', 'all'],
                       default='all', help='Test category to run')
    args = parser.parse_args()
    
    # Define test suites
    test_suites = {
        'regression': [
            'test_bug_regression.py',
        ],
        'forces': [
            'test_uma_force_accumulation.py',
            'test_batch_consistency.py',
        ],
        'hybrid': [
            'test_classical_centroid_force.py',
            'test_uma_hybrid_rpmd_cavity.py',
            'test_hybrid_rpmd_cavity.py',  # Original test
        ],
        'stress': [
            'test_uma_rpmd_stress.py',
        ],
        'edge': [
            'test_edge_cases.py',
        ],
        'statistics': [
            'test_quantum_statistics.py',
        ],
    }
    
    # Select tests to run
    if args.category == 'all':
        selected_tests = []
        for category_tests in test_suites.values():
            selected_tests.extend(category_tests)
    else:
        selected_tests = test_suites.get(args.category, [])
    
    if args.quick:
        # Skip stress tests and long statistics tests in quick mode
        selected_tests = [t for t in selected_tests if 'stress' not in t and 'statistics' not in t]
    
    print("="*80)
    print("UMA/RPMD Test Suite")
    print("="*80)
    print(f"Running {len(selected_tests)} test files")
    print(f"Category: {args.category}")
    print(f"Quick mode: {args.quick}")
    print()
    
    results = []
    total_duration = 0.0
    
    for test_file in selected_tests:
        passed, duration, error = run_test_file(test_file)
        results.append((test_file, passed, duration, error))
        total_duration += duration
        
        if passed:
            status = "PASS"
        elif passed is None:
            status = "⊘ SKIP"
        else:
            status = " FAIL"
        
        print(f"\n{status}: {test_file} ({duration:.2f}s)")
        if error:
            print(f"  Error: {error}")
    
    # Print summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    passed_count = sum(1 for _, p, _, _ in results if p is True)
    failed_count = sum(1 for _, p, _, _ in results if p is False)
    skipped_count = sum(1 for _, p, _, _ in results if p is None)
    
    for test_file, passed, duration, error in results:
        if passed:
            status = ""
        elif passed is None:
            status = "⊘"
        else:
            status = ""
        print(f"{status} {test_file:45s} ({duration:.2f}s)")
    
    print(f"\n{'='*80}")
    print(f"Total: {len(results)} tests")
    print(f"  Passed:  {passed_count}")
    print(f"  Failed:  {failed_count}")
    print(f"  Skipped: {skipped_count}")
    print(f"Total time: {total_duration:.2f}s")
    print(f"{'='*80}")
    
    if failed_count == 0:
        print("\nALL TESTS PASSED")
        return 0
    else:
        print(f"\n {failed_count} TEST(S) FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
