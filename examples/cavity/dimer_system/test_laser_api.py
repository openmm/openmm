#!/usr/bin/env python3
"""
Simple test to verify laser API is available after compilation.
This tests the Python bindings without full OpenMM import.
"""
import sys
import os

# Add build directory to path
build_python = os.path.join(os.path.dirname(__file__), 'build', 'python')
if os.path.exists(build_python):
    sys.path.insert(0, build_python)

try:
    # Try to import just the SWIG module directly
    import _openmm
    print("_openmm module loaded")
    
    # Test if CavityForce has the new methods
    cf = _openmm.CavityForce(0, 0.01, 0.001, 1.0)
    print("CavityForce created")
    
    # Test new methods exist
    methods_to_check = [
        'setCavityDriveAmplitude',
        'getCavityDriveAmplitude',
        'setCavityDriveFrequency',
        'getCavityDriveFrequency',
        'setCavityDrivePhase',
        'getCavityDrivePhase',
        'setCavityDriveEnvelope',
        'getCavityDriveEnvelopeType',
        'setCavityDriveEnabled',
        'getCavityDriveEnabled',
        'setDirectLaserAmplitude',
        'getDirectLaserAmplitude',
        'setDirectLaserFrequency',
        'getDirectLaserFrequency',
        'setDirectLaserPhase',
        'getDirectLaserPhase',
        'setDirectLaserEnvelope',
        'getDirectLaserEnvelopeType',
        'setDirectLaserCouplingEnabled',
        'getDirectLaserCouplingEnabled',
        'getCavityDriveEnergy',
        'getDirectLaserEnergy'
    ]
    
    missing = []
    for method in methods_to_check:
        if not hasattr(cf, method):
            missing.append(method)
    
    if missing:
        print(f"FAIL: Missing methods: {missing}")
        sys.exit(1)
    else:
        print(f"All {len(methods_to_check)} laser methods are available!")
        
    # Test setting values
    cf.setCavityDriveAmplitude(0.01)
    cf.setCavityDriveFrequency(0.01)
    cf.setCavityDriveEnvelope("gaussian", 25.0, 10.0)
    cf.setCavityDriveEnabled(True)
    
    cf.setDirectLaserAmplitude(0.005)
    cf.setDirectLaserFrequency(0.00913)
    cf.setDirectLaserEnvelope("square", 0.0, 50.0)
    cf.setDirectLaserCouplingEnabled(False)
    
    print("All setter methods work!")
    print("\nSUCCESS: Laser driving API is fully implemented and accessible!")
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
