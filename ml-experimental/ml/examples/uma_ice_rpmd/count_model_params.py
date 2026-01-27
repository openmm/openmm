#!/usr/bin/env python3
"""
Script to determine the actual parameter counts for UMA/eSEN models (skip uma-m-1p1).
"""

import torch
import sys
import gc

def count_parameters(model):
    """Count total and trainable parameters in a model."""
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    return total_params, trainable_params

def format_number(num):
    """Format large numbers in a readable way."""
    if num >= 1e9:
        return f"{num/1e9:.2f}B"
    elif num >= 1e6:
        return f"{num/1e6:.2f}M"
    elif num >= 1e3:
        return f"{num/1e3:.2f}K"
    else:
        return str(num)

def analyze_model(model_name):
    """Analyze a single model and return its parameter counts."""
    from fairchem.core.calculate import pretrained_mlip
    
    try:
        # Load the model on CPU to save memory
        predict_unit = pretrained_mlip.get_predict_unit(model_name, device="cpu")
        model = predict_unit.model
        
        # Count parameters
        total_params, trainable_params = count_parameters(model)
        
        # Try to get active parameters if available
        active_params = None
        if hasattr(model, 'num_params'):
            active_params = model.num_params
        elif hasattr(model, 'active_params'):
            active_params = model.active_params
        
        result = {
            'total': total_params,
            'trainable': trainable_params,
            'active': active_params,
            'success': True
        }
        
        # Clean up
        del predict_unit
        del model
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        return result
        
    except Exception as e:
        return {'success': False, 'error': str(e)}

if __name__ == '__main__':
    models = [
        ('uma-s-1', 'UMA Small v1.0'),
        ('uma-s-1p1', 'UMA Small v1.1'),
        # Skip uma-m-1p1 - too large for RAM
        ('esen-sm-direct-all-omol', 'eSEN Small-MD (OMol25)'),
        ('esen-sm-conserving-all-omol', 'eSEN Small-Conserving (OMol25)'),
        ('esen-md-direct-all-omol', 'eSEN Medium-MD (OMol25)'),
        ('esen-sm-conserving-all-oc25', 'eSEN Small-Conserving (OC25)'),
        ('esen-md-direct-all-oc25', 'eSEN Medium-MD (OC25)'),
    ]
    
    print("=" * 105)
    print("Model Parameter Analysis")
    print("=" * 105)
    print()
    print(f"{'Model':<42} {'Description':<35} {'Total Params':<15} {'Active Params':<15}")
    print("-" * 105)
    
    for model_name, description in models:
        print(f"{model_name:<42} {description:<35} ", end="", flush=True)
        
        result = analyze_model(model_name)
        
        if result['success']:
            total_str = format_number(result['total'])
            active_str = format_number(result['active']) if result['active'] else "N/A"
            print(f"{total_str:<15} {active_str:<15}")
        else:
            print(f"{'ERROR':<15} {'':<15}")
            print(f"  → {result['error']}")
    
    print()
    print(f"{'uma-m-1p1':<42} {'UMA Medium v1.1':<35} {'1.40B*':<15} {'50M* (est.)':<15}")
    print("  * From documentation (too large to load)")
    print()
    print("=" * 105)
    print("Notes:")
    print("- For UMA models with MoLE architecture, 'Active Params' shows the subset used per prediction")
    print("- 'Total Params' includes all parameters in the model")
    print("- uma-m-1p1 values are from documentation as model is too large to load")
    print("=" * 105)
