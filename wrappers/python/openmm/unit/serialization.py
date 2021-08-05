"""
Module openmm.unit.serialization

Contains tools for converting Quantities to dicts that can be serialized by
various methods, included JSON and YAML, as well as several tools to
facilitate integration with the Python standard library package json.
"""
import json
import numpy as np
from .unit import Unit
from . import unit_definitions
from .quantity import Quantity, is_quantity

# These are implemented as separate functions so that they can be used by
# external developers in their own serialization routines. They can also
# dump to other formats (YAML, TOML, etc.)

def to_dict(quantity, as_bytes=True):
    """Convert a Quantity to a dict suitable for serialization.

    Parameters
    ----------
    quantity : Quantity
    """
    units = {p.name: int(power)
             for p, power in quantity.unit.iter_base_or_scaled_units()}

    if isinstance(quantity._value, np.ndarray):
        # treat numpy as a special case
        v = quantity._value
        value = {'shape': list(v.shape),
                 'dtype': str(v.dtype)}
        if as_bytes:
            value.update({
                'ndarray': 'bytes',
                'value': quantity._value.tobytes().decode('latin-1'),
            })
            ndarray = 'bytes'
        else:
            value.update({
                'ndarray': 'list',
                'value': quantity._value.tolist(),
            })
    else:
        value = quantity._value

    return {'__openmm.unit__': units,
            '__value__': value}

def from_dict(dct):
    """Convert the dict made by ``to_dict`` back into a Quantity.
    """
    units = Unit({})
    for u_name, u_power in dct['__openmm.unit__'].items():
        units *= getattr(unit_definitions, u_name) ** u_power

    if isinstance(dct['__value__'], dict):
        d = dct['__value__']
        # treat the numpy special case
        shape = d['shape']
        dtype = d['dtype']
        if d['ndarray'] == 'list':
            val = np.array(d['value'], dtype=dtype).reshape(shape)
        elif d['ndarray'] == 'bytes':
            as_bytes = d['value'].encode('latin-1')
            val = np.frombuffer(as_bytes, dtype=dtype).reshape(shape)
        else:
            raise RuntimeError("TODO")
    else:
        val = dct['__value__']

    return val * units

#####################
### JSON-SPECIFIC ###
#####################

# This sections contains helpers to interface the Python standard library
# json package.

def object_hook(dct):
    """Can be used for the object_hook parameter in json.loads
    """
    if '__openmm.unit__' in dct:
        return from_dict(dct)
    return dct

class UnitJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if is_quantity(obj):
            return to_dict(obj, as_bytes=False)
        return json.JSONEncoder.default(self, obj)

class UnitBytesJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Quantity):
            return to_dict(obj, as_bytes=True)
        return json.JSONEncoder.default(self, obj)

class UnitJSONDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super(UnitJSONDecoder, self).__init__(
            object_hook=object_hook, *args, **kwargs
        )
