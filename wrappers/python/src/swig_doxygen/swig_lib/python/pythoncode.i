%pythoncode %{

try:
    import numpy
except ImportError:
    numpy = None

import copy
import sys
import math
import functools
import operator
RMIN_PER_SIGMA=math.pow(2, 1/6.0)
RVDW_PER_SIGMA=math.pow(2, 1/6.0)/2.0
if sys.version_info[0] == 2:
    _string_types = (basestring,)
else:
    _string_types = (bytes, str)

import simtk.unit as unit
from simtk.openmm.vec3 import Vec3


%}

%pythonappend OpenMM::Context::Context %{
    self._system = args[0]
    self._integrator = args[1]
%}

%pythonprepend OpenMM::AmoebaAngleForce::addAngle %{
    try:
        length = args[3]
        if isinstance(args, tuple):
            args = list(args)
    except (NameError, UnboundLocalError):
        if unit.is_quantity(length):
            length = length.value_in_unit(unit.degree)
    else:
        if unit.is_quantity(length):
            args[3] = length.value_in_unit(unit.degree)
%}

%pythonprepend OpenMM::AmoebaAngleForce::setAngleParameters %{
    try:
        length = args[4]
        if isinstance(args, tuple):
            args = list(args)
    except (NameError, UnboundLocalError):
        if unit.is_quantity(length):
            length = length.value_in_unit(unit.degree)
    else:
        if unit.is_quantity(length):
            args[4] = length.value_in_unit(unit.degree)
%}

%pythonprepend OpenMM::AmoebaTorsionTorsionForce::setTorsionTorsionGrid %{
    def deunitize_grid(grid):
        if isinstance(grid, tuple):
            grid = list(grid)
        for i, row in enumerate(grid):
            if isinstance(row, tuple):
                row = list(row)
                grid[i] = row
            for i, column in enumerate(row):
                if isinstance(column, tuple):
                    column = list(column)
                    row[i] = column
                # Data is angle, angle, energy, de/dang1, de/dang2, d^2e/dang1dang2
                if unit.is_quantity(column[0]):
                    column[0] = column[0].value_in_unit(unit.degree)
                if unit.is_quantity(column[1]):
                    column[1] = column[1].value_in_unit(unit.degree)
                if unit.is_quantity(column[2]):
                    column[2] = column[2].value_in_unit(unit.kilojoule_per_mole)
                if len(column) > 3 and unit.is_quantity(column[3]):
                    column[3] = column[3].value_in_unit(unit.kilojoule_per_mole/unit.radians)
                if len(column) > 4 and unit.is_quantity(column[4]):
                    column[4] = column[4].value_in_unit(unit.kilojoule_per_mole/unit.radians)
                if len(column) > 5 and unit.is_quantity(column[5]):
                    column[5] = column[5].value_in_unit(unit.kilojoule_per_mole/unit.radians**2)
        return grid
    try:
        grid = copy.deepcopy(args[1])
        if isinstance(args, tuple):
            args = list(args)
    except (NameError, UnboundLocalError):
        try:
            # Support numpy arrays
            grid = grid.tolist()
        except AttributeError:
            grid = copy.deepcopy(grid)
        grid = deunitize_grid(grid)
    else:
        args[1] = deunitize_grid(grid)
%}
