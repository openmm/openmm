from __future__ import print_function, division, absolute_import

__author__ = "Raul P. Pelaez"
__version__ = "1.0"

from openmm.app.internal.xtc_utils import (
    read_xtc,
    read_xtc_frames,
    xtc_write_frame,
    get_xtc_nframes,
    get_xtc_natoms,
)
import numpy as np
import os
from openmm import Vec3
from openmm.unit import nanometers, femtoseconds, is_quantity, norm
import math
import tempfile
import shutil


class XTCFile(object):

    """XTCFile provides methods for creating XTC files.
    To use this class, create a XTCFile object, then call writeModel() once for each model in the file.
    """

    def __init__(self, file, topology, dt, firstStep=0, interval=1, append=False):
        """Create a XTC file, or open an existing file to append.

        Parameters
        ----------
        file : file
            A file to write to
        topology : Topology
            The Topology defining the molecular system being written
        dt : time
            The time step used in the trajectory
        firstStep : int=0
            The index of the first step in the trajectory
        interval : int=1
            The frequency (measured in time steps) at which states are written
            to the trajectory
        append : bool=False
            If True, open an existing XTC file to append to.  If False, create a new file.
        """
        self._filename = file.name
        file.close()
        self._topology = topology
        self._firstStep = firstStep
        self._interval = interval
        self._modelCount = 0
        if is_quantity(dt):
            dt = dt.value_in_unit(femtoseconds)
        self._dt = dt
        if append:
            self._modelCount = get_xtc_nframes(self._filename.encode("utf-8"))
            natoms = get_xtc_natoms(self._filename.encode("utf-8"))
            if natoms != len(list(topology.atoms())):
                raise ValueError(
                    f"The number of atoms in the topology ({len(list(topology.atoms()))}) does not match the number of atoms in the XTC file ({natoms})"
                )

    def _getNumFrames(self):
        return get_xtc_nframes(self._filename.encode("utf-8"))

    def writeModel(self, positions, unitCellDimensions=None, periodicBoxVectors=None):
        """Write out a model to the XTC file.

        The periodic box can be specified either by the unit cell dimensions
        (for a rectangular box), or the full set of box vectors (for an
        arbitrary triclinic box).  If neither is specified, the box vectors
        specified in the Topology will be used. Regardless of the value
        specified, no dimensions will be written if the Topology does not
        represent a periodic system.

        Parameters
        ----------
        positions : list
            The list of atomic positions to write
        unitCellDimensions : Vec3=None
            The dimensions of the crystallographic unit cell.
        periodicBoxVectors : tuple of Vec3=None
            The vectors defining the periodic box.
        """
        if len(list(self._topology.atoms())) != len(positions):
            raise ValueError("The number of positions must match the number of atoms")
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometers)
        if any(math.isnan(norm(pos)) for pos in positions):
            raise ValueError(
                "Particle position is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan"
            )
        if any(math.isinf(norm(pos)) for pos in positions):
            raise ValueError(
                "Particle position is infinite.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan"
            )

        self._modelCount += 1
        if (
            self._interval > 1
            and self._firstStep + self._modelCount * self._interval > 1 << 31
        ):
            # This will exceed the range of a 32 bit integer.  To avoid crashing or producing a corrupt file,
            # update the file to say the trajectory consisted of a smaller number of larger steps (so the
            # total trajectory length remains correct).
            self._firstStep //= self._interval
            self._dt *= self._interval
            self._interval = 1
            with tempfile.NamedTemporaryFile() as temp:
                nframes = self._getNumFrames()
                for i in range(nframes):
                    read_frame = np.array([i]).astype(np.int32)
                    _coords, _box, _time, _step = read_xtc_frames(
                        self._filename.encode("utf-8"), read_frame
                    )
                    _step = self._firstStep + i * self._interval
                    _time = _step * self._dt
                    xtc_write_frame(
                        temp.name.encode("utf-8"),
                        _coords[:, :, 0],
                        _box[:, :, 0],
                        _time,
                        _step,
                    )
                shutil.copyfile(temp.name, self._filename)
        boxVectors = self._topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            if periodicBoxVectors is not None:
                boxVectors = periodicBoxVectors
                if is_quantity(boxVectors):
                    boxVectors = boxVectors.value_in_unit(nanometers)
            elif unitCellDimensions is not None:
                if is_quantity(unitCellDimensions):
                    unitCellDimensions = unitCellDimensions.value_in_unit(nanometers)
                boxVectors = (
                    Vec3(unitCellDimensions[0], 0, 0),
                    Vec3(0, unitCellDimensions[1], 0),
                    Vec3(0, 0, unitCellDimensions[2]),
                )
            boxVectors = np.array(
                [[vec.x, vec.y, vec.z] for vec in boxVectors], dtype=np.float32
            )
        else:
            boxVectors = np.zeros((3, 3)).astype(np.float32)
        step = self._modelCount * self._interval + self._firstStep
        time = step * self._dt
        xtc_write_frame(
            self._filename.encode("utf-8"),
            np.array(positions).astype(np.float32),
            boxVectors,
            np.float32(time),
            np.int32(step),
        )
