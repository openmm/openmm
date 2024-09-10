from __future__ import print_function, division, absolute_import

__author__ = "Raul P. Pelaez"
__version__ = "1.0"

from openmm.app.internal.xtc_utils import (
    xtc_rewrite_with_new_timestep,
    xtc_write_frame,
    get_xtc_nframes,
    get_xtc_natoms,
)
import os
from openmm import Vec3
from openmm.unit import nanometers, picoseconds, is_quantity, norm
import tempfile
import shutil


class XTCFile(object):

    """XTCFile provides methods for creating XTC files.
    To use this class, create a XTCFile object, then call writeModel() once for each model in the file.
    """

    def __init__(self, fileName, topology, dt, firstStep=0, interval=1, append=False):
        """Create a XTC file, or open an existing file to append.

        Parameters
        ----------
        fileName : str
            A file name to write to
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
        if not isinstance(fileName, str):
            raise TypeError("fileName must be a string")
        self._filename = fileName
        self._topology = topology
        self._firstStep = firstStep
        self._interval = interval
        self._modelCount = 0
        if is_quantity(dt):
            dt = dt.value_in_unit(picoseconds)
        self._dt = dt
        if append:
            if not os.path.isfile(self._filename):
                raise FileNotFoundError(f"The file '{self._filename}' does not exist.")
            self._modelCount = get_xtc_nframes(self._filename.encode("utf-8"))
            natoms = get_xtc_natoms(self._filename.encode("utf-8"))
            if natoms != len(list(topology.atoms())):
                raise ValueError(
                    f"The number of atoms in the topology ({len(list(topology.atoms()))}) does not match the number of atoms in the XTC file ({natoms})"
                )
        else:
            if os.path.isfile(self._filename) and os.path.getsize(self._filename) > 0:
                raise FileExistsError(f"The file '{self._filename}' already exists.")

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
        if self._topology.getNumAtoms() != len(positions):
            raise ValueError("The number of positions must match the number of atoms")
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometers)
        import numpy as np
        if np.isnan(positions).any():
            raise ValueError(
                "Particle position is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan"
            )
        if np.isinf(positions).any():
            raise ValueError(
                "Particle position is infinite.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan"
            )
        if (np.abs(positions * 1000) > (2 ** 31 - 1)).any():
            raise ValueError("Particle position is too large for XTC format")

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
            with tempfile.TemporaryDirectory() as temp:
                fname = os.path.join(temp, "temp.xtc")
                xtc_rewrite_with_new_timestep(
                    self._filename.encode("utf-8"),
                    fname.encode("utf-8"),
                    self._firstStep,
                    self._interval,
                    self._dt,
                )
                shutil.copyfile(fname, self._filename)
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
        step = (self._modelCount - 1) * self._interval + self._firstStep
        time = step * self._dt
        xtc_write_frame(
            self._filename.encode("utf-8"),
            np.array(positions, dtype=np.float32),
            boxVectors,
            np.float32(time),
            np.int32(step),
        )
