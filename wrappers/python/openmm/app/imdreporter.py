"""
imdreporter.py: Sends data about a simulation to another program

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import openmm as mm
import openmm.unit as unit
import socket
import struct
import threading

class IMDReporter(object):
    """IMDReporter uses the IMDv3 protocol (https://imdclient.readthedocs.io/en/latest/protocol_v3.html)
    to send information about a running simulation to another process.

    To use it, create an IMDReporter, then add it to the Simulation's list of reporters.

    When you create an IMDReporter, it opens a socket and begins listening for connections.  Once a
    client connects, it begins sending information at regular intervals, such as energy or particle
    positions.  Use constructor arguments to select what information to send.

    The reporter can operate in either blocking or non-blocking mode.  In blocking mode, the simulation
    will only run while a client is connected.  If the reporter is asked to generate a report and no
    client is currently connected, it will block until one connects.  In non-blocking mode, reports are
    simply skipped if no client is connected and the simulation continues running.

    Once a client connects, it can ask the reporter to change between blocking and non-blocking modes.
    This is useful, for example, if you want an external program to monitor the initial part of a
    simulation, but then disconnect and allow the simulation to continue running.
    """

    def __init__(self, reportInterval: int, port: int = 8888, positions: bool = False,
                 velocities: bool = False, forces: bool = False, energy: bool = False,
                 enforcePeriodicBox: bool | None = None, blocking: bool = True):
        """Create an IMDReporter.

        Parameters
        ----------
        reportInterval: int
            The interval (in time steps) at which to send information
        port: int
            The port on which to listen for connections
        positions: bool
            Whether to send particle positions to the client
        velocities: bool
            Whether to send particle velocities to the client
        forces: bool
            Whether to send particle forces to the client
        energy: bool
            Whether to send energy and temperature information to the client
        enforcePeriodicBox: optional bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        blocking: bool
            Whether the reporter should begin in blocking or non-blocking mode
        """
        self._reportInterval = reportInterval
        self._defaultReportInterval = reportInterval
        self._port = port
        self._positions = positions
        self._velocities = velocities
        self._forces = forces
        self._energy = energy
        self._enforcePeriodicBox = enforcePeriodicBox
        self._blocking = blocking
        self._paused = False
        self._connected = False
        self._closed = False
        self._disconnect = False
        self._initialized = False
        self._reports = []
        self._includes = []
        if positions:
            self._includes.append('positions')
        if velocities:
            self._includes.append('velocities')
        if forces:
            self._includes.append('forces')
        if energy:
            self._includes.append('energy')
        self._receivedPacket = bytearray(8)
        self._receivedView = memoryview(self._receivedPacket)
        self._receivedBytes = 0
        self._lock = threading.RLock()
        self._condition = threading.Condition(self._lock)
        self._thread = threading.Thread(target = lambda: self._runServerThread())
        self._thread.start()

    def close(self):
        """Close the socket and release all resources associated with this reporter."""
        with self._condition:
            self._closed = True
            self._condition.notify()
        self._thread.join()

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        dict
            A dictionary describing the required information for the next report
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return {'steps':steps, 'periodic':self._enforcePeriodicBox, 'include':self._includes}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._initialized:
            self._initializeConstants(simulation)
            self._initialized = True
        with self._condition:
            while not self._connected or self._paused:
                if not self._connected and not self._blocking:
                    # No client is connected and we're in non-blocking mode, so just return.
                    self._reports = []
                    return
                self._condition.wait()
            # Build a dict containing the information to report.
            integrator = simulation.context.getIntegrator()
            report = {'state': state, 'dt': integrator.getStepSize().value_in_unit(unit.picosecond)}
            if self._energy:
                if hasattr(integrator, 'computeSystemTemperature'):
                    temperature = integrator.computeSystemTemperature()
                else:
                    temperature = (2*state.getKineticEnergy()/(self._dof*unit.MOLAR_GAS_CONSTANT_R))
                report['temperature'] = temperature.value_in_unit(unit.kelvin)
            # Add it to the queue and wake up the other thread.
            self._reports.append(report)
            self._condition.notify()

    def _initializeConstants(self, simulation):
        """Initialize a set of constants required for the reports

        Parameters
        - simulation (Simulation) The simulation to generate a report for
        """
        system = simulation.system
        if self._energy:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0*unit.dalton:
                    dof += 3
            for i in range(system.getNumConstraints()):
                p1, p2, distance = system.getConstraintParameters(i)
                if system.getParticleMass(p1) > 0*unit.dalton or system.getParticleMass(p2) > 0*unit.dalton:
                    dof -= 1
            if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                dof -= 3
            self._dof = dof

    def _runServerThread(self):
        """Listen for connections from clients and communicate with them.  This method runs on its own thread."""
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(('', self._port))
            s.listen()
            while not self._closed:
                conn, addr = s.accept()
                with conn:
                    # A client has connected.  Send the initial information.
                    self._sendHandshake(conn)
                    self._sendSessionInfo(conn)
                    conn.setblocking(False)
                    self._disconnect = False
                    while not self._disconnect and not self._closed:
                        # The main loop.  Alternate between sending reports and listening for messages from the client.
                        try:
                            with self._lock:
                                while len(self._reports) > 0:
                                    self._sendReport(self._reports[0], conn)
                                    del self._reports[0]
                            while self._processPacket(conn):
                                pass
                        except (ConnectionError, BlockingIOError):
                            # This usually indicates the client disconnected unexpectedly.  Treat it as if we had
                            # received a disconnect message.
                            self._disconnect = True
                        with self._condition:
                            # Wait for something to do.  If a report is generated, this will be interrupted and
                            # return immediately.  If a message comes from the client we don't get any notification,
                            # so use a short timeout.
                            self._condition.wait(0.5)
                with self._lock:
                    self._connected = False

    def _sendHandshake(self, socket):
        """Send the handshake packet to establish a connection."""
        socket.sendall(struct.pack('!i', 4)) # Handshake
        socket.sendall(struct.pack('i', 3)) # IMD version in native byte order

    def _sendSessionInfo(self, socket):
        """Send the packet describing what information will be included in each report."""
        socket.sendall(struct.pack('!ii', 10, 7)) # Session info
        socket.sendall(struct.pack('bbbbbbb',
                                   1, # Always send time
                                   1 if self._energy else 0,
                                   1 if self._positions else 0,
                                   1 if self._positions else 0,
                                   0, # Always say coordinates are unwrapped, since even wrapped ones can extend outside the box
                                   1 if self._velocities else 0,
                                   1 if self._forces else 0))

    def _sendReport(self, report, socket):
        """Send a report to the client."""
        state = report['state']
        socket.sendall(struct.pack('!ii', 12, 1)) # Time
        socket.sendall(struct.pack('ddq',
                                   report['dt'],
                                   state.getTime().value_in_unit(unit.picoseconds),
                                   state.getStepCount()))
        if self._energy:
            potential = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
            socket.sendall(struct.pack('!ii', 1, 1)) # Energy
            socket.sendall(struct.pack('ifffffffff',
                                       min(state.getStepCount(), 2**31-1),
                                       report['temperature'],
                                       potential+kinetic,
                                       potential,
                                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        import numpy as np
        if self._positions:
            vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom).astype(np.float32)
            socket.sendall(struct.pack('!ii', 13, 1)) # Box
            socket.sendall(vectors.data)
            positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom).astype(np.float32)
            socket.sendall(struct.pack('!ii', 2, positions.shape[0])) # Coordinates
            socket.sendall(positions.data)
        if self._velocities:
            velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.angstrom/unit.picosecond).astype(np.float32)
            socket.sendall(struct.pack('!ii', 14, velocities.shape[0])) # Velocities
            socket.sendall(velocities.data)
        if self._forces:
            forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.angstrom).astype(np.float32)
            socket.sendall(struct.pack('!ii', 15, forces.shape[0])) # Forces
            socket.sendall(forces.data)

    def _processPacket(self, socket):
        """Receive a packet from the client and process it."""
        while self._receivedBytes < len(self._receivedView):
            try:
                self._receivedBytes += socket.recv_into(self._receivedView[self._receivedBytes:])
            except BlockingIOError:
                # There's no more data to read.
                return False
        self._receivedBytes = 0
        packetType, data = struct.unpack('!ii', self._receivedPacket)
        with self._condition:
            match packetType:
                case 0: # Disconnect
                    self._disconnect = True
                    self._condition.notify()
                case 3: # Go
                    self._connected = True
                    self._condition.notify()
                case 5: # Kill
                    # We choose to ignore this request, since there's no clean way for the report to stop
                    # the simulation.  Instead just disconnect.
                    self._disconnect = True
                    self._condition.notify()
                case 6: # MD Communication
                    raise ValueError('IMDReporter: The "MD Communication" message is not supported')
                case 7: # Pause
                    self._paused = True
                case 8: # Transmission rate
                    if data < 1:
                        self._reportInterval = self._defaultReportInterval
                    else:
                        self._reportInterval = data
                case 11: # Resume
                    self._paused = False
                    self._condition.notify()
                case 16: # Wait
                    self._blocking = (data != 0)
                    self._condition.notify()
        return True