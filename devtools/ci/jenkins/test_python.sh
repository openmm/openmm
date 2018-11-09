#!/bin/sh

# Build & test Python
make PythonInstall
python -m simtk.testInstallation
cd python/tests && py.test -v
