#!/bin/sh

# Test Python
python -m simtk.testInstallation
cd python/tests && py.test -v
