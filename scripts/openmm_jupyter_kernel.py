#!/usr/bin/env python3
import os
import sys

os.environ.setdefault("OPENMM_PLUGIN_DIR", "/usr/local/openmm/lib/plugins")
os.environ.setdefault("OPENMM_LIB_PATH", "/usr/local/openmm/lib")
lib = os.environ["OPENMM_LIB_PATH"]
os.environ["LD_LIBRARY_PATH"] = lib + (":" + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")

from ipykernel.kernelapp import IPKernelApp
IPKernelApp.launch_instance()
