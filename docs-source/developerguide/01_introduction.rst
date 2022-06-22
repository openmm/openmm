.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

Introduction
############

This guide describes the internal architecture of the OpenMM library.  It is
targeted at developers who want to add features to OpenMM, either by modifying
the core library directly or by writing plugins.  If you just want to write
applications that use OpenMM, you do not need to read this guide; the User's
Manual tells you everything you need to know.  This guide is intended for
people who want to contribute to OpenMM itself.

It is organized as follows:

* Chapter :numref:`the-core-library` describes the architecture of the core OpenMM library.  It
  discusses how the high level and low level APIs relate to each other, and the
  flow of execution between them.
* Chapter :numref:`writing-plugins` describes in detail how to write a plugin.  It focuses on the two
  most common types of plugins: those which define new Forces, and those which
  implement new Platforms.
* Chapter :numref:`the-reference-platform` discusses the architecture of the reference Platform, providing
  information relevant to writing reference implementations of new features.
* Chapter :numref:`the-cpu-platform` discusses the architecture of the CPU Platform, providing
  information relevant to writing CPU implementations of new features.
* Chapter :numref:`the-opencl-platform` discusses the architecture of the OpenCL Platform, providing
  information relevant to writing OpenCL implementations of new features.
* Chapter :numref:`the-cuda-platform` discusses the architecture of the CUDA Platform, providing
  information relevant to writing CUDA implementations of new features.
* Chapter :numref:`common-compute` describes the Common Compute framework, which lets you
  write a single implementation of a feature that can be used for both OpenCL and CUDA.


This guide assumes you are already familiar with the public API and how to use
OpenMM in applications.  If that is not the case, you should first read the
User's Manual and work through some of the example programs.  Pay especially
close attention to the “Introduction to the OpenMM Library” chapter, since it
introduces concepts that are important in understanding this guide.


