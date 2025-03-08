OpenMM extra utility scripts
============================

This directory contains standalone utility scripts for use with OpenMM.

* `optimizepme.py`: Optimizes parameters for PME simulations by running a series
  of simulations with different PME parameters and choosing the combination of
  parameters giving the best performance.
* `rigid.py`: Implements rigid bodies by adding constraints between some
  particles and converting others to virtual sites.
