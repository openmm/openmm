# Install Python examples.

INSTALL(FILES simulateAmber.py input.inpcrd input.prmtop DESTINATION examples/python-examples)
INSTALL(FILES simulateCharmm.py ala_ala_ala.pdb ala_ala_ala.psf charmm22.par charmm22.rtf DESTINATION examples/python-examples)
INSTALL(FILES simulateGromacs.py input.gro input.top DESTINATION examples/python-examples)
INSTALL(FILES simulatePdb.py input.pdb DESTINATION examples/python-examples)
INSTALL(FILES argon-chemical-potential.py DESTINATION examples/python-examples)
