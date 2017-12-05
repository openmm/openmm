.. include:: header.rst

Introduction
############

OpenMM consists of two parts:

#. A set of libraries that lets programmers easily add molecular simulation
   features to their programs
#. An “application layer” that exposes those features to end users who just want
   to run simulations


This guide is divided into three sections:

* Part I (Chapters :ref:`the-openmm-application-layer-introduction`\ -\ :ref:`creating-force-fields`\ )
  describes the application layer.  It is relevant to all users, but especially relevant to people
  who want to use OpenMM as a stand-alone application for running simulations.
* Part II (Chapters :ref:`the-openmm-library-introduction`\ -\ :ref:`drude-plugin`\ )
  describes how to use the OpenMM libraries within your own applications.  It is primarily
  relevant to programmers who want to write simulation applications.
* Part III (Chapters :ref:`the-theory-behind-openmm-introduction`\ -\ :ref:`other-features`\ )
  describes the mathematical theory behind the features found in OpenMM.  It is relevant to all users.


Online Resources
****************

You can find more documentation and other material at our website
http://openmm.org.   Among other things there is a discussion forum,
a wiki, and videos of lectures on using OpenMM.


Referencing OpenMM
******************

Any work that uses OpenMM should cite the following publication:

P. Eastman, J. Swails, J. D. Chodera, R. T. McGibbon, Y. Zhao, K. A. Beauchamp,
L.-P. Wang, A. C. Simmonett, M. P. Harrigan, C. D. Stern, R. P. Wiewiora,
B. R. Brooks, and V. S. Pande. "OpenMM 7: Rapid development of high performance
algorithms for molecular dynamics." PLOS Comp. Biol. 13(7): e1005659. (2017)

We depend on academic research grants to fund the OpenMM development efforts;
citations of our publication will help demonstrate the value of OpenMM.


Acknowledgments
***************

OpenMM software and all related activities, such as this manual, are funded by
the Simbios National Center for Biomedical Computing through the National
Institutes of Health Roadmap for Medical Research, Grant U54 GM072970, and by
National Institutes of Health grant R01-GM062868.

