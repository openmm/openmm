#!/bin/bash

rm -f *.xml
rm -f OpenMM_headers*
rm -f *.cxx

SWIG_REV_FILE_OPENMM=RevisionNumber_OpenMM.txt
SWIG_REV_FILE_PYOPENMM=RevisionNumber_pyopenmm.txt

PYTHON_PACKAGE_DIR=../../simtk/chem/openmm


if [ -n "$OPENMM_SVN_PATH" ] ; then
  SVN_INFO=$(svn info $OPENMM_SVN_PATH 2>/dev/null | awk '$0~/^Revision:/{print $2}')
  if [ -n "$SVN_INFO" ] ; then
    echo  $SVN_INFO >| $SWIG_REV_FILE_OPENMM
  fi
fi

if [ -z "$OPENMM_INCLUDE_PATH" ] ; then
  OPENMM_INCLUDE_PATH=../../OpenMM/include
fi

SVN_INFO=$(svn info 2>/dev/null | awk '$0~/^Revision:/{print $2}')
if [ -n "$SVN_INFO" ] ; then
  echo  $SVN_INFO >| $SWIG_REV_FILE_PYOPENMM
fi

cd doxygen
rm -rf xml/ html/
echo "Calling doxygen >| doxygen.out 2>| doxygen.err"
doxygen >| doxygen.out 2>| doxygen.err
if [ "$?" -ne "0" ]; then
  echo "ERROR: doxygen did not run"
  echo "       See doxygen/doxygen.out and doxygen/doxygen.err!"
  echo "Exiting!"
  exit 1
fi
cd xml
echo "Calling xsltproc combine.xslt index.xml >| ../../OpenMM_headers.xml"
xsltproc combine.xslt index.xml >| ../../OpenMM_headers.xml
cd ../../


#Build ref platform only
echo "Calling swigInputBuilder.py -i OpenMM_headers.xml -c swigInputConfig.py -o OpenMM_headers.i -d OpenMM_docstring.i -a swig_lib/python/pythonprepend.i  -z swig_lib/python/pythonappend.i "
python swigInputBuilder.py -i OpenMM_headers.xml \
                           -c swigInputConfig.py \
                           -o OpenMM_headers.i \
                           -d OpenMM_docstring.i \
                           -a swig_lib/python/pythonprepend.i \
                           -z swig_lib/python/pythonappend.i \
                        >| swigInputBuilder.out 

USING_SWIG_VERSION=$(swig -version | awk '/^SWIG Version/{print $3}')
USING_SWIG_VERSION_NUM=$(echo $USING_SWIG_VERSION | awk -F. '{print 1000*(1000*$1+$2)+$3}')
if (( $USING_SWIG_VERSION_NUM < 2000000 )) ; then
   echo "ERROR: Using swig $USING_SWIG_VERSION.  Must use 2.0.0 or better"
   exit 1
fi

echo "calling swig -python -c++ -Wall -outdir $PYTHON_PACKAGE_DIR -o OpenMMSwig.cxx OpenMM.i"
swig -python -c++ -Wall \
     -outdir $PYTHON_PACKAGE_DIR \
     -o OpenMMSwig.cxx \
     OpenMM.i


echo "Done: swig -python -c++\n"

