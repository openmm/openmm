#!/bin/bash

rm -rf xml/ html/

doxygen
cd xml
xsltproc combine.xslt index.xml >| ../all.xml

