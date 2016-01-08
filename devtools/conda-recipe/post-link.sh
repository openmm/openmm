if [[ "$OSTYPE" == "darwin"* ]]; then
    # http://answers.opencv.org/question/4134/
    # make the lib/libOpenMM*.dylib and lib/plugins/libOpenMM*.dylib
    # point to each other using absolute paths
    find "${PREFIX}/lib" -type f -name "libOpenMM*.dylib" -print0 | while IFS="" read -r -d "" dylibpath; do
        otool -L $dylibpath | grep libOpenMM | tr -d ':' | while read -a libs; do
            if [ "${file}" != "${libs[0]}" ]; then
                install_name_tool -change ${libs[0]} "${PREFIX}/lib/"`basename ${libs[0]}` $dylibpath
            fi
       done
    done

    # make site-packages/simtk/openmm/_openmm.so point to the correct
    # libraries in lib/
    find $PREFIX/lib/python*/site-packages/simtk/openmm -type f -name "_openmm.so" -print0 | while IFS="" read -r -d "" dylibpath; do
        otool -L $dylibpath | grep libOpenMM | tr -d ':' | while read -a libs; do
            if [ "${file}" != "${libs[0]}" ]; then
                install_name_tool -change ${libs[0]} "${PREFIX}/lib/"`basename ${libs[0]}` $dylibpath
            fi
       done
    done
fi

cat <<EOF > _temp.py
import os.path
import simtk
dir = os.path.dirname(simtk.__file__)
fn = os.path.join(dir, 'openmm', 'version.py')
f = open(fn, 'a')
f.write("\nopenmm_library_path = '%s/lib'\n" % os.environ['PREFIX'])
EOF
$PREFIX/bin/python _temp.py
rm _temp.py
