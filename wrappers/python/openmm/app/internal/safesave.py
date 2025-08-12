"""
safesave.py: Helper module to ensure atomic overwrite/backup of existing files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Evan Pretti
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import itertools
import os

def _getTempFilename(prefix):
    """
    Returns the name of a temporary file starting with a given prefix that is
    guaranteed not to exist already.  Upon successful return of this function,
    an empty file with the name returned will have been created.

    Parameters
    ----------
    prefix : str
        The prefix of the temporary file name to create.  If a path with
        multiple components, the directory in which to create the file must
        exist.

    Returns
    -------
    str
        The temporary file name created.
    """

    for index in itertools.count():
        name = f'{prefix}.{index}.tmp'
        try:
            with open(name, 'x'):
                return name
        except FileExistsError:
            pass

def save(data, filename):
    """
    Saves data to a specified file.  If the file exists, it will be overwritten
    atomically, or if this is not possible, a backup copy of the existing data
    will be created during overwriting and deleted once it is successful.

    Parameters
    ----------
    data : bytes or str
        The data to write.  If bytes, the file will be opened in binary mode; if
        str, in text mode.
    filename : str
        The filename to write to.
    """

    if isinstance(data, bytes):
        mode = 'wb'
    elif isinstance(data, str):
        mode = 'w'
    else:
        raise ValueError('Expected bytes or str')

    tempFilename1 = _getTempFilename(filename)
    with open(tempFilename1, mode) as file:
        file.write(data)

    try:
        # If the target file already exists, rename() should overwrite
        # atomically on POSIX and fail with a FileExistsError on Windows.
        os.rename(tempFilename1, filename)
    except FileExistsError:
        # Make a backup copy since replace() on Windows may not be atomic.
        tempFilename2 = _getTempFilename(filename)
        os.replace(filename, tempFilename2)
        os.replace(tempFilename1, filename)
        os.remove(tempFilename2)
