# python 3 params_cffi_build.py <params.h> <lazerdir>

import cffi
import sys
import os
import re

assert len(sys.argv) <= 3
file = sys.argv[1]
name = os.path.basename(file).split('.', 1)[0]

lazerdir = [os.path.abspath('..')]
if len(sys.argv) > 2:
    lazerdir = [os.path.abspath(sys.argv[2])]

with open(file) as f:
    data = f.read()

pointers = re.findall(r'\sstatic const lin_params_t (\S+)\s=', data)

cdefs = "const void *get_params (const wchar_t *name);\n"

csource = f"""
#include <wchar.h>
#include "lazer.h"
#include "{file}"

const void *get_params (const wchar_t *name) {{
"""
for ptr in pointers:
    csource += f"""
    if (wcscmp(name, L\"{ptr}\") == 0) {{
        return {ptr};
    }}
    """

csource += f"""
    return NULL;
}}
"""

ffibuilder = cffi.FFI()
ffibuilder.cdef(cdefs)
ffibuilder.set_source(f"_{name}_cffi", csource, include_dirs=lazerdir,
                      library_dirs=lazerdir, runtime_library_dirs=lazerdir)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
