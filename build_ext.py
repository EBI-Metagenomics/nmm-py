import os
from os.path import join
from typing import List

from cffi import FFI

import imm

ffibuilder = FFI()
libs = ["nmm"]

ffibuilder.include(imm.ffibuilder)

folder = os.path.dirname(os.path.abspath(__file__))

with open(join(folder, "nmm", "nmm.h"), "r") as f:
    ffibuilder.cdef(f.read())

extra_link_args: List[str] = []
if "NMM_EXTRA_LINK_ARGS" in os.environ:
    extra_link_args += os.environ["NMM_EXTRA_LINK_ARGS"].split(os.pathsep)

ffibuilder.set_source(
    "nmm._ffi",
    r"""
    #include "nmm/nmm.h"
    """,
    libraries=libs,
    extra_link_args=extra_link_args,
    language="c",
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
