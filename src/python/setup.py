from ez_setup import use_setuptools
use_setuptools()

import glob
from setuptools import setup, Extension 

MAX_KMER_SIZE = 31

cortex_include_dirs = ["alignment", "basic", "commands", "global", "graph", "graph_paths",
    "kmer", "paths", "tools"]
libdir = "../../libs"
cortex_objects = glob.glob("../../build/*/*")

extra_libs = [libdir + "/string_buffer/libstrbuf.a",
    libdir + "/seq-align/src/libalign.a"]

# htslib is compiled in as a shared library as I can't be 
# bothered recompiling it as PIC code. We therefore need 
# to set LD_LIBRARY_PATH at run time accordingly.
cortex_ext = Extension("_cortex",
    ["_cortexmodule.c", 
            # these have not been compiled, so it's easiest to include them here. 
            libdir + "/cJSON/cJSON.c", 
            libdir + "/misc/city.c",
            libdir + "/misc/mem_size.c"], 
    include_dirs=["../" + d for d in cortex_include_dirs] + [libdir, 
            libdir + "/htslib/htslib", libdir + "/seq_file"],
    define_macros=[("MAX_KMER_SIZE", MAX_KMER_SIZE), ("_USESAM", 1),
        ("CTXCHECKS", 1), ("MIN_KMER_SIZE", 3), ("_FORTIFY_SOURCE", 2)],
    extra_compile_args=["-std=c99", "-Wno-strict-prototypes"],
    extra_objects=cortex_objects + extra_libs,
    library_dirs=[libdir + "/htslib"],
    libraries=["hts"]
)
setup(
    name = "cortex",
    version = "0.1",
    ext_modules = [cortex_ext],
    packages = ["cortex"], 
)
