from setuptools import setup, Extension 

MAX_KMER_SIZE = 31

cortex_include_dirs = ["alignment", "basic", "commands", "global", "graph", "graph_paths",
    "kmer", "paths", "tools"]
libdir = "../../libs"
cortex_ext = Extension("_cortex",
    ["_cortexmodule.c"],
    include_dirs=["../" + d for d in cortex_include_dirs] + [libdir, libdir + "/htslib/htslib"],
    define_macros=[("MAX_KMER_SIZE", MAX_KMER_SIZE)],
    extra_compile_args=["-std=c99"] 

)
setup(
    name = "cortex",
    version = "0.1",
    ext_modules = [cortex_ext],
    packages = ["cortex"], 
)
