from Cython.Build import build_ext
from setuptools import Extension, setup

setup(
    name="edlib",
    description="Lightweight, super fast library for sequence alignment using edit (Levenshtein) distance.",
    long_description="",
    version="1.3.9",
    url="https://github.com/Martinsos/edlib",
    author="Martin Sosic",
    author_email="sosic.martin@gmail.com",
    license="MIT",
    keywords="edit distance levenshtein align sequence bioinformatics",
    # Build instructions
    ext_modules=[
        Extension(
            "utils",
            ["utils.pyx", "edlib/src/edlib.cpp"],
            include_dirs=["edlib/include"],
            depends=["edlib/include/edlib.h"],
            language="c++",
            extra_compile_args=["-O3", "-std=c++11"],
        )
    ],
    cmdclass={"build_ext": build_ext},
)
