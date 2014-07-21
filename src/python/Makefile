# simple makefile for development.

SRC=_cortexmodule.c
CFLAGS=-Wpointer-arith

ext3: ${SRC}
	CFLAGS=${CFLAGS} python3 setup.py build_ext --inplace

ext2: ${SRC}
	python setup.py build_ext --inplace

ctags:
	ctags *.c *.py test/*.py

