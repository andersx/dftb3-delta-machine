all: intel-openmp

gnu:
	f2py -c -m fkernels fkernels.f90

gnu-openmp:
	f2py -c -m fkernels fkernels.f90 --f90flags='-fopenmp' -lgomp

intel:
	f2py -c -m fkernels fkernels.f90 --fcompiler=intelem

intel-openmp:
	f2py -c -m fkernels fkernels.f90 --fcompiler=intelem --f90flags='-qopenmp' -liomp5

clean:
	rm fkernels.so *.pyc *.png
