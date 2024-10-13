.PHONY: all clean

all: _swoosh_params_cffi.o

_swoosh_params_cffi.o: swoosh_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py swoosh_params.h ../..

swoosh_params.h: swoosh_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/swoosh/swoosh_params.py > ../python/swoosh/swoosh_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
