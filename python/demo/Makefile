.PHONY: all clean

all: _demo_params_cffi.o

_demo_params_cffi.o: demo_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py demo_params.h ../..

demo_params.h: demo_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/demo/demo_params.py > ../python/demo/demo_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
