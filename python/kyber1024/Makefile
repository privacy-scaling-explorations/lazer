.PHONY: all clean

all: _kyber1024_params_cffi.o

_kyber1024_params_cffi.o: kyber1024_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py kyber1024_params.h ../..

kyber1024_params.h: kyber1024_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/kyber1024/kyber1024_params.py > ../python/kyber1024/kyber1024_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
