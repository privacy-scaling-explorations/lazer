.PHONY: all clean

all: _blindsig_params_cffi.o

_blindsig_params_cffi.o: blindsig_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py blindsig_params.h ../..

blindsig_params.h: blindsig_p1_params.py blindsig_p2_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/blindsig/blindsig_p1_params.py > ../python/blindsig/blindsig_params.h
	cd ../../scripts && sage lin-codegen.sage ../python/blindsig/blindsig_p2_params.py >> ../python/blindsig/blindsig_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
