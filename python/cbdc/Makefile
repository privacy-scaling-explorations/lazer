.PHONY: all clean

all: _cbdc_params_cffi.o

_cbdc_params_cffi.o: cbdc_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py cbdc_params.h ../..

cbdc_params.h: cbdc_p1_params.py cbdc_p2_params.py cbdc_popen_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/cbdc/cbdc_p1_params.py > ../python/cbdc/cbdc_params.h
	cd ../../scripts && sage lin-codegen.sage ../python/cbdc/cbdc_p2_params.py >> ../python/cbdc/cbdc_params.h
	cd ../../scripts && sage lin-codegen.sage ../python/cbdc/cbdc_popen_params.py >> ../python/cbdc/cbdc_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
