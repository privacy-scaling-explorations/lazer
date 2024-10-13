.PHONY: all clean

all: _anon_cred_params_cffi.o

_anon_cred_params_cffi.o: anon_cred_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py anon_cred_params.h ../..

anon_cred_params.h: anon_cred_p1_params.py anon_cred_p2_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/anon_cred/anon_cred_p1_params.py > ../python/anon_cred/anon_cred_params.h
	cd ../../scripts && sage lin-codegen.sage ../python/anon_cred/anon_cred_p2_params.py >> ../python/anon_cred/anon_cred_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
