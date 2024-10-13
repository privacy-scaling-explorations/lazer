.PHONY: all clean

all: _lazer_cffi.o

_lazer_cffi.o: FORCE
	python3 lazer_cffi_build.py ..

FORCE: ;

clean:
	rm -rf __pycache__
	rm -f *.so *.c *.o

