lazer
=====
A library for lattice-based zero-knowledge proofs.


Dependencies
------------
The following hardware and software is required to build
and run lazer:

- Linux amd64 / x86-64 system
- kernel version >= 4.18
- avx512 and aes instruction set extensions
- gcc compiler >= 13.2
- make >= 4.2, cmake >= 3.26
- sagemath >= 10.2
- python3 >= 3.10, its development package and the cffi package 

The following software is required to build the documentation:

 - sphinx >= 5.3 and the sphinxcontrib-bibtex package

The package was tested on Ubuntu 20.04 LTS and RHEL 8.10
with an 11th Gen Intel Core i7-11850H @ 2.50GHz and
with the dependencies' versions equal to the listed above.
We assume that it works with greater versions as well,
as long as they are backward compatible. It *may*
work with lesser versions as well.
Note that during compilation, the cpu_features package
will be cloned from GitHub, so also git and an internet
connection are required.


Building the library
--------------------
To build the lazer C library, from the base directory, run:

`make all`

Optional: use make's `-j` option to speed up compilation.

To build lazer's python module, change to the `python` subdirectory and run:

`make`

(If this step fails, check that the python development package is installed.)


Building and viewing the documentation
--------------------------------------
To build the html documentation, change to the `docs` subdirectory and run:

`make html`

(If you have multiple python versions installed, make sure sphinx uses the same version that was used to build the python module.)

To view the documentation open `build/html/index.html` in a browser e.g.:

`firefox build/html/index.html`


Building the C demos
--------------------
The `demos` subdirectory has the C demos:

- a blind signature implementation in `blindsig` 
- a PoK of Kyber1024 secrets (paper section 3.2,3.3) in `kyber1024`

To build a demo, go to the corresponding subdirectory and run:

`make`

This creates an executable called `<name>-demo`.
Run the demo via:

`./<name>-demo`


Building and running the python demos for linear relations with norms
---------------------------------------------------------------------
In the `python` subdirectory are all the python demos mentioned in the paper:

- an anonymous credentials (paper section 6.2) in `anon_cred` 
- a blind signature implementation in `blindsig` 
- a PoK of Kyber1024 secrets (paper section 6.1) in `kyber1024`
- the proof required in the Swoosh NIKE in `swoosh` 
- a proof for the general lattice relation As=t in `demo`

To build a demo, go to the corresponding subdirectory and run:

`make`

Each demo is implemented in a python script with the same name as its directory i.e., `<demo>.py`. To run a demo,
go to the corresponding subdirectory `<demo>` and run:

`python3 <demo>.py`


Running the python demo for LaBRADOR
------------------------------------
An aggregate signature (paper section 6.3) implementation is in the `python` subdirectory.

Run it via:

`python3 agg_sig.py`


Generating the demos' proof parameters from a specification
-----------------------------------------------------------

A demo's proof parameters are specified in files named `*params.py` in the demo's subdirectory. For convenience, the code gererated from those specifications is included in the package (the `*params.h` header files) such that make only runs the code generator when the specification is changed (or when the header file was deleted).
Since the code generator is a sagemath script (`scripts/lin-codegen.sage`) that calls the lattice-estimator multiple times, the code generation process may take multiple minutes, especially for large parameter set, like for Swoosh.

The code generator can be used from the `scripts` subdirectory via:

`sage lin-codegen.sage <specification> > <headerfile>`


Instructions for artifact evaluation and result reproduction
------------------------------------------------------------

1. Set up a system that meets lazer's requirements decribed above.
2. Obtain and build the library and demos as described above.
3. The code for kyber proof in C from sections 3.2, 3.3 was extended by a main function and sample inpts and is available as a runnable demo in `demos/kyber1024/kyber1024-demo`.
4. The code for kyber proof in python from section 6.1 is available as a runnable demo in `python/kyber1024/kyber1024.py`.
5. The code for the anonymous credentials from section 6.2 is avalable as a runnable demo in `python/anon_cred/anon_cred.py`.
6. The code for the aggregate signature from section 6.3 is available as a runnable demo in `python/agg_sig.py`.
7. Detailed descriptions of the demos are in the documentation.
8. If you want to do anything beyond verifying the paper's results, check out the documentation of the python module's interface.


