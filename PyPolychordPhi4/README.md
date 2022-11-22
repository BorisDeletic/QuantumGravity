requires: 
* [anesthetic](https://github.com/williamjameshandley/anesthetic) >=	2.0.0-beta.12
* [PolyChord](https://github.com/PolyChord/PolyChordLite) >= 1.20.1

Example script in `run_phi4.py`

Run as `mpirun -np [Your CPU core count] python run_phi4.py`


----

Added `run_polychord.cpp` snippet. From another project interfacing to a particle physics integrator library called STRIPPER

Can be compiled standalone (modulo the missing Integrands header from the STRIPPER lib which is private), the interfaces headers are linked against the PolyChordLite lib dir. Should give example of interfacing a likelihood call from an external c++ lib to polychord. Style is not so different to PyPolyChord in that it's a setup a likelihood call, a prior and hit run_pypolychord. 
