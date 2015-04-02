# PYQE Interfaces for Quantum Espresso

# Requirements
Python package PYQE has been installed to your machine. See INSTALL in
home directory for information on how to do this (hint: it is
installed like every other python package).

# Direct Interfaces (Currently Only One PWBase) 
PYQE was designed to have two levels of interfaces. The first level, direct
interfaces were designed to preserve all the features and
configurations available in the input files for a given
tool/command. Currently only `pw.x` the plane wave code is implemented
however with the current framework it will be very easy to add
additional commands and fit them into one cohesive package. 

The pw.x interface is `pyqe.espresso.PWBase`. An example of its usage
can be found in `example/qe/si.py`. This file is meant to mimic the
example in quantum espresso `PW/examples/example01/`.

## Why use python interface instead of input files?
There are many added benefits to using the PYQE package over a regualr
input file.  
 - PYQE is able to validate the input before actually doing a run by
   doing range, domain, and can check if the configuration makes
   sense (e.g. does the directory outdir exist?).  
 - Automatically read output files and save files.  
 - Anything valid in python is valid in pyqe!  
It should also be mentioned that since pyqe can read in input and
output files. All of the analysis of the dft runs can easily be done
with python.

# ASE interfaces 
Often times the direct interface approach is cumbersome in that the
python files can be very long. The python Atomic Simulation
Environment (ASE) is a package that provides a unifying approach of
using MD and DFT codes. ASE allows for the definition of `Calculators`
which provide an interface for interacting with DFT codes in a
standard way. For more information consult the ase
[webpage](https://wiki.fysik.dtu.dk/ase/index.html). PYQE implements
an ASE calculator interface `pyqe.ase.qe.QE`. This interface provides
a simple interface to using Quantum Espresso. *This interface is
recommended*. The PYQE ASE Calculator makes use of all of the Base
direct implementations (PWBase) only 1 currently.

# Parallel Execution
When PYQE executes a run it runs pw.x as:
```
<prefix> pw.x < stdin > stdout <postfix>
```

*prefix* and *postfix* can be modified using `pyqe.config.prefix` and
`pyqe.config.postfix`. See `pyqe.config` for more adjustable
parameters. stdin and stdout are handled by pyqe and is how pyqe can
easily read the results from the run.

Thus the following would configure a parallel run on 16 processors. If
you are using a Job Scheduler make sure to request the proper number
of processors!

```
import pyqe.config as config
config.prefix = ['mpirun', '-np', '16']
```

