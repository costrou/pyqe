# PYQE Interfaces for Quantum Espresso

# Requirements
Python package PYQE has been installed to your machine. See INSTALL in
home directory for information on how to do this (hint: it is
installed like every other python package).

# QE PWBase/ASE interfaces 
When creating a python interface for Quantum Espresso effort was made
to preserve all the features and configurations available in the input
files for pw.x. However there are many added benefits to using the
pyqe package over a regular input file. PYQE is able to validate the
input file before actually doing a run along with automatically
reading in the output file and save files. An example for PW/example01
with Silicon using the PWBase interface `pyqe.espresso.PWBase` can be
seen in `qe/si.py`.

Often times this approach is cumbersome in that the python files can
be very long. PYQE also implements an ASE calculator interface
`pyqe.ase.qe.QE`. Therefore you can easily view your systems using ASE
and use all of its additional functionality. This interface provides a
simple interface to using QE. *This interface is recommended*. 

# Parallel Execution
When PYQE executes a run it runs pw.x as:
    <prefix> pw.x < stdin > stdout <postfix>

<prefix> and <postfix> can be modified using `pyqe.config.prefix` and
`pyqe.config.postfix`. See config for more adjustable
parameters. stdin and stdout are handled by pyqe.

Thus the following would configure a parallel run on 16 processors.

```
import pyqe.config as config
config.prefix = ['mpirun', '-np', '16']
```

