#!/usr/bin/python3
from distutils.core import setup

import sys
# import os

long_description = """
pyqe is a python package providing an open source python api for
controlling Quantum Espresso in the Python Language.
"""

if sys.version_info < (3, 0, 0, 'final', 0):
    raise SystemExit('Python 3.0 or later is required!')

# packages = []
# for dirname, diranmes, filenames in os.walk('pyqe'):
#     if '__init__.py' in filenames:
#         packages.append(dirname.replace('/', '.'))

package_dir = {'pyqe': 'pyqe'}

setup(name='python-qe',
      version="0.0.1",
      description='Quantum Espresso Python Toolkit',
      url='http://git.aves.io/costrouc/pyqe/wikis/home',
      maintainer='Christopher Ostrouchov',
      maintainer_email='chris.ostrouchov+pyqe@gmail.com',
      license='LGPLv2.1+',
      platforms=['linux'],
      package_dir=package_dir,
      long_description=long_description)
