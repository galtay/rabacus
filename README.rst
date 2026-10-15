Description
======================
Rabacus is a `Python <http://www.python.org>`_  package for performing
analytic radiative transfer calculations in simple geometries relevant to
cosmology and astrophysics. It also contains tools to calculate cosmological
quantities such as the power spectrum and mass function.


Prerequisites
======================

Rabacus requires Python 3.12 or later and a Fortran compiler.  Much of
Rabacus is written in Fortran 90 and wrapped for Python using
`f2py <https://numpy.org/doc/stable/f2py/>`_, which is part of numpy.
The extension is built with `meson-python
<https://mesonbuild.com/meson-python/>`_.

We recommend `uv <https://docs.astral.sh/uv/>`_ for managing the Python
environment and the gnu fortran compiler
`gfortran <https://gcc.gnu.org/wiki/GFortran>`_.

On macOS, ``gfortran`` is part of the Homebrew ``gcc`` package, ::

  brew install gcc

On Debian based systems (such as Ubuntu), ::

  sudo apt-get install gfortran


Installation
======================

Development install with uv
---------------------------

Clone the repository and let uv create the environment, build the Fortran
extension, and install Rabacus in editable mode, ::

  git clone https://github.com/galtay/rabacus.git
  cd rabacus
  uv sync

The editable install rebuilds the Fortran extension automatically the next
time ``rabacus`` is imported after a source file changes.  To run the test
suite, ::

  uv run pytest            # fast tests (a few seconds)
  uv run pytest -m slow    # Python vs. Fortran solver comparisons

Install with pip
----------------

Rabacus can also be installed into any environment with pip, ::

  pip install .

OpenMP
------

The Fortran code contains OpenMP directives and is compiled with OpenMP
support when the compiler provides it.  To build without OpenMP, ::

  pip install . -Csetup-args=-Dopenmp=disabled

Choosing a Fortran compiler
---------------------------

Meson uses the compiler named by the ``FC`` environment variable, or the
first Fortran compiler it finds on the ``PATH``.  For example, to use
the Intel compiler, ::

  FC=ifx uv sync --reinstall-package rabacus


Testing install
------------------------

Detailed examples of using rabacus are available by following the link
to the users guide below.  However, we present a short example with the
expected output below as a way to quickly test that a new installation
has basic functionality.  We first import rabacus and then create an
object that gives access to the meta galactic radiation background
described in
`Haardt & Madau 2012
<http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`_.  Finally, we
ask for the photo-heating rate of He I at a redshift of 3.0.  ::

  import rabacus as ra
  hm12 = ra.HM12_Photorates_Table()
  z = 3.0
  print(hm12.He1h(z))

The expected output from a working rabacus installation is given
below.  Note that there may be differences in the last significant
figure due to different processor architectures. ::

  3.39163517433e-12 eV/s



Author
=====================
Rabacus was written by Gabriel Altay and any questions can be directed
to gabriel.altay@gmail.com



Project URLs
=====================

* version control (https://github.com/galtay/rabacus)
* PyPI (https://pypi.org/project/rabacus)


License
=====================
Rabacus is released under the MIT license (see ``LICENSE``).  The file
``rabacus/f2py/legendre_polynomial.f90`` is third-party code by John
Burkardt distributed under the GNU LGPL.
