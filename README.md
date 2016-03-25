Description
======================
Rabacus is a `Python <http://www.python.org>`_  package for performing 
analytic radiative transfer calculations in simple geometries relevant to 
cosmology and astrophysics. It also contains tools to calculate cosmological
quantities such as the power spectrum and mass function.  


Prerequisites
======================

The Rabacus package requires three other Python packages and a fortran 
compiler, 

- `Scipy <http://www.scipy.org/scipylib/index.html>`_
- `Numpy <http://www.numpy.org>`_ (version 1.7 or later)
- `Quantities <https://pythonhosted.org/quantities>`_
- Fortran compiler in your path.    


Installing prerequisites with pip
------------------------------------

A simple way to install Python packages is using the package manager 
`pip <https://pypi.python.org/pypi/pip>`_.  To check if you have
pip installed on your system, type ``pip`` at the command line, :: 

  pip

If this produces usage instructions then congratulations, you have pip
installed.  If not, instructions for installing pip can be found 
`here <http://www.pip-installer.org/en/latest/installing.html>`_.
To check if the python packages are installed on your system, attempt 
to import them from the python command prompt, ::

  >>> import scipy
  >>> import numpy as np
  >>> import quantities as pq

If any of these import commands produce an error message you will need 
to install the proper software before installing Rabacus. 
Once you have access to pip, you can install any missing prerequisites 
using the following commands, ::

  sudo pip install scipy
  sudo pip install numpy 
  sudo pip install quantities

If you do not have root access on your system you can pass the ``--user`` flag 
which will install the packages into a hidden folder called ``.local`` in 
your home directory, ::

  pip install --user scipy
  pip install --user numpy 
  pip install --user quantities


Installing prerequisites on Debian (Ubuntu)
---------------------------------------------

On Debian based systems (such as Ubuntu) you may prefer to install these 
prerequisites using the APT tool, ::
 
  sudo apt-get install python-scipy python-numpy python-quantities

To increase the speed of execution, much of Rabacus is written in 
Fortran 90 and then wrapped using the f2py tool that is part of numpy.  For 
the installation to be successful, a fortran compiler must be in your 
executable path.  If you don't already have one, I recommend the 
gnu fortran compiler `gfortran <http://gcc.gnu.org/wiki/GFortran>`_.     
On Debian based systems (such as Ubuntu) you can install this 
compiler using the APT tool, ::

  sudo apt-get install gfortran

Installation
======================

With the prequisites installed on your system, you are ready to
install the Rabacus package itself.  

Setting ``F90``  environment variable
-------------------------------------

Rabacus makes use of OpenMP directives in the Fortran code base and so
we have to make sure the code is compiled correctly.  In order to do
this, you have to let the build system know what Fortran 90 compiler
you are going to be using.  The simplest way to do this is to set the
environment variable ``F90`` before following the installation
instructions below.  Rabacus has been tested with the intel compiler and
the gnu gfortran compiler.  For other compilers you will have to
follow the `Manual Install` instructions below.

To use the gfortran compiler, type the following at the command line
(in Bash) ::

  export F90=gfortran

To use the intel compiler set :: 

  export F90=ifort


Single command install
-------------------------

If you have made the appropriate sacrifices to the computer gods, you
should be able to install an OpenMP enabled version of Rabacus with a
single comand line call to pip, ::

  sudo pip install rabacus

As was the case for the prerequisites, if you do not have root access
on your system you can pass the ``--user`` flag which will install
Rabacus into a hidden folder called ``.local`` in your home directory,
::

  pip install --user rabacus


If the last two lines printed to the screen are, :: 

  Successfully installed rabacus 
  Cleaning up...

then congratulations you have a working copy of Rabacus. To double
check, begin an ipython session and attempt an import, ::

  import rabacus as ra

Packages installed with pip can be uninstalled in the same way, ::

  pip uninstall rabacus


Manual install
------------------------

If the above process fails for any reason we can always download
Rabacus and manually invoke the setup script.  The first step is to
download and untar the Rabacus tar.gz file from the PyPI site
(https://pypi.python.org/pypi/rabacus) and change into the main
Rabacus directory, ::

  gunzip rabacus-x.x.x.tar.gz
  tar xvf rabacus-x.x.x.tar 
  cd rabacus-x.x.x

Now we have direct access to the ``setup.py`` file which gives us a
lot more freedom but it comes at the cost of slightly more complexity.
First it's a good idea to see which fortran compilers are detected on
your machine.  The following command will list all of the fortan
compilers found on your system and all the compilers available for
your system but not found.  ::

  f2py -c --help-fcompiler

For example on my machine I get the following, ::

  Fortran compilers found:
    --fcompiler=gnu95    GNU Fortran 95 compiler (4.8.1-10)
    --fcompiler=intelem  Intel Fortran Compiler for 64-bit apps (14.0.2.144)
  Compilers available for this platform, but not found:
    --fcompiler=absoft   Absoft Corp Fortran Compiler
    --fcompiler=compaq   Compaq Fortran Compiler
    --fcompiler=g95      G95 Fortran Compiler
    --fcompiler=gnu      GNU Fortran 77 compiler
    --fcompiler=intel    Intel Fortran Compiler for 32-bit apps
    --fcompiler=intele   Intel Fortran Compiler for Itanium apps
    --fcompiler=lahey    Lahey/Fujitsu Fortran 95 Compiler
    --fcompiler=nag      NAGWare Fortran 95 Compiler
    --fcompiler=pathf95  PathScale Fortran Compiler
    --fcompiler=pg       Portland Group Fortran Compiler
    --fcompiler=vast     Pacific-Sierra Research Fortran 90 Compiler

Now we decide which of the fortran compilers to use and which flags to
pass the build command.  Suppose you wanted to use the Intel compiler.
Edit the ``setup.py`` file such that the variable ``f90_flags`` is a
list of compile flags and ``omp_lib`` is a list containing the linking
flags.  For example, ::

  f90_flags = ["-openmp", "-fPIC", "-xHost", "-O3", "-ipo", 
               "-funroll-loops", "-heap-arrays", "-mcmodel=medium"]  

  omp_lib = ["-liomp5"]

These variables are already defined near the top of the ``setup.py``
file and will need to be overwritten.  Once this is done, we give the
build command to the ``setup.py`` script, ::

  python setup.py build --fcompiler=intelem

After the package is built, give the install command to actually
install it, ::

  sudo python setup.py install --record rabacus_install_files.txt

The last part of the command is to allow for easy uninstall.  This
process just involves deleting all installed files which will be
listed in the file ``rabacus_install_files.txt``. This can be
accomplished using the following command, ::

  cat rabacus_install_files.txt | xargs sudo rm -rf

The install can also be done locally for those without root permission
on their system by passing the ``--user`` flag to the install command,
::

  python setup.py install --user --record rabacus_install_files.txt

Note that if you previously did an install of Rabacus that required
the ``sudo`` command you will likely need to delete the
``rabacus.egg-info`` directory and some directories inside the
``build`` directory as they will need to be modified but will be owned
by ``root``.  If you are only doing a local install then this shoudn't
be necessary.  This procedure should work for any fortran compiler
supported by f2py (i.e. any compiler in the list returned when using
the ``--help-fcompiler`` flag.


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
  print hm12.He1h(z)

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

* PyPI (https://pypi.python.org/pypi/rabacus) 
* documentation (http://pythonhosted.org//rabacus)
* version control (https://bitbucket.org/galtay/rabacus)


