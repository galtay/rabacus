project(
  'rabacus',
  'c', 'fortran',
  version: '0.10.0',
  license: 'MIT',
  meson_version: '>=1.3.0',
  default_options: ['buildtype=release', 'b_ndebug=if-release'],
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Include paths are made relative to the source dir so that a virtual
# environment living inside the source tree (e.g. .venv) is accepted.
incdirs = run_command(py,
  ['-c', '''
import os, sys, numpy, numpy.f2py
for d in (numpy.get_include(), numpy.f2py.get_include()):
    try:
        d = os.path.relpath(d, sys.argv[1])
    except ValueError:
        pass
    print(d)
''', meson.current_source_dir()],
  check: true,
).stdout().strip().split('\n')
incdir_numpy = incdirs[0]
incdir_f2py = incdirs[1]

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Fortran sources wrapped into the rabacus.f2py.rabacus_fc extension.
# Meson resolves Fortran module dependencies, so order does not matter.
fortran_sources = files(
  'rabacus/f2py/types.f90',
  'rabacus/f2py/utils.f90',
  'rabacus/f2py/zhang_jin.f90',
  'rabacus/f2py/zhang_jin_f.f90',
  'rabacus/f2py/slatec.f90',
  'rabacus/f2py/special_functions.f90',
  'rabacus/f2py/legendre_polynomial.f90',
  'rabacus/f2py/m_mrgrnk.f90',
  'rabacus/f2py/physical_constants.f90',
  'rabacus/f2py/geometry.f90',
  'rabacus/f2py/hui_gnedin_97.f90',
  'rabacus/f2py/hm12.f90',
  'rabacus/f2py/chem_cool_rates.f90',
  'rabacus/f2py/ion_solver.f90',
  'rabacus/f2py/verner_96.f90',
  'rabacus/f2py/photo_xsections.f90',
  'rabacus/f2py/spectra.f90',
  'rabacus/f2py/source_point.f90',
  'rabacus/f2py/source_plane.f90',
  'rabacus/f2py/source_background.f90',
  'rabacus/f2py/sphere_base.f90',
  'rabacus/f2py/sphere_stromgren.f90',
  'rabacus/f2py/sphere_bgnd.f90',
  'rabacus/f2py/slab_base.f90',
  'rabacus/f2py/slab_plane.f90',
  'rabacus/f2py/slab_bgnd.f90',
)

f2py_wrappers = custom_target(
  'rabacus_fc_f2py',
  input: fortran_sources,
  output: ['rabacus_fcmodule.c', 'rabacus_fc-f2pywrappers2.f90'],
  command: [py, '-m', 'numpy.f2py', '@INPUT@',
            '-m', 'rabacus_fc', '--lower', '--build-dir', '@OUTDIR@'],
)

# OpenMP is applied to the Fortran code only; the f2py C wrapper does
# not need it (and Apple clang does not support -fopenmp).
omp_dep = dependency('openmp', language: 'fortran', required: get_option('openmp'))

fortran_lib = static_library(
  'rabacus_fortran',
  fortran_sources,
  dependencies: omp_dep,
  pic: true,
)

py.extension_module(
  'rabacus_fc',
  [f2py_wrappers, incdir_f2py / 'fortranobject.c'],
  include_directories: inc_np,
  link_with: fortran_lib,
  dependencies: [py_dep, omp_dep.partial_dependency(link_args: true)],
  link_language: 'fortran',
  install: true,
  subdir: 'rabacus/f2py',
)

install_subdir(
  'rabacus',
  install_dir: py.get_install_dir(),
  exclude_directories: ['__pycache__'],
)
