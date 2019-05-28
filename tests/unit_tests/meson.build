gtest_main_dep = dependency('gtest', main: true, 
                            fallback : ['gtest', 'gtest_main_dep'],  # magic of meson wrap
                            required: true)


gtest_dep = dependency('gtest', main:false, 
                       fallback : ['gtest', 'gtest_dep'],  # magic of meson wrap
                       required: true)

exec_construct_hamilt = executable('construct_hamilt.x',
                                  ['test_construct_hamilt.cpp', env, time, lattice, basis, hamilt],
                                  include_directories: inc,
                                  dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

_tests = ['construct_hamilt']

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
  tmp = get_variable('exec_' + _t)
  test(_t,
       mpiexec,
       args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
       workdir: meson.current_source_dir(),
       suite: ['unit_tests', _t])
endforeach

subdir('environment')

if get_option('DEVEL') == true
   subdir('lattice')
   subdir('spinOPft')
endif

subdir('sublattice')

subdir('basis')

subdir('get_elem_diag')
subdir('get_coup_elems')
subdir('build_diag')
subdir('build_off_diag')

if get_option('DISORDER') == true
   subdir('tools_random')
   subdir('get_disorder_elem')
   subdir('build_disorder')
endif

subdir('prealloc_info')

subdir('solver')

subdir('correlation')

subdir('phys')