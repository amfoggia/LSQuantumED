exec_build_off_diag60 = executable('build_off_diag60.x',
                                  ['test_build_off_diag60.cpp', env, time, lattice, basis, hamilt],
                                  include_directories: inc,
                                  dependencies: [petsc_dep, slepc_dep, mpi_dep])

exec_build_off_diag61 = executable('build_off_diag61.x',
                                  ['test_build_off_diag61.cpp', env, time, lattice, basis, hamilt],
                                  include_directories: inc,
                                  dependencies: [petsc_dep, slepc_dep, mpi_dep])

exec_build_off_diag62 = executable('build_off_diag62.x',
                                  ['test_build_off_diag62.cpp', env, time, lattice, basis, hamilt],
                                  include_directories: inc,
                                  dependencies: [petsc_dep, slepc_dep, mpi_dep])

exec_build_off_diag63 = executable('build_off_diag63.x',
                                  ['test_build_off_diag63.cpp', env, time, lattice, basis, hamilt],
                                  include_directories: inc,
                                  dependencies: [petsc_dep, slepc_dep, mpi_dep])

exec_build_off_diag6m1 = executable('build_off_diag6m1.x',
                                   ['test_build_off_diag6m1.cpp', env, time, lattice, basis, hamilt],
                                   include_directories: inc,
                                   dependencies: [petsc_dep, slepc_dep, mpi_dep])

exec_build_off_diag6m2 = executable('build_off_diag6m2.x',
                                   ['test_build_off_diag6m2.cpp', env, time, lattice, basis, hamilt],
                                   include_directories: inc,
                                   dependencies: [petsc_dep, slepc_dep, mpi_dep])

exec_build_off_diag6m3 = executable('build_off_diag6m3.x',
                                   ['test_build_off_diag6m3.cpp', env, time, lattice, basis, hamilt],
                                   include_directories: inc,
                                   dependencies: [petsc_dep, slepc_dep, mpi_dep])

compare = find_program('../../../compare.sh')
_prefix = meson.current_source_dir()

_compare_tests = ['build_off_diag']
_cases = ['60', '61', '62', '6m1', '6m2']

foreach _t: _compare_tests
  foreach _num: _cases
     tmp = get_variable('exec_' + _t + _num)
     test(_t + _num, 
          compare, 
          args:[tmp.full_path(), '2',
                _prefix + '/../../reference_data/' + _t + _num + '.ref',
                _prefix + '/' + _t + _num + '.dat',
                '-malloc_dump'],
          workdir: meson.current_source_dir(),
	  suite: ['compare_tests', _t])
  endforeach	
endforeach

_cases = ['63', '6m3']

foreach _t: _compare_tests
  foreach _num: _cases
     tmp = get_variable('exec_' + _t + _num)
     test(_t + _num, 
          compare, 
          args:[tmp.full_path(), '1',
                _prefix + '/../../reference_data/' + _t + _num + '.ref',
                _prefix + '/' + _t + _num + '.dat',
                '-malloc_dump'],
          workdir: meson.current_source_dir(),
	  suite: ['compare_tests', _t])
  endforeach	
endforeach