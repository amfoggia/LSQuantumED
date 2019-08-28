exec_build_diag125h = executable('build_diag125.x',
                                ['test_build_diag125.cpp', env, time, lattice, basis, hamilt],
                                include_directories: inc,
		                cpp_args:['-I'+boost],
                                dependencies: [slepc_dep, mpi_dep])

exec_build_diag82h = executable('build_diag82.x',
                               ['test_build_diag82.cpp', env, time, lattice, basis, hamilt],
                               include_directories: inc,
                               cpp_args:['-I'+boost],
                               dependencies: [slepc_dep, mpi_dep])

compare = find_program('../../../compare.sh')
_prefix = meson.current_source_dir()

_compare_tests = ['build_diag']
_cases = ['125h', '82h']

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