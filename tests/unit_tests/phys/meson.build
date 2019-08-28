exec_dsf1 = executable('dsf1.x',
                      ['test_dsf1.cpp', env, time, lattice, basis, hamilt, solver, correl],
                      include_directories: inc,
		      cpp_args:['-I'+boost],
                      dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_dsf2 = executable('dsf2.x',
                      ['test_dsf2.cpp', env, time, lattice, basis, hamilt, solver, correl],
                      include_directories: inc,
                      cpp_args:['-I'+boost],
                      dependencies: [gtest_dep, slepc_dep, mpi_dep])

compare = find_program('../../compare.sh')
_prefix = meson.current_source_dir()

_compare_tests = ['dsf']
_cases = ['1','2']

foreach _t: _compare_tests
  foreach _num: _cases
     tmp = get_variable('exec_' + _t + _num)
     test(_t + _num, 
          compare, 
          args:[tmp.full_path(), '2',
                _prefix + '/../reference_data/dsf_d0_' + _num + '.ref',
                _prefix + '/dsf' + _num + '/dsf_d0.dat',
                '-malloc_dump'],
	  is_parallel: true,
          workdir: meson.current_source_dir(),
          suite: ['compare_tests', _t])
  endforeach    
endforeach