exec_build_off_diag125 = executable('build_off_diag125.x',
                                   ['test_build_off_diag125.cpp', env, time, lattice, basis, hamilt],
                                   include_directories: inc,
                                   dependencies: [petsc_dep, slepc_dep])

exec_build_off_diag82 = executable('build_off_diag82.x',
                                  ['test_build_off_diag82.cpp', env, time, lattice, basis, hamilt],
                                  include_directories: inc,
                                  dependencies: [petsc_dep, slepc_dep])

compare = find_program('../../../compare.sh')
_prefix = meson.current_source_dir()

_compare_tests = ['build_off_diag']
_cases = ['125', '82']

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