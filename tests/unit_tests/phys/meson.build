exec_dsf1 = executable('dsf1.x',
                      ['test_dsf1.cpp', env, time, lattice, basis, hamilt, solver, correl, phys],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

# exec_dsf2 = executable('dsf2.x',
#                       ['test_dsf2.cpp', env, time, lattice, basis, hamilt, solver, correl, phys],
#                       include_directories: inc,
#                       dependencies: [gtest_dep, petsc_dep, slepc_dep])


compare = find_program('../../compare.sh')
_prefix = meson.current_source_dir()

_compare_tests = ['dsf']
_cases = ['1']

foreach _t: _compare_tests
  foreach _num: _cases
     tmp = get_variable('exec_' + _t + _num)
     test(_t + _num, 
          compare, 
          args:[tmp.full_path(), '2',
                _prefix + '/../reference_data/dsf_d0.ref',
                _prefix + '/dsf_d0.dat',
                '-malloc_dump'],
          workdir: meson.current_source_dir(),
          suite: ['compare_tests', _t])
  endforeach    
endforeach


# _tests = ['dsf']
# _cases = ['1']

# _prefix = meson.current_source_dir()

# mpiexec = find_program('mpiexec')

# foreach _t: _tests
#   foreach _num: _cases
#     tmp = get_variable('exec_' + _t + _num)
#     test(_t,
#       mpiexec,
#       args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
#       workdir: meson.current_source_dir(),
#       suite: ['unit_tests', _t+_num])
#     endforeach
# endforeach