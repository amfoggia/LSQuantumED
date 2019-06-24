project('QuantumMagnetsLSSimulation', 'cpp', default_options : ['cpp_std=c++17'])

# Use second option when running on Marconi
#inc = include_directories(['.','include','include/spinOperators','include/tools','/opt/lib/boost_1_67_0'])
inc = include_directories(['.','include','include/spinOperators','include/tools'])
#inc = include_directories(['.','include','include/spinOperators','include/tools','/cineca/prod/opt/libraries/boost/1.66.0/intelmpi--2018--binary/include'])

conf_data = configuration_data()
conf_data.set_quoted('version', meson.project_version())
conf_data.set('DEVEL', get_option('DEVEL'))
conf_data.set('TIME_CODE', get_option('TIME_CODE'))
conf_data.set('DISORDER', get_option('DISORDER'))
configure_file(input: 'include/meson_config.h.in',
	       output: 'meson_config.h',
	       configuration: conf_data)

petsc_dep = dependency('PETSc')
slepc_dep = dependency('SLEPc')

# Comment when running on marconi
mpi_dep = dependency('mpich')

subdir('src')

main_src = [env,
	    basis,
	    lattice,
	    sublat,
	    hamilt,
	    solver,
	    rand_dis,
	    time,
	    correl,
	    phys]

if get_option('DISORDER') == false
   exec_main_NO_dis = executable('main_NO_disorder.x',
                                ['main_NO_disorder.cpp', main_src],
                                include_directories: inc,
                                dependencies: [petsc_dep, slepc_dep, mpi_dep])
endif

if get_option('DISORDER') == true
   exec_main_dis = executable('main_disorder.x',
   		             ['main_disorder.cpp', main_src],
              	             include_directories: inc,
                             dependencies: [petsc_dep, slepc_dep, mpi_dep])
   exec_simple_main = executable('simple_main.x',
   		                ['simple_main.cpp', main_src],
              	                include_directories: inc,
                                dependencies: [petsc_dep, slepc_dep, mpi_dep])
endif

exec_main = executable('main.x',
                      ['main.cpp', main_src],
                      include_directories: inc,
                      dependencies: [petsc_dep, slepc_dep, mpi_dep])

subdir('tests')
