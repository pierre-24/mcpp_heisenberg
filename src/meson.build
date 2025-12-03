lib_srcs = files(
    'hamiltonian.cpp',
    'logging.cpp',
    'geometry.cpp',
    'mc.cpp',
)

includes = include_directories('../include')

# -- DEPENDENCIES
lib_deps = []

# add fmt
fmt_proj = subproject('fmt')
fmt_dep = fmt_proj.get_variable('fmt_dep')
lib_deps += fmt_dep

# add plog
plog_proj = subproject('sergiusthebest-plog')
plog_dep = plog_proj.get_variable('plog_dep')
lib_deps += plog_dep

build_args += '-DPLOG_CAPTURE_FILE'

# add highfive
highfive_proj = subproject('highfive')
highfive_dep = highfive_proj.get_variable('highfive_dep')
lib_deps += highfive_dep

# add armadillo
armadillo_proj = subproject('armadillo', default_options: ['use_lapack=true', 'use_blas=true'])
armadillo_dep = armadillo_proj.get_variable('armadillo_dep')
lib_deps += armadillo_dep

if cxx.get_id() == 'clang'
    lib_deps +=  declare_dependency(
        compile_args: '-fopenmp',
        link_args: '-fopenmp',
    )
else
    lib_deps += dependency('openmp', language: 'cpp')
endif

# linear algebra backend
foreach pkg: get_option('la_packages')
    lib_deps += dependency(pkg, required: true)
endforeach

foreach lib: get_option('la_libraries')
    lib_deps += cc.find_library(lib, required : true)
endforeach

# -- MAKE LIBRARY
mcpp_heisenberg_lib = library(
    meson.project_name(),
    lib_srcs,
    install : true,
    cpp_args : build_args,
    include_directories : includes,
    dependencies: lib_deps,
)

lmcpp_heisenberg_dep = declare_dependency(
    link_with: mcpp_heisenberg_lib,
    include_directories: includes,
    dependencies: lib_deps
)

# -- EXECUTABLE
mcpp_heisenberg = executable(
    meson.project_name() + '_run',
    lib_srcs + 'main.cpp',
    include_directories : includes,
    dependencies: lib_deps + lmcpp_heisenberg_dep,
    cpp_args : build_args,
    install : true
)
