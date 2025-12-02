exec_srcs = files(
    'main.cpp',
)

includes = include_directories('../include')

# -- DEPENDENCIES
exec_deps = []

# add fmt
fmt_proj = subproject('fmt')
fmt_dep = fmt_proj.get_variable('fmt_dep')
exec_deps += fmt_dep

# add armadillo
armadillo_proj = subproject('armadillo', default_options: ['use_lapack=true', 'use_blas=true'])
armadillo_dep = armadillo_proj.get_variable('armadillo_dep')
exec_deps += armadillo_dep

if cxx.get_id() == 'clang'
    exec_deps +=  declare_dependency(
        compile_args: '-fopenmp',
        link_args: '-fopenmp',
    )
else
    exec_deps += dependency('openmp', language: 'cpp')
endif

# linear algebra backend
foreach pkg: get_option('la_packages')
    exec_deps += dependency(pkg, required: true)
endforeach

foreach lib: get_option('la_libraries')
    exec_deps += cc.find_library(lib, required : true)
endforeach

# -- EXECUTABLE
mcpp_heisenberg = executable(
    meson.project_name(),
    exec_srcs,
    cpp_args : build_args,
    include_directories : includes,
    dependencies: exec_deps,
    install : true
)
