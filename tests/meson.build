# dependencies
tests_deps = [lmcpp_heisenberg_dep]

gtest_proj = subproject('gtest')
gtest_dep = gtest_proj.get_variable('gtest_dep')

tests_deps += gtest_dep

# create tests suite
tests = [
    'configs',
    'hamiltonian',
    'geometry',
    'mc',
]

tests_srcs = []

foreach t: tests
    tests_srcs += files('tests_' + t + '.cpp')

    test_exec = executable(
        'tests_suite_' + t,
        sources: [files('tests_' + t + '.cpp'), 'main.cpp'],
        include_directories : include_directories('.'),
        dependencies: tests_deps,
        install: false,
        build_by_default: false,
        cpp_args : build_args,
    )

    test('tests_suite_' + t,
         test_exec,
         workdir : meson.current_source_dir(),
         suite: 'lib_code',
    )
endforeach

# test app
test_app_deps = [lmcpp_heisenberg_run_dep, gtest_dep]

test_app_exec = executable(
    'tests_suite_app',
    sources: ['tests_app.cpp', 'main.cpp'],
    include_directories : include_directories('.'),
    dependencies: test_app_deps,
    install: false,
    build_by_default: false,
    cpp_args : build_args,
)

test('tests_suite_app',
     test_app_exec,
     workdir : meson.current_source_dir(),
     suite: 'app_code',
)