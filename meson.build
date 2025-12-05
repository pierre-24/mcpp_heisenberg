project(
    'mcpp_heisenberg', 'cpp',
    version : '0.1.2',
    default_options : ['warning_level=3', 'cpp_std=c++20']
)

cxx = meson.get_compiler('cpp')

build_args = [
    '-DPROJECT_NAME="' + meson.project_name() + '"',
    '-DPROJECT_VERSION="' + meson.project_version() + '"',
]

date_program = find_program('date', required: false)
if date_program.found()
    build_args += '-DPROJECT_BUILD_DATE="' + run_command(date_program).stdout().replace('\n', '') + '"'
endif

git_program = find_program('git', required: false)
if git_program.found()
    build_args += '-DPROJECT_BUILD_COMMIT="' + run_command(git_program, 'rev-parse', 'HEAD').stdout().substring(0, 7) + '"'
endif

subdir('lib')

subdir('app')

# -- TESTS
if get_option('tests')
    subdir('tests')
endif