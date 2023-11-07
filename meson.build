project(
  'olptm', 
  'c', 
  version: '0.0.0', 
  license: 'MIT',
)

add_global_arguments(language: 'c')

cc = meson.get_compiler('c')
math_dep = cc.find_library('m', required: true)
inc = include_directories('include')

# subdir('include')
# subdir('src')
# subdir('demo')

ltm_src = ['src/ltm.c']
install_headers('include/ltm.h', subdir: 'olptm')
install_headers('include/constants.h', subdir: 'olptm')
ltmlib = shared_library('ltm', ltm_src, include_directories: inc, dependencies: math_dep, install: true)
demoexe = executable('demo', 'demo/demo.c', include_directories: inc, link_with: ltmlib)
