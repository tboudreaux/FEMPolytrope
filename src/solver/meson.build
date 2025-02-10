exFile = 'laneEmden.cpp'
exName = 'laneEmden'

executable(
  exName,
  exFile,
  link_with: [libcoeff, libutils],
  dependencies: [mfem_dep],
  include_directories: ['../coeff/', '../utils/'],
  install_rpath: '@loader_path/../../src'
)
