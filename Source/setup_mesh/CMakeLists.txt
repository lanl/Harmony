
add_executable(setup_mesh
  setup_mesh.cpp
)
add_compile_options(-lm)
#target_link_libraries(<math.h> m)
target_link_libraries(setup_mesh m)
