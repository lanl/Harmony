
add_executable( create_mesh
  create_mesh.c
)
add_compile_options(-lm)
#target_link_libraries(<math.h> m)
target_link_libraries(create_mesh m)
