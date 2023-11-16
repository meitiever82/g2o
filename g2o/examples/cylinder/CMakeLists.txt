add_executable(fit_cylinder
  fit_cylinder.cpp
)

set_target_properties(fit_cylinder PROPERTIES OUTPUT_NAME fit_cylinder)
target_link_libraries(fit_cylinder core types_slam3d solver_dense)
