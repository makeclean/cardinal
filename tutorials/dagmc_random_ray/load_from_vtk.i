[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = 'plot_1.vtk'
  []
[]

[AuxVariables]
  [flux1]
    family = MONOMIAL
    order = CONSTANT
    initial_value_from_file = 'flux_group_0'
  []
[]
