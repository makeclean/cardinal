[Mesh]
  [structured_mesh]
    type = CartesianMeshGenerator
    dim = 3
    dx = '1.0 1.0 3.0 5.0 90.0'
    dy = '1.0 1.0 3.0 5.0 90.0'
    dz = '1.0 1.0 3.0 5.0 90.0'
  []
  allow_renumbering = false
[]


[Problem]
  type = OpenMCCellAverageProblem
  verbose = true
  cell_level = 0

  power = 1000.0
  volume_calculation = vol

  [Tallies]
    [heat_source]
      type = StructuredMeshTally
      dimensions = '3'
      mesh_type = rectilinear
      lower_left = '0 0 0'
      x_grid = '0.0 1.0 2.0 5.0 10.0 100.0'
      y_grid = '0.0 1.0 2.0 5.0 10.0 100.0'
      z_grid = '0.0 1.0 2.0 5.0 10.0 100.0'

      name = heat_source
      score = 'flux'
      check_tally_sum = false
      normalize_by_global_tally = false
    []
  []
[]

[UserObjects]
  [vol]
    type = OpenMCVolumeCalculation
    n_samples = 200000
  []
[]

[Postprocessors]
  [k]
    type = KEigenvalue
    value_type = 'tracklength'
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  csv = true
[]
