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

[AuxVariables]
  [total_flux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [sum_aux]
    type = ArrayVarReductionAux
    variable = total_flux
    array_variable = flux
    value_type = sum
  []
[]
    
[Problem]
  type = OpenMCCellAverageProblem
  verbose = true
  cell_level = 0
  
  source_strength = 1.0e18
  volume_calculation = vol

  [Tallies]
    [Flux]
      type = StructuredMeshTally
      dimensions = '3'
      mesh_type = rectilinear
      filters = 'Energy'
      lower_left = '0 0 0'
      x_grid = '0.0 1.0 2.0 5.0 10.0 100.0'
      y_grid = '0.0 1.0 2.0 5.0 10.0 100.0'
      z_grid = '0.0 1.0 2.0 5.0 10.0 100.0'
      name = flux
      estimator = 'tracklength'
      score = 'flux'
      add_energy_array = true
      check_tally_sum = false
      normalize_by_global_tally = false
      output = unrelaxed_tally 
    []
  []
  [Filters]
    [Energy]
      type = EnergyFilter
      # CASMO 2 group structure for testing. May result in some missed particles
      energy_boundaries = '0.0 6.25e-1 2.0e7'
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
  # Group fluxes obtained by integrating each component of the 'flux' array auxvariable.
  [Total_Flux_1]
    type = ElementIntegralArrayVariablePostprocessor
    variable = flux
    component = 0
  []
  [Total_Flux_2]
    type = ElementIntegralArrayVariablePostprocessor
    variable = flux
    component = 1
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  csv = true
[]
