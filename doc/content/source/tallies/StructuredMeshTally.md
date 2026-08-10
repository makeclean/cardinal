# StructuredMeshTally
  id=ssm

!alert note
Structured mesh tallies are not currently supported when running neutronics with OpenMC's
random ray solver.

## Description

The `StructuredMeshTally` class wraps an OpenMC tally with a **structured** (regular or
rectilinear) mesh spatial filter. Like the [MeshTally](MeshTally.md), the OpenMC structured
mesh is generated "on the fly" from text parameters, and the tally results are written into
monomial (constant) MOOSE auxiliary variables on the `[Mesh]` block. The `[Mesh]` block must
therefore describe the same grid (same bounds, cell counts, and x-fastest element ordering).
Because the OpenMC mesh is generated from text parameters, it need not be provided as a file.

The structured mesh is described by per-axis node coordinates. Both the OpenMC mesh and the
native libMesh mesh are built from the same coordinates, so the OpenMC tally bins line up
one-to-one with the elements of the `[Mesh]` (in bin ordering, i.e. x-fastest, then y,
then z, matching `openmc::StructuredMesh::get_bin_from_indices`).

Two mesh types are supported via the `mesh_type` parameter:

- `regular` (default): a uniform Cartesian mesh specified by `lower_left` and either
  `upper_right` or `width`, together with the number of cells `nx`/`ny`/`nz` along each axis.
  This maps to an `openmc::RegularMesh`.
- `rectilinear`: a non-uniform Cartesian mesh specified by the grid line coordinates
  `x_grid`/`y_grid`/`z_grid`. This maps to an `openmc::RectilinearMesh`.

The number of spatial directions is given by `dimensions` (1, 2, or 3). Units should match
OpenMC's length scale (centimeters); see [OpenMCCellAverageProblem.md#scaling].

### Writing results into an array auxvariable

When the tally uses a single energy filter (and no other external filters), setting the
`add_energy_array` parameter to `true` (default `false`) additionally writes each score into a
MOOSE **array auxvariable** named after the score, with one component per energy bin. This is
in addition to the usual scalar auxvariables (one per `(score, energy-bin)`), which are always
created. The array components are labeled with the energy bin names (e.g. `g1`, `g2`, ...).

For example, a `flux` score with a 2-group `EnergyFilter` produces the scalar auxvariables
`flux_g1` and `flux_g2`, plus a `flux` array auxvariable with two components. Components of the
array auxvariable can be integrated with `ElementIntegralArrayVariablePostprocessor`, or output
and inspected in Exodus.

## Example Input File Syntax

!syntax parameters /Problem/Tallies/StructuredMeshTally

!syntax inputs /Problem/Tallies/StructuredMeshTally