/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTICE                        */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#ifdef ENABLE_OPENMC_COUPLING
#include "StructuredMeshTally.h"

#include "OpenMCCellAverageProblem.h"
#include "MooseUtils.h"

registerMooseObject("CardinalApp", StructuredMeshTally);

InputParameters
StructuredMeshTally::validParams()
{
  auto params = MeshTallyBase::validParams();
  params.addClassDescription(
      "A tally which maps results onto an OpenMC structured mesh (regular or rectilinear) "
      "and writes them into the elements of the [Mesh] block.");

  params.addParam<MooseEnum>("mesh_type",
                             MooseEnum("regular rectilinear", "regular"),
                             "Whether to use a regular (uniform) or rectilinear "
                             "(non-uniform) OpenMC mesh.");

  params.addRangeCheckedParam<unsigned int>(
      "dimensions", "dimensions >= 1 & dimensions <= 3",
      "Number of spatial dimensions of the structured mesh (1, 2, or 3).");

  params.addRequiredParam<std::vector<Real>>("lower_left",
                                             "Lower-left corner coordinates of the mesh.");

  params.addParam<std::vector<Real>>("upper_right",
                                     "Upper-right corner coordinates of the mesh. For a regular "
                                     "mesh, either 'upper_right' or 'width' must be provided.");

  params.addParam<std::vector<Real>>(
      "width", "Uniform width of the mesh cells along each axis. For a regular mesh, "
               "either 'upper_right' or 'width' must be provided.");

  params.addParam<unsigned int>("nx", "Number of cells in the x-direction.");
  params.addParam<unsigned int>("ny", "Number of cells in the y-direction.");
  params.addParam<unsigned int>("nz", "Number of cells in the z-direction.");

  params.addParam<std::vector<Real>>(
      "x_grid", "Sorted x-coordinates of the grid lines for a rectilinear mesh.");
  params.addParam<std::vector<Real>>(
      "y_grid", "Sorted y-coordinates of the grid lines for a rectilinear mesh.");
  params.addParam<std::vector<Real>>(
      "z_grid", "Sorted z-coordinates of the grid lines for a rectilinear mesh.");

  return params;
}

StructuredMeshTally::StructuredMeshTally(const InputParameters & parameters)
  : MeshTallyBase(parameters),
    _dimension(getParam<unsigned int>("dimensions")),
    _mesh_type(getParam<MooseEnum>("mesh_type") == "rectilinear" ? structured_mesh::MeshType::RECTILINEAR
                                                                  : structured_mesh::MeshType::REGULAR)
{
  // The structured mesh is fully independent of the [Mesh], so we reuse it as-is.
  _mesh = std::make_unique<structured_mesh::StructuredMesh>(
      buildCoordinates(), _dimension, _mesh_type);
}

std::array<std::vector<Real>, 3>
StructuredMeshTally::buildCoordinates() const
{
  std::array<std::vector<Real>, 3> coords;

  if (_mesh_type == structured_mesh::MeshType::RECTILINEAR)
  {
    // The grid lines fully describe the mesh; cell counts are derived from them.
    const std::vector<std::vector<Real>> grids{getParam<std::vector<Real>>("x_grid"),
                                               getParam<std::vector<Real>>("y_grid"),
                                               getParam<std::vector<Real>>("z_grid")};
    for (unsigned int i = 0; i < 3; ++i)
      if (i < _dimension)
      {
        if (grids[i].size() < 2)
        {
          std::string names[]{"x_grid", "y_grid", "z_grid"};
          paramError(names[i], "Each used grid must be provided with at least two points.");
        }
        coords[i] = grids[i];
      }
      else
        coords[i] = {};
  }
  else
  {
    const auto & lower_left = getParam<std::vector<Real>>("lower_left");
    if (lower_left.size() != _dimension)
      paramError("lower_left",
                 "The 'lower_left' must have a number of entries equal to 'dimensions'.");

    std::vector<Real> upper_right;
    std::vector<unsigned int> n;
    const std::vector<unsigned int> counts = {isParamValid("nx") ? getParam<unsigned int>("nx") : 0,
                                              isParamValid("ny") ? getParam<unsigned int>("ny") : 0,
                                              isParamValid("nz") ? getParam<unsigned int>("nz") : 0};

    if (isParamValid("upper_right"))
      upper_right = getParam<std::vector<Real>>("upper_right");
    else if (isParamValid("width"))
    {
      const auto & width = getParam<std::vector<Real>>("width");
      if (width.size() != _dimension)
        paramError("width", "The 'width' must have a number of entries equal to 'dimensions'.");

      upper_right = lower_left;
      for (unsigned int i = 0; i < _dimension; ++i)
        upper_right[i] = lower_left[i] + counts[i] * width[i];
    }
    else
      paramError("upper_right", "Either 'upper_right' or 'width' must be provided.");

    if (upper_right.size() != _dimension)
      paramError("upper_right",
                 "The 'upper_right' must have a number of entries equal to 'dimensions'.");

    for (unsigned int a = 0; a < _dimension; ++a)
    {
      unsigned int ncells = counts[a];
      if (ncells == 0)
        paramError("nx", "A cell count must be provided for each used axis.");

      coords[a].resize(ncells + 1);
      for (unsigned int j = 0; j <= ncells; ++j)
        coords[a][j] = lower_left[a] + (upper_right[a] - lower_left[a]) * j / ncells;
    }
    for (unsigned int a = _dimension; a < 3; ++a)
      coords[a].clear();
  }

  return coords;
}

std::pair<unsigned int, openmc::Filter *>
StructuredMeshTally::spatialFilter()
{
  // Build the OpenMC structured mesh and the native libMesh mesh.
  _mesh->buildOpenMCMesh();
  _mesh->setID(-1);
  _mesh->buildLibMeshMesh(_openmc_problem.getMooseMesh().getMesh().comm());

  // Ensure that the [Mesh] block is the same grid so the results can be copied into
  // auxiliary variables on the problem mesh.
  checkMeshGridMatchesProblemMesh();

  createMeshFilter(_mesh->meshIndex());

  return std::make_pair(openmc::model::tally_filters.size() - 1, _mesh_filter);
}

void
StructuredMeshTally::checkMeshGridMatchesProblemMesh()
{
  // The results are written onto the elements of the problem's [Mesh] block, so that
  // mesh must be the same replicated grid as the structured mesh (same bounds, cell
  // counts, and x-fastest element ordering), analogous to the requirement placed on the
  // 'mesh_template' for an unstructured MeshTally.
  const auto & moose_mesh = _openmc_problem.getMooseMesh().getMesh();
  if (!moose_mesh.is_replicated())
    mooseError("Structured mesh tallies require a replicated mesh in the [Mesh] block!");

  for (unsigned int e = 0; e < _mesh->nBins(); ++e)
  {
    auto elem_ptr = _openmc_problem.getMooseMesh().queryElemPtr(e);
    if (!elem_ptr)
      continue;

    const auto pt = _mesh->binCentroid(e);
    Point centroid_bin = {pt(0), pt(1), pt(2)};

    // The OpenMC mesh is in units of cm; scale the [Mesh] before comparing.
    Point centroid_mesh = elem_ptr->vertex_average() * _openmc_problem.scaling();

    for (unsigned int j = 0; j < 3; ++j)
      if (!MooseUtils::absoluteFuzzyEqual(centroid_mesh(j), centroid_bin(j)))
        paramError("mesh_type",
                   "Centroid for bin " + Moose::stringify(e) + " in the structured mesh (cm): " +
                       _openmc_problem.printPoint(centroid_bin) +
                       "\ndoes not match centroid for element " + Moose::stringify(e) +
                       " in the [Mesh] (cm): " + _openmc_problem.printPoint(centroid_mesh) +
                       "!\n\nThe [Mesh] block must be identical to the structured mesh grid, "
                       "with matching cell counts and element ordering (x-fastest).");
  }
}
#endif
