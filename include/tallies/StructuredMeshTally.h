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

#pragma once

#include "MeshTallyBase.h"
#include "StructuredMesh.h"

#include "openmc/tallies/filter_mesh.h"

/**
 * A tally which maps results onto an OpenMC structured mesh (regular or rectilinear).
 *
 * Unlike the unstructured MeshTally, the OpenMC structured mesh is generated "on the fly"
 * from text-based parameters, and is independent of the mesh in the [Mesh] block. The
 * tally results are read back and written into MONOMIAL MOOSE auxiliary variables
 * on the [Mesh] block (which must correspond to the same grid).
 */
class StructuredMeshTally : public MeshTallyBase
{
public:
  static InputParameters validParams();

  StructuredMeshTally(const InputParameters & parameters);

  virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override;

  /// Get a reference to the structured mesh (both OpenMC and libMesh representations)
  structured_mesh::StructuredMesh & structuredMesh() { return *_mesh; }

  /// Get a reference to the structured mesh (both OpenMC and libMesh representations)
  const structured_mesh::StructuredMesh & structuredMesh() const { return *_mesh; }

protected:
  /**
   * Check that the problem's [Mesh] block is identical to the structured mesh grid
   * (same bounds, cell counts, and x-fastest element ordering) so that tally bin e maps
   * one-to-one onto element e when results are copied into the MOOSE auxiliary variables.
   */
  void checkMeshGridMatchesProblemMesh();

  /// OpenMC mesh index of the mesh added by this tally.
  unsigned int openmcMeshIndex() const override { return _mesh->meshIndex(); }

  /// Number of bins in the mesh.
  unsigned int nBins() const override { return _mesh->nBins(); }

  /// Volume in cm^3 of the mesh bin.
  Real binVolume(unsigned int bin) const override { return _mesh->binVolume(bin); }

  /// Map a mesh bin onto the element of the [Mesh] which stores its value (identity).
  unsigned int binToElemId(unsigned int bin) const override { return bin; }

  /// Build the per-axis node coordinates from the input parameters.
  std::array<std::vector<Real>, 3> buildCoordinates() const;

private:
  /// The structured mesh (both the OpenMC and libMesh representations).
  std::unique_ptr<structured_mesh::StructuredMesh> _mesh;

  /// Number of spatial dimensions.
  const unsigned int _dimension;

  /// Mesh type (regular or rectilinear).
  const structured_mesh::MeshType _mesh_type;
};