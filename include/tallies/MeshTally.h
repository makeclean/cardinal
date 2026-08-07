/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
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
#include "OpenMCCellAverageProblem.h"

namespace libMesh
{
class ReplicatedMesh;
}

class MeshTally : public MeshTallyBase
{
public:
  static InputParameters validParams();

  MeshTally(const InputParameters & parameters);

  /**
   * A function to generate the mesh filter needed by this object.
   * @return a pair where the first entry is the filter index in the global filter array and the
   * second entry is an OpenMC unstructured mesh filter
   */
  virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override;

  /**
   * A function which gathers the sums and means from all tallies linked to this tally. MeshTally
   * overrides this function to gather global tallies for distributed mesh tallies.
   */
  virtual void gatherLinkedSum() override;

  /**
   * A function to return if this object is adding a global tally. MeshTally modifies this behavior
   * to add a single global tally for distributed mesh tallies (which then communicate with
   * tally linkages).
   */
  virtual bool addingGlobalTally() const override { return _needs_global_tally && _instance == 0; }

protected:
  /**
   * Check the setup of the mesh template and translations. Because a simple copy transfer
   * is used to write a mesh tally onto the [Mesh], we require that the
   * meshes are identical - both in terms of the element ordering and the actual dimensions of
   * each element. This function performs as many checks as possible to ensure that the meshes
   * are indeed identical.
   */
  void checkMeshTemplateAndTranslations();

  /**
   * Mesh template file to use for creating mesh tallies in OpenMC; currently, this mesh
   * must be identical to the mesh used in the [Mesh] block because a simple copy transfer
   * is used to extract the tallies and put on the application's mesh in preparation for
   * a transfer to another MOOSE app. If not set, this indicates that tallying will be
   * performed directly on the [Mesh].
   * TODO: allow the mesh to not be identical, both in terms of using different units
   * and more general differences like not having a particular phase present
   */
  const std::string * _mesh_template_filename = nullptr;

  /// The translation to apply to the mesh template.
  Point _mesh_translation;

  /// The index into an array of mesh translations.
  const unsigned int _instance;

  /// The index of the mesh added by this tally.
  unsigned int _mesh_index;

  /// OpenMC unstructured mesh instance for use with mesh tallies
  openmc::UnstructuredMesh * _mesh_template;

  /// Whether we're using an indirection layer to map between the OpenMC mesh tally and the MOOSE mesh.
  const bool _use_dof_map;

  /**
   * For use with block restriction only. A copy of the mesh is made which only contains elements in
   * the blocks the user wishes to tally on. This is necessary at the moment as the point locators
   * used in OpenMC to find collision sites are not passed a set of block IDs to filter elements.
   * TODO: Fix this in OpenMC
   */
  std::unique_ptr<libMesh::ReplicatedMesh> _libmesh_mesh_copy;
  /// A mapping between the OpenMC bins (active block restricted elements) and all elements.
  std::vector<unsigned int> _bin_to_element_mapping;

  /// OpenMC mesh index of the mesh added by this tally.
  unsigned int openmcMeshIndex() const override { return _mesh_index; }

  /// Number of bins in the mesh.
  unsigned int nBins() const override { return _mesh_filter->n_bins(); }

  /// Volume in cm^3 of the mesh bin.
  Real binVolume(unsigned int bin) const override;

  /// Map a mesh bin onto the element of the [Mesh] which stores its value.
  unsigned int binToElemId(unsigned int bin) const override;
};
