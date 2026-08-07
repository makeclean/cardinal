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
/*                                                                      */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                      */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#pragma once

#include "TallyBase.h"

#include "openmc/tallies/filter_mesh.h"

/**
 * A base class for tallies which score on an OpenMC mesh and write the results back into
 * the elements of the MOOSE problem mesh - i.e. the OpenMC mesh bins line up one-to-one
 * with the elements of the [Mesh] block.
 *
 * This class contains the machinery common to [MeshTally] and [StructuredMeshTally]: the
 * mesh-appropriate estimator defaults and random-ray checks in the constructor, creation
 * of the OpenMC mesh filter, removal of the OpenMC mesh on reset, and the volumetric
 * normalization and element-writing loop performed when storing tally results.
 */
class MeshTallyBase : public TallyBase
{
public:
  static InputParameters validParams();

  MeshTallyBase(const InputParameters & parameters);

  /// A function to reset the tally. Deletes the OpenMC mesh added by this object.
  virtual void resetTally() override;

protected:
  /**
   * Create the OpenMC mesh filter for the mesh this object tallies on and store it
   * in _mesh_filter.
   * @param[in] mesh_index index of the OpenMC mesh the filter should wrap
   * @return the newly-created OpenMC mesh filter
   */
  openmc::MeshFilter * createMeshFilter(unsigned int mesh_index);

  /// Index into openmc::model::meshes of the mesh added by this object.
  virtual unsigned int openmcMeshIndex() const = 0;

  /// Number of bins in the mesh.
  virtual unsigned int nBins() const = 0;

  /// Volume in cm^3 of the mesh bin.
  virtual Real binVolume(unsigned int bin) const = 0;

  /**
   * Map a mesh bin onto the element of the problem's [Mesh] which stores its value.
   * For a structured mesh this is the identity mapping; for a translated or
   * block-restricted unstructured mesh it can be an indirection.
   */
  virtual unsigned int binToElemId(unsigned int bin) const = 0;

  /// The OpenMC mesh filter for this mesh tally.
  openmc::MeshFilter * _mesh_filter = nullptr;
};