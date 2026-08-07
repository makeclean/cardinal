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

#ifdef ENABLE_OPENMC_COUPLING
#include "MeshTallyBase.h"

#include "OpenMCCellAverageProblem.h"

InputParameters
MeshTallyBase::validParams()
{
  auto params = TallyBase::validParams();
  params.addClassDescription(
      "Base class for mesh tallies, which score on an OpenMC mesh and write the results "
      "back into the elements of the problem mesh.");
  return params;
}

MeshTallyBase::MeshTallyBase(const InputParameters & parameters) : TallyBase(parameters)
{
  // The random ray solver requires tracklength estimators, which are not supported for
  // mesh tallies because the regions are OpenMC mesh cells rather than flat source regions.
  if (_openmc_problem.runRandomRay())
    mooseError("Mesh tallies are not supported when using the random ray solver!");

  const bool nu_scatter =
      std::find(_tally_score.begin(), _tally_score.end(), "nu-scatter") != _tally_score.end();

  // Mesh tallies don't support tracklength estimators.
  if (isParamValid("estimator"))
  {
    if (_estimator == openmc::TallyEstimator::TRACKLENGTH)
      paramError("estimator",
                 "Tracklength estimators are currently incompatible with mesh tallies!");
  }
  else
    _estimator = nu_scatter ? openmc::TallyEstimator::ANALOG : openmc::TallyEstimator::COLLISION;
}

void
MeshTallyBase::resetTally()
{
  TallyBase::resetTally();

  // Erase the OpenMC mesh.
  openmc::model::meshes.erase(openmc::model::meshes.begin() + openmcMeshIndex());
}

openmc::MeshFilter *
MeshTallyBase::createMeshFilter(unsigned int mesh_index)
{
  _mesh_filter = dynamic_cast<openmc::MeshFilter *>(openmc::Filter::create("mesh"));
  _mesh_filter->set_mesh(mesh_index);
  return _mesh_filter;
}

Real
MeshTallyBase::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 const std::vector<OMCTensor> & tally_vals,
                                 bool norm_by_src_rate)
{
  Real total = 0.0;

  for (unsigned int ext_bin = 0; ext_bin < _num_ext_filter_bins; ++ext_bin)
  {
    for (unsigned int e = 0; e < nBins(); ++e)
    {
      Real unnormalized_tally = tally_vals[local_score](ext_bin * nBins() + e);

      // Divide each tally by the volume it corresponds to in MOOSE because we will
      // apply it as a volumetric tally (per unit volume). Because the OpenMC mesh is
      // constructed in units of cm, we need to adjust by the scaling factor.
      Real volumetric_tally = unnormalized_tally;
      volumetric_tally *= norm_by_src_rate
                              ? _openmc_problem.tallyMultiplier(_tally_score[local_score],
                                                                _local_mean_tally[local_score]) /
                                    binVolume(e) * _openmc_problem.scaling() *
                                    _openmc_problem.scaling() * _openmc_problem.scaling()
                              : 1.0;
      total += _ext_bins_to_skip[ext_bin] ? 0.0 : unnormalized_tally;

      auto var = var_numbers[local_score * _num_ext_filter_bins + ext_bin];
      fillElementalAuxVariable(var, {binToElemId(e)}, volumetric_tally);
    }
  }

  return total;
}
#endif