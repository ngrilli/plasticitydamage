#ifndef THERMALEXPANSIONHEATSOURCE_H
#define THERMALEXPANSIONHEATSOURCE_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class ThermalExpansionHeatSource;

template <>
InputParameters validParams<ThermalExpansionHeatSource>();

/**
 * This kernel calculates the heat source term corresponding to crack friction
 */
class ThermalExpansionHeatSource : public HeatSource
{
public:
  ThermalExpansionHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  const Real _thermal_expansion_heat;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // THERMALEXPANSIONHEATSOURCE_H
