#ifndef CRACKPROPAGATIONHEATSOURCE_H
#define CRACKPROPAGATIONHEATSOURCE_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class CrackPropagationHeatSource;

template <>
InputParameters validParams<CrackPropagationHeatSource>();

/**
 * This kernel calculates the heat source term corresponding to crack friction
 */
class CrackPropagationHeatSource : public HeatSource
{
public:
  CrackPropagationHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  const VariableValue & _c;
  const VariableValue & _c_old;

  const unsigned int _c_var;

  /// Contribution of umdamaged strain energy to damage evolution
  const MaterialProperty<Real> & _G0_pos;

  /// Variation of undamaged strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dG0_pos_dstrain;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // CRACKPROPAGATIONHEATSOURCE_H
