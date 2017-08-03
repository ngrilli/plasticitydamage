#ifndef PLASTICHEATINGSOURCEMIEHE2016_H
#define PLASTICHEATINGSOURCEMIEHE2016_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PlasticHeatingSourceMiehe2016;

template <>
InputParameters validParams<PlasticHeatingSourceMiehe2016>();

/**
 * This kernel calculates the heat source term corresponding to plastic deformationa
 * as in Miehe 2016
 */
class PlasticHeatingSourceMiehe2016 : public HeatSource
{
public:
  PlasticHeatingSourceMiehe2016(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  /// Contribution of plastic strain energy to damage evolution
  const MaterialProperty<Real> & _W0p;
  const MaterialProperty<Real> & _W0p_old;

  /// Variation of plastic strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dW0p_dstrain;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  /// couples damage phase field variable
  const VariableValue & _cval;
  const unsigned int _c_var;

};

#endif // PLASTICHEATINGSOURCEMIEHE2016_H
