#ifndef PLASTICHEATINGSOURCEMIEHE2016_H
#define PLASTICHEATINGSOURCEMIEHE2016_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PlasticHeatingSourceMiehe2016;

template <>
InputParameters validParams<PlasticHeatingSourceMiehe2016>();

/**
 * This kernel calculates the heat source term corresponding to plastic deformation
 * as in Miehe 2016: mechanical term: Cauchy stress * plastic strain_rate
 * using finite strain formalism
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

  /// Cauchy stress
  const MaterialProperty<RankTwoTensor> & _stress;

  /// deformation gradient
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;

  /// plastic deformation gradient
  const MaterialProperty<RankTwoTensor> & _fp;
  const MaterialProperty<RankTwoTensor> & _fp_old;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // PLASTICHEATINGSOURCEMIEHE2016_H
