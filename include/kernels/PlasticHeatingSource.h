#ifndef PLASTICHEATINGSOURCE_H
#define PLASTICHEATINGSOURCE_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PlasticHeatingSource;

template <>
InputParameters validParams<PlasticHeatingSource>();

/**
 * This kernel calculates the heat source term corresponding to crack friction
 */
class PlasticHeatingSource : public HeatSource
{
public:
  PlasticHeatingSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  /// Contribution of umdamaged strain energy to damage evolution
  const MaterialProperty<Real> & _W0p;
  const MaterialProperty<Real> & _W0p_old;

  /// Variation of undamaged strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dW0p_dstrain;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // PLASTICHEATINGSOURCE_H
