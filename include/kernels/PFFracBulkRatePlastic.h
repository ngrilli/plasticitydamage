/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFFRACBULKRATEPLASTIC_H
#define PFFRACBULKRATEPLASTIC_H
/**
 * Phase field based fracture model
 * This kernel computes the residual and jacobian for bulk free energy contribution to c
 * Refer to Formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation 63
 */
#include "KernelValue.h"
#include "RankTwoTensor.h"

//Forward Declarations
class PFFracBulkRatePlastic;

template<>
InputParameters validParams<PFFracBulkRatePlastic>();

class PFFracBulkRatePlastic : public KernelValue
{
public:

  PFFracBulkRatePlastic(const InputParameters & parameters);

protected:

  enum PFFunctionType
  {
    Jacobian,
    Residual
  };

  virtual Real precomputeQpResidual();
  virtual Real precomputeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  ///Critical energy release rate for fracture
  const MaterialProperty<Real> & _gc_prop;
  ///Contribution of umdamaged strain energy (elastic and plastic) to damage evolution
  const MaterialProperty<Real> & _W0e;
  const MaterialProperty<Real> & _W0p;
  ///Variation of undamaged strain energy driving damage evolution with strain (elastic and plastic)
  const MaterialProperty<RankTwoTensor> * _dW0e_dstrain;
  const MaterialProperty<RankTwoTensor> * _dW0p_dstrain;
  ///Auxiliary variable: beta = Laplacian of c
  const VariableValue & _betaval;
  const unsigned int _beta_var;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  std::string _base_name;

  ///Characteristic length, controls damage zone thickness
  Real _l;
  ///Viscosity parameter ( visco -> 0, rate independent )
  Real _visco;
  /// Threshold fracture energy
  Real _Wc;

 private:

};
#endif //PFFRACBULKRATEPLASTIC_H
