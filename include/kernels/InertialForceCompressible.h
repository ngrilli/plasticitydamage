/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INERTIALFORCECOMPRESSIBLE_H
#define INERTIALFORCECOMPRESSIBLE_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class InertialForceCompressible;

template<>
InputParameters validParams<InertialForceCompressible>();

class InertialForceCompressible : public Kernel
{
public:

  InertialForceCompressible(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const VariableValue & _density; // not anymore a material property
  const VariableValue & _u_old;
  const VariableValue & _vel_old;
  const VariableValue & _accel_old;
  const Real _beta;
  const Real _gamma;
  const Real _eta;
  const Real _alpha;

};

#endif //INERTIALFORCECOMPRESSIBLE_H
