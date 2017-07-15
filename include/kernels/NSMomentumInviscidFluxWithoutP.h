/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NSMOMENTUMINVISCIDFLUXWITHOUTP_H
#define NSMOMENTUMINVISCIDFLUXWITHOUTP_H

#include "NSKernel.h"


// ForwardDeclarations
class NSMomentumInviscidFluxWithoutP;

template<>
InputParameters validParams<NSMomentumInviscidFluxWithoutP>();



/**
 * The inviscid flux (only convective) for the
 * momentum conservation equations.
 */
class NSMomentumInviscidFluxWithoutP : public NSKernel
{
public:

  NSMomentumInviscidFluxWithoutP(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  // Coupled variables (pressure is not needed)
  // const VariableValue & _pressure;

  // Parameters
  const unsigned int _component;

private:
  // To be used from both the on and off-diagonal
  // computeQpJacobian functions.  Variable numbering
  // should be in the canonical ordering regardless of
  // Moose's numbering.
  Real computeJacobianHelper(unsigned int m);
};

#endif
