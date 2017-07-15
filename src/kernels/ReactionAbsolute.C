/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ReactionAbsolute.h"

template <>
InputParameters
validParams<ReactionAbsolute>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

ReactionAbsolute::ReactionAbsolute(const InputParameters & parameters) : Kernel(parameters) {}

Real
ReactionAbsolute::computeQpResidual()
{
  if ( _u[_qp] < 0.0 ) {
    return _test[_i][_qp] * _u[_qp];
  }
  else {
    return 0.0;
  }
}

Real
ReactionAbsolute::computeQpJacobian()
{
  if ( _u[_qp] < 0.0 ) {
    return _test[_i][_qp] * _phi[_j][_qp];
  }
  else {
    return 0.0;
  }
}
