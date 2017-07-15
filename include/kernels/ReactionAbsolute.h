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

#ifndef REACTIONABSOLUTE_H
#define REACTIONABSOLUTE_H

#include "Kernel.h"

// Forward Declaration
class ReactionAbsolute;

template <>
InputParameters validParams<ReactionAbsolute>();

class ReactionAbsolute : public Kernel
{
public:
  ReactionAbsolute(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};
#endif // REACTIONABSOLUTE_H
