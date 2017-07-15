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

#ifndef DAMAGEICXYZ_H
#define DAMAGEICXYZ_H

// MOOSE Includes
#include "InitialCondition.h"

// Forward Declarations
class damageICxyz;

template<>
InputParameters validParams<damageICxyz>();

/**
 * ExampleIC just returns a constant value.
 */
class damageICxyz : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  damageICxyz(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p) override;

private:
  Real _coefficient;
  Real _radius;
  unsigned int _xyzfilesize;

};

#endif //DAMAGEICXYZ_H
