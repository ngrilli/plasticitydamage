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

#ifndef DAMAGEIC_H
#define DAMAGEIC_H

// MOOSE Includes
#include "InitialCondition.h"

// Forward Declarations
class damageIC;

template<>
InputParameters validParams<damageIC>();

/**
 * ExampleIC just returns a constant value.
 */
class damageIC : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  damageIC(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p) override;

private:
  Real _coefficient;
  Real _xmax, _ymax, _sizeimagex, _sizeimagey;
};

#endif //DAMAGEIC_H
