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

#include "damageIC.h"

template<>
InputParameters validParams<damageIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("coefficient", "The value of the initial condition");
  params.addRequiredParam<Real>("xmax", "Max X coordinate");
  params.addRequiredParam<Real>("ymax", "Max Y coordinate");
  params.addRequiredParam<Real>("sizeimagex", "X size of the image");
  params.addRequiredParam<Real>("sizeimagey", "Y size of the image");
  return params;
}

damageIC::damageIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _coefficient(getParam<Real>("coefficient")),
    _xmax(getParam<Real>("xmax")),
    _ymax(getParam<Real>("ymax")),
    _sizeimagex(getParam<Real>("sizeimagex")),
    _sizeimagey(getParam<Real>("sizeimagey"))
{}

Real
damageIC::value(const Point & p)
{
  /**
   * _value * x
   * The Point class is defined in libMesh.  The spatial
   * coordinates x,y,z can be accessed individually using
   * the parenthesis operator and a numeric index from 0..2
   */

  Real xPixel;
  Real yPixel;

  xPixel = std::floor( p(0) * _sizeimagex / _xmax ) + 1;
  xPixel = std::max(1.0,xPixel);
  xPixel = std::min(_sizeimagex,xPixel);
  yPixel = std::floor( p(1) * _sizeimagey / _ymax ) + 1;
  yPixel = std::max(1.0,yPixel);
  yPixel = std::min(_sizeimagey,yPixel);

  FILE * damagefile;
  //FILE * outfile;

  damagefile = fopen("damageonecolumn.txt","r");
  //outfile = fopen("out.txt","a");

  int j, k;
  int damagetemp;

  j = xPixel + (yPixel - 1) * _sizeimagex; // row of the wanted damage value
  for ( k=0; k<j; k++)
  {
    std::fscanf( damagefile, "%d", &damagetemp);
  }

  //fclose(outfile);
  fclose(damagefile);

  //if (damagetemp > 0.0) {
  //    std::cout << damagetemp << std::endl;
  //}

  return damagetemp * _coefficient;
}
