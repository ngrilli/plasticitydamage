#include "CrystalPlasticitySlipResistanceDislo.h"

template<>
InputParameters validParams<CrystalPlasticitySlipResistanceDislo>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addParam<std::string>("uo_state_var_name", "Name of state variable property: Same as state variable user object specified in input file.");
  params.addClassDescription("Dislocation based constitutive models' slip resistance base class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticitySlipResistanceDislo::CrystalPlasticitySlipResistanceDislo(const InputParameters & parameters) :
    CrystalPlasticitySlipResistance(parameters),
    _mat_prop_state_var(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_state_var_name")))
{
}

bool
CrystalPlasticitySlipResistanceDislo::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{

  DenseVector<Real> rho(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    rho(i) = _mat_prop_state_var[qp][i];
    // val[i] is Taylor slip resistance
    val[i] = 100.0 * std::sqrt( rho(i) );
  }

  return true;
}
