#ifndef CRYSTALPLASTICITYSLIPRESISTANCEDISLO_H
#define CRYSTALPLASTICITYSLIPRESISTANCEDISLO_H

#include "CrystalPlasticitySlipResistance.h"

class CrystalPlasticitySlipResistanceDislo;

template<>
InputParameters validParams<CrystalPlasticitySlipResistanceDislo>();

/**
 * Dislocation based constitutive model slip resistance userobject class.
 */
class CrystalPlasticitySlipResistanceDislo : public CrystalPlasticitySlipResistance
{
 public:
  CrystalPlasticitySlipResistanceDislo(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:
  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;
};

#endif // CRYSTALPLASTICITYSLIPRESISTANCEDISLO_H
