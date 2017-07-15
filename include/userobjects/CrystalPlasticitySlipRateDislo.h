#ifndef CRYSTALPLASTICITYSLIPRATEDISLO_H
#define CRYSTALPLASTICITYSLIPRATEDISLO_H

#include "CrystalPlasticitySlipRate.h"
#include "RankTwoTensor.h"

class CrystalPlasticitySlipRateDislo;

template<>
InputParameters validParams<CrystalPlasticitySlipRateDislo>();

/**
 * Dislocation based constitutive model slip rate userobject class.
 */
class CrystalPlasticitySlipRateDislo : public CrystalPlasticitySlipRate
{
 public:
  CrystalPlasticitySlipRateDislo(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real /*dt*/, std::vector<Real> & val) const;
  virtual void calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const;

 protected:
  virtual void readFileFlowRateParams();
  virtual void getFlowRateParams();

  const MaterialProperty<std::vector<Real> > & _mat_prop_state_var;

  const MaterialProperty<RankTwoTensor> & _pk2;

  DenseVector<Real> _a0;
  DenseVector<Real> _xm;

  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction;
};

#endif // CRYSTALPLASTICITYSLIPRATEDISLO_H
