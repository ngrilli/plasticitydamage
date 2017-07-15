/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCrystalPlasticityW0p.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityW0p>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity base class: FCC system with power law flow rule implemented. Damage. Miehe 2016");
  params.addRequiredCoupledVar("c","Order parameter for damage");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");

  return params;
}

FiniteStrainCrystalPlasticityW0p::FiniteStrainCrystalPlasticityW0p(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _W0e(declareProperty<Real>("W0e")), // elastic energy
    _W0p(declareProperty<Real>("W0p")), // plastic energy
    _W0p_old(declarePropertyOld<Real>("W0p")), // plastic energy of previous increment
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dW0e_dstrain(declareProperty<RankTwoTensor>("dW0e_dstrain")),
    _dW0p_dstrain(declareProperty<RankTwoTensor>("dW0p_dstrain")),
    _etens(LIBMESH_DIM),
    _eigval(LIBMESH_DIM)
{
}

void
FiniteStrainCrystalPlasticityW0p::preSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _W0p_tmp = _W0p_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _W0p_tmp = _W0p_tmp_old = _W0p_old[_qp];
      _accslip_tmp_old = _acc_slip_old[_qp];
    }
    else
    {
      _gss_tmp = _gss_tmp_old;
      _W0p_tmp = _W0p_tmp_old;
    }
  }
}

void
FiniteStrainCrystalPlasticityW0p::solveStatevar()
{
  Real gmax, gdiff;
  unsigned int iterg;
  std::vector<Real> gss_prev(_nss);

  gmax = 1.1 * _gtol;
  iterg = 0;

  while (gmax > _gtol && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;

    update_energies();

    postSolveStress(); // Update _fp_old_inv = _fp_old

    gss_prev = _gss_tmp;

    update_slip_system_resistance(); // Update slip system resistance

    gmax = 0.0;
    for (unsigned i = 0; i < _nss; ++i)
    {
      gdiff = std::abs(gss_prev[i] - _gss_tmp[i]); // Calculate increment size

      if (gdiff > gmax)
        gmax = gdiff;
    }
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Hardness Integration error gmax", gmax, "\n");
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCrystalPlasticityW0p::postSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss[_qp] = _gss_tmp;
    _W0p[_qp] = _W0p_tmp;
    _acc_slip[_qp] = _accslip_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _W0p[_qp] = _W0p_tmp;
      _acc_slip[_qp] = _accslip_tmp;
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _W0p_tmp_old = _W0p_tmp;
      _accslip_tmp_old = _accslip_tmp;
    }
  }
}

// Update slip system resistance and elastic and plastic work
void
FiniteStrainCrystalPlasticityW0p::update_energies()
{
  RankTwoTensor cauchy_stress, WpToTrace, invFe, ee, ce, iden;
  Real detFe, detFe_old;
  Real c = _c[_qp];

  if (_max_substep_iter == 1) //No substepping
  {
    _W0p_tmp = _W0p_old[_qp];
  }
  else
  {
    _W0p_tmp = _W0p_tmp_old;
  }

  // Update elastic and plastic work
  detFe = _fe.det();
  //_fe_old[_qp] = _deformation_gradient_old[_qp] * _fp_old_inv; // works with 1 substep onl
  //detFe_old = _fe_old[_qp].det();
  invFe = _fe.inverse();

  // _pk2[_qp] is the updated piola-kirchoff
  cauchy_stress = _fe * _pk2[_qp] * _fe.transpose()/detFe;
  //_cauchy_out[_qp] = _deformation_gradient_old[_qp];
  //WeToTrace = cauchy_stress * ( _fe - _fe_old[_qp] ) * invFe * detFe;
  WpToTrace = cauchy_stress * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;

  iden.zero();
  iden.addIa(1.0);

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  _W0e[_qp] = 0.5 * _pk2[_qp].doubleContraction(ee);
  //_W0e_tmp += WeToTrace.trace();
  _W0p_tmp += WpToTrace.trace();
  _dW0e_dstrain[_qp] = cauchy_stress;
  _dW0p_dstrain[_qp] = cauchy_stress;
  _dstress_dc[_qp] = -cauchy_stress * (2.0 * (1.0 - c));

}

void
FiniteStrainCrystalPlasticityW0p::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  Real trD;
  Real c = _c[_qp];
  Real xfac = Utility::pow<2>(1.0-c) + _kdamage;

  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();
  ee = 0.5 * ( ce - iden );

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr(i);

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  //ce.symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

  //for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //{
  //  _eigval[i] = ( _eigval[i] - 1.0 ) / 2.0;
  //}

  // Tensors of outerproduct of eigen vectors
  //for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //    for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
  //      _etens[i](j, k) = _eigvec(j, i) * _eigvec(k, i);

  //RankTwoTensor ee0pos, ee0neg;
  //ee0pos.zero();
  //ee0neg.zero();
  //for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //{
  //  ee0pos += _etens[i] * ( std::abs(_eigval[i]) + _eigval[i] ) / 2.0;
  //  ee0neg += _etens[i] * ( std::abs(_eigval[i]) - _eigval[i] ) / 2.0;
  //}

  pk2_new = _elasticity_tensor[_qp] * ee;

  //pk2_new.addIa(-1.0/3.0 * pk2_new.trace());
  //pk2_new.addIa( (_Bulk_Modulus_Ref/_n_Murnaghan) * ( 1.0 - std::pow( 1.0/_fe.det() , _n_Murnaghan ) ) );

  //trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  //trD /= _deformation_gradient_old[_qp].det();

  //pk2_new.addIa( _C0 * trD * trD );
  //pk2_new.addIa( _C1 * abs(trD) );

  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau. Override to modify.
void
FiniteStrainCrystalPlasticityW0p::getSlipIncrements()
{
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i)) *
                    copysign(1.0, _tau(i)) * _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _slip_incr(i) = _slip_incr_tol * _dt * copysign(1.0, _tau(i));
      //_err_tol = true;
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
#endif
      return;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dslipdtau(i) = _a0(i) / _xm(i) *
                    std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i) - 1.0) / _gss_tmp[i] *
                    _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _dslipdtau(i) = 0.0;
    }
  }
}
