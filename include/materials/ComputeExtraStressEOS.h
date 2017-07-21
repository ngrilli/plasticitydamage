/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEEXTRASTRESSEOS_H
#define COMPUTEEXTRASTRESSEOS_H

#include "ComputeExtraStressBase.h"

/**
 * ComputeExtraStressEOS computes a volumetric extra stress including a Birch-Murnaghan EOS
 * and bulk viscosity damping for shock propagation
 * that is substituted to the volumetric stress
 * calculated by the constitutive model.
 */
class ComputeExtraStressEOS : public ComputeExtraStressBase
{
public:
  ComputeExtraStressEOS(const InputParameters & parameters);

protected:
  virtual void computeQpExtraStress();

  // exponent in Birch-Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus in Birch-Murnaghan EOS
  const Real _Bulk_Modulus_Ref;

  // Von Neumann damping coefficient
  const Real _C0;

  // Landshoff damping coefficient
  const Real _C1;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> & _stress_old;

  RankTwoTensor _extra_stress_tensor;
};

#endif // COMPUTEEXTRASTRESSEOS_H
