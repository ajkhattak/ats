/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Exchange flux between multiple continua.

/*!

Evaluates the following exchange flux model:

.. math::
   q_{exchange} = k_r K \frac{\Gamma}{\delta} (p_M - p_m)

where :math:`p` is the pressure of the Macro and micro porespaces,
respectively, K is some measure of an absolute permeability, :math:`\Gamma [-]`
is the exchange coefficient, :math:`\delta [m]` is a unit of distance
characterizing the typical distance between pore, and :math:`k_r` is the
relative permeability, which is upwinded based on the larger of the two
pressures.

Note that the expected domain for this is the micropore domain, but may be
changed on the input line.

.. _micropore-macropore-flux-evaluator-spec:
.. admonition:: micropore-macropore-flux-evaluator-spec

   * `"micropore domain`" ``[string]`` **DOMAIN** Defaults to the domain of the flux's
     variable name.

   * `"macropore domain`" ``[string]`` **macropore**

   * `"micropore macropore flux model parameters`" ``[micropore-macropore-flux-model-spec]``

   KEYS:
   - `"micropore pressure`" **pressure**
   - `"macropore pressure`" **MACROPORE_DOMAIN-pressure**
   - `"micropore relative permeability`" **relative_permeability**
   - `"macropore relative permeability`" **MACROPORE_DOMAIN-relative_permeability**
   - `"permeability`" **permeability**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MicroporeMacroporeFluxModel;

class MicroporeMacroporeFluxEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  MicroporeMacroporeFluxEvaluator(Teuchos::ParameterList& plist);
  MicroporeMacroporeFluxEvaluator(const MicroporeMacroporeFluxEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<MicroporeMacroporeFluxModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key pm_key_;
  Key pM_key_;
  Key krM_key_;
  Key krm_key_;
  Key K_key_;

  Teuchos::RCP<MicroporeMacroporeFluxModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,MicroporeMacroporeFluxEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace


