/*
  The ideal gas equation of state evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_GENERAL_EOS_IDEAL_GAS_EVALUATOR_HH_
#define AMANZI_GENERAL_EOS_IDEAL_GAS_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace General {
namespace Relations {

class EosIdealGasModel;

class EosIdealGasEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  EosIdealGasEvaluator(Teuchos::ParameterList& plist);
  EosIdealGasEvaluator(const EosIdealGasEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<EosIdealGasModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Key pres_key_;

  Teuchos::RCP<EosIdealGasModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EosIdealGasEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
