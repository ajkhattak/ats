/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Fully three-phase (air, water, ice) permafrost energy equation, with only
mobile water.

Inherits TwoPhase instead of EnergyBase to pick up the enthalpy from TwoPhase.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_INTERFROST_ENERGY_HH_
#define PKS_ENERGY_INTERFROST_ENERGY_HH_

#include "energy_three_phase.hh"

namespace Amanzi {
namespace Energy {

class InterfrostEnergy : public ThreePhase {

public:

  InterfrostEnergy(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    ThreePhase(FElist, plist, S, solution) {}

  // Virtual destructor
  virtual ~InterfrostEnergy() {}
  
  // -- accumulation term
  virtual void UpdatePreconditioner(double t,
          Teuchos::RCP<const TreeVector> up, double h);

 protected:
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& f);
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

private:
  static RegisteredPKFactory<InterfrostEnergy> reg_;

  friend  class MPCCoupledFlowEnergy;
};

} // namespace Energy
} // namespace Amanzi

#endif
