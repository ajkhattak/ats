/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */

#ifndef MPC_SUBSURFACE_HH_
#define MPC_SUBSURFACE_HH_

#include "TreeOperator.hh"
#include "pk_physical_bdf_default.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCDelegateEWCSubsurface;

namespace Operators {
class PDE_Diffusion;
class PDE_DiffusionWithGravity;
class PDE_Accumulation;
class Operator;
class UpwindTotalFlux;
class UpwindArithmeticMean;
class Upwinding;
}

namespace Flow {
class Richards;
}

class MPCSubsurface : public StrongMPC<PK_PhysicalBDF_Default> {

 public:

  MPCSubsurface(Teuchos::ParameterList& pk_tree_list,
                const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln) :
      PK(pk_tree_list, global_list, S, soln),
      StrongMPC<PK_PhysicalBDF_Default>(pk_tree_list, global_list, S, soln),
      update_pcs_(0)
  {
    dump_ = plist_->get<bool>("dump preconditioner", false);
  }

  // -- Initialize owned (dependent) variables.
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // update the predictor to be physically consistent
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> up0,
          Teuchos::RCP<TreeVector> up);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
    UpdatePreconditioner(t, up, h, true);
  }

  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h, bool assemble);
  
  // preconditioner application
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  // virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  //     ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
  //                      Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

  Teuchos::RCP<Operators::TreeOperator> preconditioner() { return preconditioner_; }

  void ComputeDivCorrection( const Teuchos::RCP<const CompositeVector>& flux,
                             const Teuchos::RCP<const CompositeVector>& k,
                             const Teuchos::RCP<const CompositeVector>& dk,
                             const Teuchos::RCP<Epetra_MultiVector>& res);
  
 protected:

  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_BLOCK_DIAGONAL = 1,
    PRECON_PICARD = 2,
    PRECON_EWC = 3,
    PRECON_NO_FLOW_COUPLING = 4,    
  };

  Teuchos::RCP<Operators::TreeOperator> preconditioner_;
  Teuchos::RCP<Operators::TreeOperator> linsolve_preconditioner_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // preconditioner methods
  PreconditionerType precon_type_;

  // Additional precon terms
  //   equations are given by:
  // 1. conservation of WC: dWC/dt + div q = 0
  // 2. conservation of E:  dE/dt + div K grad T + div hq = 0

  // dWC / dT off-diagonal block
  Teuchos::RCP<Operators::Operator> dWC_dT_block_;
  // -- d ( div q ) / dT  terms
  Teuchos::RCP<Operators::PDE_DiffusionWithGravity> ddivq_dT_;
  Teuchos::RCP<Operators::Upwinding> upwinding_dkrdT_;
  // -- d ( dWC/dt ) / dT terms
  Teuchos::RCP<Operators::PDE_Accumulation> dWC_dT_;

  // dE / dp off-diagonal block
  Teuchos::RCP<Operators::Operator> dE_dp_block_;
  // -- d ( div K grad T ) / dp terms
  Teuchos::RCP<Operators::PDE_Diffusion> ddivKgT_dp_;
  Teuchos::RCP<Operators::Upwinding> upwinding_dKappa_dp_;
  // -- d ( div hq ) / dp terms
  Teuchos::RCP<Operators::PDE_DiffusionWithGravity> ddivhq_dp_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_hkr_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_dhkr_dp_;
  // -- d ( dE/dt ) / dp terms
  Teuchos::RCP<Operators::PDE_Accumulation> dE_dp_;

  // dE / dT on-diagonal block additional terms that use q info
  // -- d ( div hq ) / dT terms
  Teuchos::RCP<Operators::PDE_DiffusionWithGravity> ddivhq_dT_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_dhkr_dT_;

  
  // friend sub-pk Richards (need K_, some flags from private data)
  Teuchos::RCP<Flow::Richards> richards_pk_;

  Key domain_name_;
  Key temp_key_;
  Key pres_key_;
  Key e_key_;
  Key wc_key_;
  Key tc_key_;
  Key uw_tc_key_;
  Key kr_key_;
  Key uw_kr_key_;
  Key enth_key_;
  Key hkr_key_;
  Key uw_hkr_key_;
  Key energy_flux_key_;
  Key mass_flux_key_;
  Key mass_flux_dir_key_;
  Key rho_key_;
  Key ddivq_dT_key_, ddivKgT_dp_key_;

  bool is_fv_;
  
  // EWC delegate
  Teuchos::RCP<MPCDelegateEWCSubsurface> ewc_;

  // cruft for easier global debugging
  bool dump_;
  int update_pcs_;
  Teuchos::RCP<Debugger> db_;
  
private:
  // factory registration
  static RegisteredPKFactory<MPCSubsurface> reg_;

};


} // namespace

#endif


