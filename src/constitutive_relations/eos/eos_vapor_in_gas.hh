/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for a combination of air and vapor pressure.  Mass density is not
  available, not because it can't be calculated, but because it depends upon
  omega.  It's not really needed, and if it were, would not fit the EOS
  interface without serious work.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_VAPOR_IN_GAS_HH_
#define AMANZI_RELATIONS_EOS_VAPOR_IN_GAS_HH_

#include "Teuchos_ParameterList.hpp"
#include "eos.hh"
#include "dbc.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSVaporInGas : public EOS {

public:
  EOSVaporInGas(Teuchos::ParameterList& eos_plist);

  double MassDensity(std::vector<double>& params) { AMANZI_ASSERT(0); return 0.0; }
  double DMassDensityDT(std::vector<double>& params) { AMANZI_ASSERT(0); return 0.0; }
  double DMassDensityDp(std::vector<double>& params) { AMANZI_ASSERT(0); return 0.0; }

  double MolarDensity(std::vector<double>& params);
  double DMolarDensityDT(std::vector<double>& params);
  double DMolarDensityDp(std::vector<double>& params);

  bool IsConstantMolarMass() { return false; }
  double MolarMass() { AMANZI_ASSERT(0); return 0.0; }

protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  Teuchos::RCP<EOS> gas_eos_;

 private:
  static Utils::RegisteredFactory<EOS,EOSVaporInGas> factory_;

};

} // namespace
} // namespace

#endif
