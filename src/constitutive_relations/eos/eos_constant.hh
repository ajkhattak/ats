/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Simple EOS for constant density.
  Defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_CONSTANT_HH_
#define AMANZI_RELATIONS_EOS_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSConstant : public EOSConstantMolarMass {

public:
  explicit EOSConstant(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(std::vector<double>& params) override { return rho_; }
  virtual double DMassDensityDT(std::vector<double>& params) override { return 0.0; }
  virtual double DMassDensityDp(std::vector<double>& params) override { return 0.0; }

  virtual double DMolarDensityDT(std::vector<double>& params) override { return 0.0; }
  virtual double DMolarDensityDp(std::vector<double>& params) override { return 0.0; }

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;

  static Utils::RegisteredFactory<EOS,EOSConstant> factory_;

};

} // namespace
} // namespace

#endif
