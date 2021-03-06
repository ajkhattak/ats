/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived WeakMPC class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "pk_mpcsubcycled_ats.hh"

namespace Amanzi {

RegisteredPKFactory<PK_MPCSubcycled_ATS> PK_MPCSubcycled_ATS::reg_("pk_mpcsubcycled_ats");

} // namespace
