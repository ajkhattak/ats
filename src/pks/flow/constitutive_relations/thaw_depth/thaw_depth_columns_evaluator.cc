/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes the thaw depth 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "thaw_depth_columns_evaluator.hh"

namespace Amanzi {
namespace Flow {



ThawDepthColumnsEvaluator::ThawDepthColumnsEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  Key dset_name = plist.get<std::string>("domain set name", "column");
  
  domain_ = Keys::getDomain(my_key_); //surface_column domain
  int col_id = Keys::getDomainSetIndex<int>(domain_);
  
  Key domain_ss = Keys::getDomainInSet(dset_name, col_id);
  
  temp_key_ = Keys::readKey(plist, domain_ss, "temperature", "temperature");
  dependencies_.insert(temp_key_);

  trans_width_ =  plist_.get<double>("transition width [K]", 0.2);
}
  

ThawDepthColumnsEvaluator::ThawDepthColumnsEvaluator(const ThawDepthColumnsEvaluator& other)
  : SecondaryVariableFieldEvaluator(other),
    temp_key_(other.temp_key_)
{}
  
Teuchos::RCP<FieldEvaluator>
ThawDepthColumnsEvaluator::Clone() const
{
  return Teuchos::rcp(new ThawDepthColumnsEvaluator(*this));
}


void
ThawDepthColumnsEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  double trans_temp = 273.15 + 0.5*trans_width_;

  // search through the column and find the deepest unfrozen cell

  std::string domain_ss = Keys::getDomain(temp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);
  AmanziGeometry::Point z_centroid(top_z_centroid);

  const auto& temp_c = *S->GetFieldData(temp_key_)
    ->ViewComponent("cell", false);
    
  int col_cells = temp_c.MyLength();
  for (int i=0; i!=col_cells; ++i) {
    if (temp_c[0][i] >= trans_temp) {
      z_centroid = S->GetMesh(domain_ss)->face_centroid(i+1);
    }
  }
  
  res_c[0][0] = top_z_centroid[2] - z_centroid[2];
  
}
  
void
ThawDepthColumnsEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
ThawDepthColumnsEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ThawDepthColumnsEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{

  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
