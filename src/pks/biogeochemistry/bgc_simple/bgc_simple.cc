/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Chonggang Xu

   Simple implementation of CLM's Century model for carbon decomposition and a
   simplified 2-PFT (sedge, moss) vegetation model for creating carbon.

   CURRENT ASSUMPTIONS:
     1. parallel decomp not in the vertical
     2. fields are not ordered along the column, and so must be copied
     3. all columns have the same number of cells
   ------------------------------------------------------------------------- */

#include "MeshPartition.hh"

#include "bgc_simple_funcs.hh"

#include "bgc_simple.hh"

namespace Amanzi {
namespace BGC {

BGCSimple::BGCSimple(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                     Teuchos::ParameterList& FElist,
                     const Teuchos::RCP<TreeVector>& solution) :
    PKPhysicalBase(plist, FElist, solution),
    PKDefaultBase(plist, FElist, solution),
    ncells_per_col_(-1) {}

// is a PK
// -- Setup data
void BGCSimple::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);

  // my mesh is the subsurface mesh, but we need the surface mesh, index by column, as well
  surf_mesh_ = S->GetMesh("surface");
  soil_part_name_ = plist_->get<std::string>("soil partition name");

  // Create the additional, non-managed data structures
  int nPools = plist_->get<int>("number of carbon pools", 7);

  // -- SoilCarbonParameters
  Teuchos::ParameterList& sc_params = plist_->sublist("soil carbon parameters");
  std::string mesh_part_name = sc_params.get<std::string>("mesh partition");
  const Functions::MeshPartition& mp = *S->GetMeshPartition(mesh_part_name);
  const std::vector<std::string>& regions = mp.regions();

  for (std::vector<std::string>::const_iterator region=regions.begin();
       region!=regions.end(); ++region) {
    sc_params_.push_back(Teuchos::rcp(
        new SoilCarbonParameters(nPools, sc_params.sublist(*region))));
  }

  // -- PFTs -- old and new!
  Teuchos::ParameterList& pft_params = plist_->sublist("pft parameters");
  std::vector<std::string> pft_names;
  for (Teuchos::ParameterList::ConstIterator lcv=pft_params.begin();
       lcv!=pft_params.end(); ++lcv) {
    std::string pft_name = lcv->first;
    pft_names.push_back(pft_name);
  }

  int ncols = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  pfts_old_.resize(ncols);
  pfts_.resize(ncols);
  for (unsigned int col=0; col!=ncols; ++col) {
    int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, col);
    ColIterator col_iter(*mesh_, f);
    std::size_t ncol_cells = col_iter.size();

    // unclear which this should be:
    // -- col area is the true face area
    double col_area = mesh_->face_area(f);
    // -- col area is the projected face area
    // double col_area = surf_mesh_->cell_volume(col);

    if (ncells_per_col_ < 0) {
      ncells_per_col_ = ncol_cells;
    } else {
      ASSERT(ncol_cells == ncells_per_col_);
    }

    pfts_old_[col].resize(pft_names.size());
    pfts_[col].resize(pft_names.size());

    for (int i=0; i!=pft_names.size(); ++i) {
      std::string pft_name = pft_names[i];
      Teuchos::ParameterList& pft_plist = pft_params.sublist(pft_name);
      pfts_old_[col][i] = Teuchos::rcp(new PFT(pft_name, ncol_cells));
      pfts_old_[col][i]->Init(pft_plist,col_area);
      pfts_[col][i] = Teuchos::rcp(new PFT(*pfts_old_[col][i]));
    }
  }

  // -- soil carbon pools
  soil_carbon_pools_.resize(ncols);
  for (unsigned int col=0; col!=ncols; ++col) {
    soil_carbon_pools_[col].resize(ncells_per_col_);
    ColIterator col_iter(*mesh_, surf_mesh_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);

    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      // col_iter[i] = cell id, mp[cell_id] = index into partition list, sc_params_[index] = correct params
      soil_carbon_pools_[col][i] = Teuchos::rcp(new SoilCarbon(sc_params_[mp[col_iter[i]]]));
    }
  }

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, nPools);

  // requirements: diagnostic, no evaluator
  S->RequireField("total_biomass", name_)->SetMesh(surf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, pft_names.size());

  // requirement: total decomp (diagnostic)
  S->RequireField("co2_decomposition", name_)->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // requirement: temp of each cell
  S->RequireFieldEvaluator("temperature");
  S->RequireField("temperature")->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirement: pressure
  S->RequireFieldEvaluator("pressure");
  S->RequireField("pressure")->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: surface cell volume
  S->RequireField("surface_cell_volume")->SetMesh(surf_mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_cell_volume");

  // requirements: Met data
  S->RequireFieldEvaluator("incoming_shortwave_radiation");
  S->RequireField("incoming_shortwave_radiation")->SetMesh(surf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("air_temperature");
  S->RequireField("air_temperature")->SetMesh(surf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("relative_humidity");
  S->RequireField("relative_humidity")->SetMesh(surf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("wind_speed");
  S->RequireField("wind_speed")->SetMesh(surf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("co2_concentration");
  S->RequireField("co2_concentration")->SetMesh(surf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // parameters
  lat_ = plist_->get<double>("latitude [degrees]", 60.);
  wind_speed_ref_ht_ = plist_->get<double>("wind speed reference height [m]", 2.0);
  cryoturbation_coef_ = plist_->get<double>("cryoturbation mixing coefficient [cm^2/yr]", 5.0);
  cryoturbation_coef_ /= 365.25e4; // convert to m^2/day

}

// -- Initialize owned (dependent) variables.
void BGCSimple::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S);

  // diagnostic variables
  S->GetFieldData("co2_decomposition", name_)->PutScalar(0.);
  S->GetField("co2_decomposition", name_)->set_initialized();
  S->GetFieldData("total_biomass", name_)->PutScalar(0.);
  S->GetField("total_biomass", name_)->set_initialized();

  // init root carbon
  Teuchos::RCP<Epetra_SerialDenseVector> col_temp =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> col_depth =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> col_dz =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  S->GetFieldEvaluator("temperature")->HasFieldChanged(S, name_);
  const Epetra_Vector& temp = *(*S->GetFieldData("temperature")->ViewComponent("cell",false))(0);

  int ncols = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int col=0; col!=ncols; ++col) {
    FieldToColumn_(col, temp, col_temp.ptr());
    ColDepthDz_(col, col_depth.ptr(), col_dz.ptr());

    int npft = pfts_old_[col].size();
    for (int i=0; i!=npft; ++i) {
      pfts_[col][i]->InitRoots(*col_temp, *col_depth, *col_dz);
      *pfts_old_[col][i] = *pfts_[col][i];
    }
  }
}

// -- Commit any secondary (dependent) variables.
void BGCSimple::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // Copy the PFT over, which includes all additional state required, commit
  // the step as succesful.
  int ncols = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int col=0; col!=ncols; ++col) {
    int npft = pfts_old_[col].size();
    for (int i=0; i!=npft; ++i) {
      *pfts_old_[col][i] = *pfts_[col][i];
    }
  }
}



// -- advance the model
bool BGCSimple::advance(double dt) {
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  // Copy the PFT from old to new, in case we failed the previous attempt at
  // this timestep.  This is hackery to get around the fact that PFTs are not
  // (but should be) in state.
  AmanziMesh::Entity_ID ncols = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID col=0; col!=ncols; ++col) {
    int npft = pfts_old_[col].size();
    for (int i=0; i!=npft; ++i) {
      *pfts_[col][i] = *pfts_old_[col][i];
    }
  }

  // grab the required fields
  Epetra_MultiVector& sc_pools = *S_next_->GetFieldData(key_, name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& co2_decomp = *S_next_->GetFieldData("co2_decomposition", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& biomass = *S_next_->GetFieldData("total_biomass", name_)
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& temp = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("pressure")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& pres = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

  S_inter_->GetFieldEvaluator("incoming_shortwave_radiation")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& qSWin = *S_inter_->GetFieldData("incoming_shortwave_radiation")
      ->ViewComponent("cell",false);

  S_inter_->GetFieldEvaluator("air_temperature")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& air_temp = *S_inter_->GetFieldData("air_temperature")
      ->ViewComponent("cell",false);

  S_inter_->GetFieldEvaluator("relative_humidity")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& rel_hum = *S_inter_->GetFieldData("relative_humidity")
      ->ViewComponent("cell",false);

  S_inter_->GetFieldEvaluator("wind_speed")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& wind_speed = *S_inter_->GetFieldData("wind_speed")
      ->ViewComponent("cell",false);

  S_inter_->GetFieldEvaluator("co2_concentration")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& co2 = *S_inter_->GetFieldData("co2_concentration")
      ->ViewComponent("cell",false);

  // note that this is used as the column area, which is maybe not always
  // right.  Likely correct for soil carbon calculations and incorrect for
  // surface vegetation calculations (where the subsurface's face area is more
  // correct?)
  S_inter_->GetFieldEvaluator("surface_cell_volume")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& scv = *S_inter_->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell", false);

  // Create workspace arrays (these should be removed when data is correctly oriented).
  Teuchos::RCP<Epetra_SerialDenseVector> temp_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> pres_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> dz_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> depth_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  // Create a workspace array for the result
  Epetra_SerialDenseVector co2_decomp_c(ncells_per_col_);

  // Grab the mesh partition to get soil properties
  Teuchos::RCP<const Functions::MeshPartition> mp = S_next_->GetMeshPartition(soil_part_name_);

  // loop over columns and apply the model
  for (AmanziMesh::Entity_ID col=0; col!=ncols; ++col) {
    // update the various soil arrays
    FieldToColumn_(col, *temp(0), temp_c.ptr());
    FieldToColumn_(col, *pres(0), pres_c.ptr());
    ColDepthDz_(col, depth_c.ptr(), dz_c.ptr());

    // copy over the soil carbon arrays
    ColIterator col_iter(*mesh_, surf_mesh_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);
    // -- serious cache thrash... --etc
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      AmanziGeometry::Point centroid = mesh_->cell_centroid(col_iter[i]);
      //      std::cout << "Col iter col=" << col << ", index i=" << i << ", cell=" << col_iter[i] << " at " << centroid << std::endl;
      for (int p=0; p!=soil_carbon_pools_[col][i]->nPools; ++p) {
        soil_carbon_pools_[col][i]->SOM[p] = sc_pools[p][col_iter[i]];
      }
    }

    // Create the Met data struct
    MetData met;
    met.qSWin = qSWin[0][col];
    met.tair = air_temp[0][col];
    met.windv = wind_speed[0][col];
    met.wind_ref_ht = wind_speed_ref_ht_;
    met.relhum = rel_hum[0][col];
    met.CO2a = co2[0][col];
    met.lat = lat_;

    // call the model
    BGCAdvance(S_inter_->time(), dt, scv[0][col], cryoturbation_coef_, met,
               *temp_c, *pres_c, *depth_c, *dz_c,
               pfts_[col], soil_carbon_pools_[col], co2_decomp_c);

    // copy back
    // -- serious cache thrash... --etc
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      for (int p=0; p!=soil_carbon_pools_[col][i]->nPools; ++p) {
        sc_pools[p][col_iter[i]] = soil_carbon_pools_[col][i]->SOM[p];
      }

      // and integrate the decomp
      co2_decomp[0][col_iter[i]] += co2_decomp_c[i];
    }

    for (int f=0; f!=pfts_[col].size(); ++f) {
      biomass[f][col] = pfts_[col][f]->totalBiomass;
    }

  } // end loop over columns

  return false;
}


// helper function for pushing field to column
void BGCSimple::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
        Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy) {
  if (col_vec == Teuchos::null) {
    col_vec = Teuchos::ptr(new Epetra_SerialDenseVector(ncells_per_col_));
  }

  ColIterator col_iter(*mesh_, surf_mesh_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    (*col_vec)[i] = vec[col_iter[i]];
  }
}

// helper function for collecting column dz and depth
void BGCSimple::ColDepthDz_(AmanziMesh::Entity_ID col,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
  AmanziMesh::Entity_ID f_above = surf_mesh_->entity_get_parent(AmanziMesh::CELL, col);
  ColIterator col_iter(*mesh_, f_above, ncells_per_col_);

  AmanziGeometry::Point surf_centroid = mesh_->face_centroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->cell_centroid(col_iter[i])[2];

    // dz
    // -- find face_below
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->face_normal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
    ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }

}




} // namespace
} // namespace