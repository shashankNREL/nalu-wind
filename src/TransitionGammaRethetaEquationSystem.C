#include <TransitionGammaRethetaEquationSystem.h>
#include <ShearStressTransportEquationSystem.h>
#include <AlgorithmDriver.h>
#include <ComputeSSTMaxLengthScaleElemAlgorithm.h>
#include <FieldFunctions.h>
#include <master_element/MasterElement.h>
#include <master_element/MasterElementFactory.h>
#include <NaluEnv.h>
#include <SpecificDissipationRateEquationSystem.h>
#include <SolutionOptions.h>
#include <TurbKineticEnergyEquationSystem.h>
#include <Realm.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// basic c++
#include <cmath>
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TransitionModelEquationSystem - manage Gamma-Retheta
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TransitionModelEquationSystem::TransitionModelEquationSystem(
    EquationSystems& eqSystems):
    EquationSystem(eqSystems, "TransitionWrap"),
    ShearStressTransportEquationSystem(eqSystems),
    GammaEqSys_(NULL),
    ReThetaEqSys_(NULL),
    SSTEqSys_(NULL),
    gamma_(NULL),
    reTheta_(NULL)
{
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create SST, Gamma and ReTheta equation
  SSTEqSys_ = new ShearStressTransportEquationSystem(eqSystems);
  GammaEqSys_ = new GammaEquationSystem(eqSystems);
  ReThetaEqSys_ = new ReThetaEquationSystem(eqSystems);
}


TransitionModelEquationSystem::~TransitionModelEquationSystem()
{

}
//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void TransitionModelEquationSystem::initialize(){
  SSTEqSys_->initialize();
  GammaEqSys_->convergenceTolerance_ = convergenceTolerance_;
  ReThetaEqSys_->convergenceTolerance_ = convergenceTolerance_;
}
//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void TransitionModelEquationSystem::register_nodal_fields(stk::mesh::Part *part)
{
  stk::mesh::MetaData &metaData = realm_.meta_data();
  const int numStates = realm_.number_of_states();
  SSTEqSys_->register_nodal_fields(part);

  // register Gamma and Retheta nodal fields
  gamma_ = &(metaData.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "intermittency", numStates));
  stk::mesh::put_field_on_mesh(*gamma_, *part, nullptr);
  reTheta_ = &(metaData.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "Retheta", numStates));
  stk::mesh::put_field_on_mesh(*reTheta_, *part, nullptr);

}
//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void TransitionModelEquationSystem::register_interior_algorithm(stk::mesh::Part *part)
{
  SSTEqSys_->register_interior_algorithm(part);

}
void TransitionModelEquationSystem::solve_and_update_external()
{

  if (SSTEqSys_->isInit_){
    // compute projected nodal gradients
    SSTEqSys_->tkeEqSys_->compute_projected_nodal_gradient();
    SSTEqSys_->sdrEqSys_->assemble_nodal_gradient();
    SSTEqSys_->clip_min_distance_to_wall();
    
    // deal with DES option
    if ( SST_DES == realm_.solutionOptions_->turbulenceModel_ )
      SSTEqSys_->sstMaxLengthScaleAlgDriver_->execute();

    SSTEqSys_->isInit_ = false;
  }else if (realm_.has_mesh_motion()) {
    if (realm_.currentNonlinearIteration_ == 1)
      SSTEqSys_->clip_min_distance_to_wall();

    if (SST_DES == realm_.solutionOptions_->turbulenceModel_)
      SSTEqSys_->sstMaxLengthScaleAlgDriver_->execute();
  }

  // compute blending for SST model
  SSTEqSys_->compute_f_one_blending();

  // SST effective viscosity for k and omega
  SSTEqSys_->tkeEqSys_->compute_effective_diff_flux_coeff();
  SSTEqSys_->sdrEqSys_->compute_effective_diff_flux_coeff();

  // wall values
  SSTEqSys_->tkeEqSys_->compute_wall_model_parameters();
  SSTEqSys_->sdrEqSys_->compute_wall_model_parameters();
  
  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // tke and sdr assemble, load_complete and solve; Jacobi iteration
    GammaEqSys_->assemble_and_solve(GammaEqSys_->gTmp_);
    ReThetaEqSys_->assemble_and_solve(ReThetaEqSys_->reTmp_);
    SSTEqSys_->tkeEqSys_->assemble_and_solve(SSTEqSys_->tkeEqSys_->kTmp_);
    SSTEqSys_->sdrEqSys_->assemble_and_solve(SSTEqSys_->sdrEqSys_->wTmp_);
    
    // update each
    SSTEqSys_->update_and_clip();

    // compute projected nodal gradients
    SSTEqSys_->tkeEqSys_->compute_projected_nodal_gradient();
    SSTEqSys_->sdrEqSys_->assemble_nodal_gradient();
  }



}
