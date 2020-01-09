#include <RethetaEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarElemSolverAlgorithm.h>
#include <AssembleScalarElemOpenSolverAlgorithm.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <Enums.h>
#include <FieldFunctions.h>
#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSystem.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <Realm.h>
#include <Realms.h>
#include <ScalarMassElemSuppAlgDep.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>

// template for supp algs
#include <AlgTraits.h>
#include <kernel/KernelBuilder.h>
#include <kernel/KernelBuilderLog.h>

// consolidated
#include <AssembleElemSolverAlgorithm.h>
#include <kernel/ScalarMassElemKernel.h>
#include <kernel/ScalarAdvDiffElemKernel.h>
#include <kernel/ScalarUpwAdvDiffElemKernel.h>


// edge kernels
#include <edge_kernels/ScalarEdgeSolverAlg.h>

// node kernels
#include <node_kernels/NodeKernelUtils.h>
#include <node_kernels/BLTGammaNodeKernel.h>
#include <node_kernels/ScalarGclNodeKernel.h>

// nso
#include <nso/ScalarNSOElemKernel.h>
#include <nso/ScalarNSOKeElemSuppAlg.h>
#include <nso/ScalarNSOElemSuppAlgDep.h>

#include <overset/UpdateOversetFringeAlgorithmDriver.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{


RethetaEquationSystem::RethetaEquationSystem(
    EquationSystems& eqSystems)
    : EquationSystem(eqSystems, "RethetaEQS", "retheta_transition"),
      retheta_(NULL),
      dRetdx_(NULL),
      RetTmp_(NULL),
      visc_(NULL),
      tvisc_(NULL),
      evisc_(NULL),
      assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "retheta_transition", "dRetdx"))
{
  dofName_ = "Retheta_transition";
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("retheta_transition");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_RET_TRANSITION);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // determine nodal gradient form
  set_nodal_gradient("retheta_transition");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for Retheta: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

}
RethetaEquationSystem::~RethetaEquationSystem()
{

}
void RethetaEquationSystem::register_nodal_fields(
    stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  retheta_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "retheta_transition", numStates));
  stk::mesh::put_field_on_mesh(*retheta_, *part, nullptr);
  realm_.augment_restart_variable_list("retheta_transition");

  dRetdx_ = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dRetdx"));
  stk::mesh::put_field_on_mesh(*dwdx_, *part, nDim, nullptr);

  RetTmp_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "RetTmp"));
  stk::mesh::put_field_on_mesh(*RetTmp_, *part, nullptr);

  visc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*visc_, *part, nullptr);

  tvisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
  stk::mesh::put_field_on_mesh(*tvisc_, *part, nullptr);

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_sdr"));
  stk::mesh::put_field_on_mesh(*evisc_, *part, nullptr);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &RetN = retheta_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &RetNp1 = retheta_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &RetNp1, &RetN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
  
}
void RethetaEquationSystem::register_interior_algorithm(
    stk::mesh::Part *part)
{
  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &RetNp1 = retheta_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dRetdxNone = dRetdx_->field_of_state(stk::mesh::StateNone);

  // non-solver, dwdx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = NULL;
    if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
      theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &RetNp1, &dRetdxNone);
    }
    else {
      throw std::runtime_error("Retheta equation system:: Element schemes not implemented");
    }
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // solver; interior contribution (advection + diffusion)
  if (!realm_.solutionOptions_->useConsolidatedSolverAlg_) {

    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
      = solverAlgDriver_->solverAlgMap_.find(algType);
    if (itsi == solverAlgDriver_->solverAlgMap_.end()) {
      SolverAlgorithm* theAlg = NULL;
      if (realm_.realmUsesEdges_) {
        theAlg = new ScalarEdgeSolverAlg(realm_, part, this, retheta_, dRetdx_, evisc_);
      }
      else {
        throw std::runtime_error("Retheta equation system:: Element schemes not implemented");
      }
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }

    // Check if the user has requested CMM or LMM algorithms; if so, do not
    // include Nodal Mass algorithms
    std::vector<std::string> checkAlgNames = {
      "specific_dissipation_rate_time_derivative",
      "lumped_specific_dissipation_rate_time_derivative"};
    bool elementMassAlg = supp_alg_is_requested(checkAlgNames);
    auto& solverAlgMap = solverAlgDriver_->solverAlgMap_;
    process_ngp_node_kernels(
      solverAlgMap, realm_, part, this,
      [&](AssembleNGPNodeSolverAlgorithm& nodeAlg) {
        if (!elementMassAlg)
          nodeAlg.add_kernel<ScalarMassBDFNodeKernel>(realm_.bulk_data(), retheta_);

        nodeAlg.add_kernel<BLTRe0tNodeKernel>(realm_.meta_data());
      },
      [&](AssembleNGPNodeSolverAlgorithm& nodeAlg, std::string& srcName) {
        if (srcName == "gcl") {
          nodeAlg.add_kernel<ScalarGclNodeKernel>(realm_.bulk_data(), retheta_);
          NaluEnv::self().naluOutputP0() << " - " << srcName << std::endl;
        }
        else
          throw std::runtime_error("Retheta EqSys: Invalid source term: " + srcName);
      });
  }
}

void RethetaEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}
void RethetaEquationSystem::reinitialize_linear_system()
{
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_RET_TRANSITION;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }
    // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("retheta_transition");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_RET_TRANSITION);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
  
}
void RethetaEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &RetN = retheta_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &RetNp1 = retheta_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), RetN, RetNp1, realm_.get_activate_aura());
}

} // namespace nalu
} // namespace Sierra
