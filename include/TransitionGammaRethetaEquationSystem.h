/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TransitionGammaRethetaEquationSystem_h
#define TransitionGammaRethetaEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class EquationSystems;
class AlgorithmDriver;
class ShearStressTransportEquationSystem; 
class RethetaEquationSystem;
class GammaEquationSystem; 

class TransitionModelEquationSystem : public EquationSystem {

public:

  TransitionModelEquationSystem(
    EquationSystems& equationSystems);
  virtual ~TransitionModelEquationSystem();
  
  virtual void initialize();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  virtual void register_interior_algorithm(stk::mesh::Part *part);

  virtual void solve_and_update();

  void initial_work();
  void post_adapt_work();

  ShearStressTransportEquationSystem *SSTEqSys_;
  GammaEquationSystem *GammaEqSys_;
  ReThetaEquationSystem *ReThetaEqSys_;

  ScalarFieldType *gamma_;
  ScalarFieldType *reTheta_;

};

}  // namespace nalu
}  // namespace Sierra

#endif
