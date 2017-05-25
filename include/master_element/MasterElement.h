/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MasterElement_h
#define MasterElement_h

#include <AlgTraits.h>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <array>

#ifdef __INTEL_COMPILER
#define POINTER_RESTRICT restrict
#else
#define POINTER_RESTRICT __restrict__
#endif

namespace stk {
struct topology;
}

namespace sierra{
namespace nalu{

namespace Jacobian{
enum Direction
{
  S_DIRECTION = 0,
  T_DIRECTION = 1,
  U_DIRECTION = 2
};
}

struct ElementDescription;
class MasterElement;

MasterElement*
get_surface_master_element(
  const stk::topology& theTopo,
  ElementDescription* desc = nullptr,
  std::string quadType = "GaussLegendre");

MasterElement*
get_volume_master_element(
  const stk::topology& theTopo,
  ElementDescription* desc = nullptr,
  std::string quadType = "GaussLegendre");

class MasterElement
{
public:
  MasterElement();
  virtual ~MasterElement();

  virtual void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) {
    throw std::runtime_error("determinant not implemented");}

  virtual void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv) {
    throw std::runtime_error("gij not implemented");}

  virtual void nodal_grad_op(
    const int nelem,
    double *deriv,
    double * error ) {
    throw std::runtime_error("nodal_grad_op not implemented");}

  virtual void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("face_grad_op not implemented; avoid this element type at open bcs, walls and symms");}
  
  virtual const int * adjacentNodes() {
    throw std::runtime_error("adjacentNodes not implementedunknown bc");
    return NULL;}

  virtual const int * ipNodeMap(int ordinal = 0) {
      throw std::runtime_error("ipNodeMap not implemented");
      return NULL;}

  virtual void shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shape_fcn not implemented"); }

  virtual void shifted_shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shifted_shape_fcn not implemented"); }

  virtual int opposingNodes(
    const int ordinal, const int node) {
    throw std::runtime_error("adjacentNodes not implemented"); }

  virtual int opposingFace(
    const int ordinal, const int node) {
    throw std::runtime_error("opposingFace not implemented"); 
    return 0; }

  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord) {
    throw std::runtime_error("isInElement not implemented"); 
    return 1.0e6; }

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result) {
    throw std::runtime_error("interpolatePoint not implemented"); }
  
  virtual void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc) {
    throw std::runtime_error("general_shape_fcn not implement"); }

  virtual void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("general_face_grad_op not implemented");}

  virtual void general_normal(
    const double *isoParCoord,
    const double *coords,
    double *normal) {
    throw std::runtime_error("general_normal not implemented");}

  virtual void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords) {
    throw std::runtime_error("sidePcoords_to_elemPcoords");}

  virtual const int* side_node_ordinals(int sideOrdinal) {
    throw std::runtime_error("side_node_ordinals not implemented");
  }

  double isoparametric_mapping(const double b, const double a, const double xi) const;
  bool within_tolerance(const double & val, const double & tol);
  double vector_norm_sq(const double * vect, int len);

  int nDim_;
  int nodesPerElement_;
  int numIntPoints_;
  double scaleToStandardIsoFac_;

  std::vector<int> lrscv_;
  std::vector<int> ipNodeMap_;
  std::vector<int> oppNode_;
  std::vector<int> oppFace_;
  std::vector<double> intgLoc_;
  std::vector<double> intgLocShift_;
  std::vector<double> intgExpFace_;
  std::vector<double> nodeLoc_;
  std::vector<int> sideNodeOrdinals_;
  std::vector<int> sideOffset_;

  // FEM
  std::vector<double>weights_;
};

// Hex 8 subcontrol volume
class HexSCV : public MasterElement
{
public:

  HexSCV();
  virtual ~HexSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

};

// Hex 8 subcontrol surface
class HexSCS : public MasterElement
{
public:

  HexSCS();
  virtual ~HexSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);
  
  const int* side_node_ordinals(int sideOrdinal) final;

  double parametric_distance(const std::vector<double> &x);
};

class HexahedralP2Element : public MasterElement
{
public:
  using Traits = AlgTraitsHex27;

  HexahedralP2Element();
  virtual ~HexahedralP2Element() {}

  void shape_fcn(double *shpfc);
  void shifted_shape_fcn(double *shpfc);

protected:
  struct ContourData {
    Jacobian::Direction direction;
    double weight;
  };

  int tensor_product_node_map(int i, int j, int k) const;

  double gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double shifted_gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double tensor_product_weight(
    int s1Node, int s2Node, int s3Node,
    int s1Ip, int s2Ip, int s3Ip) const;

  double tensor_product_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const;

  virtual void eval_shape_functions_at_ips();
  virtual void eval_shape_functions_at_shifted_ips();

  virtual void eval_shape_derivs_at_ips();
  virtual void eval_shape_derivs_at_shifted_ips();

  void eval_shape_derivs_at_face_ips();

  void set_quadrature_rule();
  void GLLGLL_quadrature_weights();

  void hex27_shape_deriv(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  double parametric_distance(const std::array<double, 3>& x);

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  const double scsDist_;
  const int nodes1D_;
  const int numQuad_;

  // quadrature info
  std::vector<double> gaussAbscissae1D_;
  std::vector<double> gaussAbscissae_;
  std::vector<double> gaussAbscissaeShift_;
  std::vector<double> gaussWeight_;
  std::vector<double> scsEndLoc_;

  std::vector<int> stkNodeMap_;

  std::vector<double> shapeFunctions_;
  std::vector<double> shapeFunctionsShift_;
  std::vector<double> shapeDerivs_;
  std::vector<double> shapeDerivsShift_;
  std::vector<double> expFaceShapeDerivs_;

private:
  void hex27_shape_fcn(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;
};

// 3D Quad 27 subcontrol volume
class Hex27SCV : public HexahedralP2Element
{
public:
  Hex27SCV();
  virtual ~Hex27SCV() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

private:
  void set_interior_info();

  double jacobian_determinant(
    const double *POINTER_RESTRICT elemNodalCoords,
    const double *POINTER_RESTRICT shapeDerivs ) const;

  std::vector<double> ipWeight_;
};

// 3D Hex 27 subcontrol surface
class Hex27SCS : public HexahedralP2Element
{
public:
  Hex27SCS();
  virtual ~Hex27SCS() {}

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  const int * adjacentNodes();

  const int * ipNodeMap(int ordinal = 0);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  const int* side_node_ordinals(int sideOrdinal) final;
private:
  void set_interior_info();
  void set_boundary_info();

  template <Jacobian::Direction dir>
  void area_vector(const double *POINTER_RESTRICT elemNodalCoords,
    double *POINTER_RESTRICT shapeDeriv,
    double *POINTER_RESTRICT areaVector ) const;

  void gradient(
    const double *POINTER_RESTRICT elemNodalCoords,
    const double *POINTER_RESTRICT shapeDeriv,
    double *POINTER_RESTRICT grad,
    double *POINTER_RESTRICT det_j ) const;

  std::vector<ContourData> ipInfo_;
  int ipsPerFace_;
};

// Tet 4 subcontrol volume
class TetSCV : public MasterElement
{
public:

  TetSCV();
  virtual ~TetSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tet_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// Tet 4 subcontrol surface
class TetSCS : public MasterElement
{
public:

  TetSCS();
  virtual ~TetSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tet_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );
  
  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // helper
  double parametric_distance(const std::vector<double> &x);

  const int* side_node_ordinals(int sideOrdinal) final;

};

// Pyramid 5 subcontrol volume
class PyrSCV : public MasterElement
{
public:

  PyrSCV();
  virtual ~PyrSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// Pyramid 5 subcontrol surface
class PyrSCS : public MasterElement
{
public:

  PyrSCS();
  virtual ~PyrSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void pyr_derivative(
    const int npts,
    const double *intLoc,
    double *deriv);

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double *error);

  const int* side_node_ordinals(int sideOrdinal) final;

  double parametric_distance(const std::array<double,3>& x);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
};

// Wedge 6 subcontrol volume
class WedSCV : public MasterElement
{
public:
  WedSCV();
  virtual ~WedSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  void wedge_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// Wedge 6 subcontrol surface
class WedSCS : public MasterElement
{
public:
  WedSCS();
  virtual ~WedSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void wedge_derivative(
    const int npts,
    const double *intLoc,
    double *deriv);

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);
  
  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
  
  void wedge_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // helper functions to isInElement
  double parametric_distance( const double X, const double Y);
  double parametric_distance( const std::vector<double> &x);

  const int* side_node_ordinals(int sideOrdinal) final;

};

// 2D Quad 4 subcontrol volume
class Quad2DSCV : public MasterElement
{
public:
  Quad2DSCV();
  virtual ~Quad2DSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// 2D Quad 4 subcontrol surface
class Quad2DSCS : public MasterElement
{
public:
  Quad2DSCS();
  virtual ~Quad2DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
     const double *coords,
     double *gupperij,
     double *gij,
     double *deriv);

  const int * adjacentNodes();

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
  
  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  const int* side_node_ordinals(int sideOrdinal) final;
};

class QuadrilateralP2Element : public MasterElement
{
public:
  using Traits = AlgTraitsQuad9_2D;

  QuadrilateralP2Element();
  virtual ~QuadrilateralP2Element() {}

  void shape_fcn(double *shpfc);
  void shifted_shape_fcn(double *shpfc);
protected:
  struct ContourData {
    Jacobian::Direction direction;
    double weight;
  };

  void set_quadrature_rule();
  void GLLGLL_quadrature_weights();

  int tensor_product_node_map(int i, int j) const;

  double gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double shifted_gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double tensor_product_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const;

  double tensor_product_weight(int s1Node, int s1Ip) const;

  double parametric_distance(const std::array<double, 2>& x);

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void eval_shape_functions_at_ips();
  void eval_shape_functions_at_shifted_ips();

  void eval_shape_derivs_at_ips();
  void eval_shape_derivs_at_shifted_ips();

  void eval_shape_derivs_at_face_ips();

  const double scsDist_;
  const int nodes1D_;
  int numQuad_;

  //quadrature info
  std::vector<double> gaussAbscissae1D_;
  std::vector<double> gaussAbscissae_;
  std::vector<double> gaussAbscissaeShift_;
  std::vector<double> gaussWeight_;

  std::vector<int> stkNodeMap_;
  std::vector<double> scsEndLoc_;

  std::vector<double> shapeFunctions_;
  std::vector<double> shapeFunctionsShift_;
  std::vector<double> shapeDerivs_;
  std::vector<double> shapeDerivsShift_;
  std::vector<double> expFaceShapeDerivs_;
private:
  void quad9_shape_fcn(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  void quad9_shape_deriv(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;
};

// 3D Quad 27 subcontrol volume
class Quad92DSCV : public QuadrilateralP2Element
{
public:
  Quad92DSCV();
  virtual ~Quad92DSCV() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

private:
  void set_interior_info();

  double jacobian_determinant(
    const double *POINTER_RESTRICT elemNodalCoords,
    const double *POINTER_RESTRICT shapeDerivs ) const;

  std::vector<double> ipWeight_;
};

// 3D Hex 27 subcontrol surface
class Quad92DSCS : public QuadrilateralP2Element
{
public:
  Quad92DSCS();
  virtual ~Quad92DSCS() {}

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  const int * ipNodeMap(int ordinal = 0);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  const int* side_node_ordinals(int sideOrdinal) final;


private:
  void set_interior_info();
  void set_boundary_info();

  template <Jacobian::Direction direction> void
  area_vector(
    const double *POINTER_RESTRICT elemNodalCoords,
    double *POINTER_RESTRICT shapeDeriv,
    double *POINTER_RESTRICT areaVector ) const;

  std::vector<ContourData> ipInfo_;
  int ipsPerFace_;
};

// 2D Tri 3 subcontrol volume
class Tri2DSCV : public MasterElement
{
public:
  Tri2DSCV();
  virtual ~Tri2DSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  void tri_shape_fcn(
    const int &npts,
    const double *par_coord,
    double* shape_fcn);

};

// 2D Tri 3 subcontrol surface
class Tri2DSCS : public MasterElement
{
public:
  Tri2DSCS();
  virtual ~Tri2DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tri_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);
  
  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  double tri_parametric_distance(
    const std::vector<double> &x);
  
  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  const int* side_node_ordinals(int sideOrdinal) final;


};

// 3D Quad 4
class Quad3DSCS : public MasterElement
{
public:

  Quad3DSCS();
  virtual ~Quad3DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_normal(
    const double *isoParCoord,
    const double *coords,
    double *normal);

  void non_unit_face_normal(
    const double * par_coord,
    const double * elem_nodal_coor,
    double * normal_vector );
  
  double parametric_distance(const std::vector<double> &x);

  const double elemThickness_;
};

// 3D Quad 9
class Quad93DSCS : public HexahedralP2Element
{
public:
  Quad93DSCS();
  virtual ~Quad93DSCS() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_normal(
    const double *isoParCoord,
    const double *coords,
    double *normal);

private:
  void set_interior_info();
  void eval_shape_functions_at_ips() final;
  void eval_shape_derivs_at_ips() final;

  void eval_shape_functions_at_shifted_ips() final;
  void eval_shape_derivs_at_shifted_ips() final;

  void area_vector(
    const double *POINTER_RESTRICT coords,
    const double *POINTER_RESTRICT shapeDerivs,
    double *POINTER_RESTRICT areaVector) const;

  void quad9_shape_fcn(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  void quad9_shape_deriv(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  void non_unit_face_normal(
    const double *isoParCoord,
    const double *elemNodalCoord,
    double *normalVector);

  double parametric_distance(const std::vector<double> &x);

  std::vector<double> ipWeight_;
  const int surfaceDimension_;
};

// 3D Tri 3
class Tri3DSCS : public MasterElement
{
public:

  Tri3DSCS();
  virtual ~Tri3DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
     double *shpfc);

   void shifted_shape_fcn(
     double *shpfc);

   void tri_shape_fcn(
     const int &npts,
     const double *par_coord,
     double* shape_fcn);

   double isInElement(
     const double *elemNodalCoord,
     const double *pointCoord,
     double *isoParCoord);

  double parametric_distance(
    const std::vector<double> &x);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_normal(
    const double *isoParCoord,
    const double *coords,
    double *normal);
};

// edge 2d
class Edge2DSCS : public MasterElement
{
public:
  Edge2DSCS();
  virtual ~Edge2DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_normal(
    const double *isoParCoord,
    const double *coords,
    double *normal);

  double parametric_distance(const std::vector<double> &x);

  const double elemThickness_;  
};

// edge 2d
class Edge32DSCS : public QuadrilateralP2Element
{
public:
  Edge32DSCS();
  virtual ~Edge32DSCS() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

private:
  void area_vector(
    const double *POINTER_RESTRICT coords,
    const double s,
    double *POINTER_RESTRICT areaVector) const;

  std::vector<double> ipWeight_;
};

} // namespace nalu
} // namespace Sierra

#endif
