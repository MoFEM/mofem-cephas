/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Description: FIXME
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __COMPLEX_FOR_LAZY_NEUMANM_FORCES_HPP
#define __COMPLEX_FOR_LAZY_NEUMANM_FORCES_HPP

namespace ObosleteUsersModules {

struct NeummanForcesSurfaceComplexForLazy {

  struct MyTriangleSpatialFE;
  struct AuxMethodSpatial: public TriElementForcesAndSurcesCore::UserDataOperator {

    MyTriangleSpatialFE *myPtr;
    AuxMethodSpatial(const string &field_name,MyTriangleSpatialFE *_myPtr);
    PetscErrorCode doWork(int side, EntityType type, DataForcesAndSurcesCore::EntData &data);

  };

  struct AuxMethodMaterial: public TriElementForcesAndSurcesCore::UserDataOperator {

    MyTriangleSpatialFE *myPtr;
    AuxMethodMaterial(const string &field_name,MyTriangleSpatialFE *_myPtr);
    PetscErrorCode doWork(int side, EntityType type, DataForcesAndSurcesCore::EntData &data);

  };

  struct MyTriangleSpatialFE: public TriElementForcesAndSurcesCore {

    double *sCaleLhs;
    double *sCaleRhs;
    enum FORCES { CONSERVATIVE = 1, NONCONSERVATIVE = 2};
    FORCES typeOfForces;
    const double eps;
    bool uSeF;

    Mat Aij;
    Vec F;

    MyTriangleSpatialFE(FieldInterface &_mField,Mat _Aij,Vec &_F,double *scale_lhs,double *scale_rhs);

    int getRule(int order) { return max(1,order); };

    double *N;
    double *N_face;
    double *N_edge[3];
    double *diffN;
    double *diffN_face;
    double *diffN_edge[3];

    int order_face;
    int order_edge[3];
    double *dofs_x;
    double *dofs_x_edge[3];
    double *dofs_x_face;
    double *idofs_x;
    double *idofs_x_edge[3];
    double *idofs_x_face;
    int *dofs_x_indices;
    int *dofs_x_edge_indices[3];
    int *dofs_x_face_indices;

    int order_face_material;
    int order_edge_material[3];
    double *dofs_X;
    double *dofs_X_edge[3];
    double *dofs_X_face;
    double *idofs_X;
    double *idofs_X_edge[3];
    double *idofs_X_face;

    int *dofs_X_indices;


    ublas::vector<double> tLoc,tGlob;
    ublas::matrix<double> tLocNodal,tGlobNodal;
    double *t_loc;

    ublas::vector<int> dOfs_x_indices,dOfs_x_face_indices;
    ublas::vector<ublas::vector<int> > dOfs_x_edge_indices;
    ublas::vector<int> dOfs_X_indices,dOfs_X_face_indices;
    ublas::vector<ublas::vector<int> > dOfs_X_edge_indices;

    ublas::vector<double> dOfs_x,dOfs_x_face;
    ublas::vector<ublas::vector<double> > dOfs_x_edge;
    ublas::vector<double> dOfs_X,dOfs_X_face;
    ublas::vector<ublas::vector<double> > dOfs_X_edge;

    ublas::vector<double> fExtNode,fExtFace;
    ublas::vector<ublas::vector<double> > fExtEdge;
    double *Fext_edge[3];

    ublas::matrix<double> kExtNodeNode,kExtFaceNode;
    ublas::vector<ublas::matrix<double> > kExtEdgeNode;
    double *Kext_edge_node[3];

    ublas::matrix<double> kExtNodeFace,kExtFaceFace;
    ublas::vector<ublas::matrix<double> > kExtEdgeFace;
    double *Kext_edge_face[3];

    ublas::vector<ublas::matrix<double> > kExtFaceEdge,kExtNodeEdge;
    ublas::matrix<ublas::matrix<double> > kExtEdgeEdge;
    double *Kext_node_edge[3];
    double *Kext_face_edge[3];
    double *Kext_edge_edge[3][3];

    virtual PetscErrorCode calcTraction();
    virtual PetscErrorCode rHs();
    virtual PetscErrorCode lHs();

    PetscErrorCode preProcess();
    PetscErrorCode operator()();

    PetscErrorCode addForce(int ms_id);
    PetscErrorCode addPreassure(int ms_id);

    struct bCForce {
      ForceCubitBcData data;
      Range tRis;
    };
    map<int,bCForce> mapForce;
    struct bCPreassure {
      PressureCubitBcData data;
      Range tRis;
    };
    map<int,bCPreassure> mapPreassure;
    PetscErrorCode reBaseToFaceLoocalCoordSystem(ublas::matrix<double> &t_glob_nodal);

    boost::ptr_vector<MethodsForOp> methodsOp;

  };

  struct MyTriangleMaterialFE: public MyTriangleSpatialFE {

    MyTriangleMaterialFE(FieldInterface &_mField,Mat _Aij,Vec &_F,double *scale_lhs,double *scale_rhs); 

    PetscErrorCode rHs();
    PetscErrorCode lHs();

  };

  FieldInterface &mField;
  MyTriangleSpatialFE feSpatial;
  MyTriangleMaterialFE feMaterial;

  Tag thScale;

  double *sCale;
  PetscErrorCode setForceScale(double scale) {
      PetscFunctionBegin;
      *sCale = scale;
      PetscFunctionReturn(0);
  }

  MyTriangleSpatialFE& getLoopSpatialFe() { return feSpatial; }
  MyTriangleMaterialFE& getLoopMaterialFe() { return feMaterial; }

  NeummanForcesSurfaceComplexForLazy(FieldInterface &m_field,Mat _Aij,Vec _F,double *scale_lhs,double *scale_rhs);
  NeummanForcesSurfaceComplexForLazy(FieldInterface &m_field,Mat _Aij,Vec _F);

};


}

#endif //__COMPLEX_FOR_LAZY_NEUMANM_FORCES_HPP

