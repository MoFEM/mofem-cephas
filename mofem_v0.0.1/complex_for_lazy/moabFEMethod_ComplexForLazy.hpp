/** \file moabField_Core.hpp
 * \brief Core moabField::FEMethod class for user interface
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __MOABFEMETHOD_COMPLEXFORLAZY_HPP__
#define __MOABFEMETHOD_COMPLEXFORLAZY_HPP__

#include "moabField.hpp"
#include "Core_dataStructures.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** \brief The user interface for NonLineae FE method (tangent is calulated
 * using complex direvatives)
 * 
 * This class give user some data structures and methods on those that
 * structures which could be useful. Method of calulating complex direvatives
 * is numerically inefficient, however it is simple to use.
 *
 * 
*/
struct FEMethod_ComplexForLazy: public FEMethod_UpLevelStudent {

  enum analysis { spatail_analysis = 1, material_analysis = 1<<1 };
  analysis type_of_analysis;

  double lambda,mu;
  void *ptr_matctx;

  double eps;

  FEMethod_ComplexForLazy(Interface& _moab,analysis _type,
    double _lambda,double _mu,int _verbose = 0);

  vector<double> g_NTET,g_NTRI;
  const double *g_TET_W;
    
  ErrorCode rval;  
  PetscErrorCode ierr;
  ParallelComm* pcomm;

  vector<vector<DofIdx> > RowGlob;
  vector<vector<DofIdx> > ColGlob;
  vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
  vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
  vector<vector<ublas::matrix<FieldData> > > rowBMatrices;
  vector<vector<ublas::matrix<FieldData> > > colNMatrices;
  vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
  vector<vector<ublas::matrix<FieldData> > > colBMatrices;

  vector<double*> edgeNinvJac;
  vector<double*> faceNinvJac;
  double *volumeN;
  vector<double*> diff_edgeNinvJac;
  vector<double*> diff_faceNinvJac;
  double *diff_volumeNinvJac;

  //Tangent_hh_hierachical
  ublas::matrix<double> Khh,Kvolumeh;
  vector<ublas::matrix<double> > Kedgeh_data,Kfaceh_data;
  double* Kedgeh[6];
  double* Kfaceh[4];
  //Tangent_hh_hierachical_edge
  vector<ublas::matrix<double> > Khedge_data,Khh_volumeedge_data;
  double* Khedge[6];
  double* Khh_volumeedge[6];
  ublas::matrix<ublas::matrix<double> > Khh_edgeedge_data;
  ublas::matrix<ublas::matrix<double> > Khh_faceedge_data;
  double* Khh_edgeedge[6][6];
  double* Khh_faceedge[4][6];
  //Tangent_hh_hierachical_face
  vector<ublas::matrix<double> > Khface_data,Khh_volumeface_data;
  double* Khface[4];
  double* Khh_volumeface[4];
  ublas::matrix<ublas::matrix<double> > Khh_faceface_data;
  ublas::matrix<ublas::matrix<double> > Khh_edgeface_data;
  double* Khh_faceface[4][4];
  double* Khh_edgeface[6][4];
  //Tangent_hh_hierachical_volume
  ublas::matrix<double> Khvolume;
  vector<ublas::matrix<double> > Khh_edgevolume_data,Khh_facevolume_data;
  double *Khh_edgevolume[6];
  double *Khh_facevolume[4];
  ublas::matrix<double> Khh_volumevolume;
  //Fint
  ublas::vector<double> Fblock_x,Fint_h_volume;
  vector<ublas::vector<double> > Fint_h_edge_data;
  vector<ublas::vector<double> > Fint_h_face_data;
  double* Fint_h_edge[6];
  double* Fint_h_face[4];


  ublas::vector<double> dofs_x,dofs_x_volume;
  vector<ublas::vector<double> > dofs_x_edge_data,dofs_x_face_data;
  vector<double*> dofs_x_edge,dofs_x_face;
  
  PetscErrorCode GetIndices();
  PetscErrorCode GetTangent();
  PetscErrorCode GetFint();

  int face_order;
  vector<int> FaceEdgeOrder;
  vector<DofIdx> FaceIndices;
  vector<FieldData> FaceData;
  vector<DofIdx> NodeIndices;
  vector<FieldData> NodeData;
  vector<vector<DofIdx> > EdgeIndices_data;
  vector<vector<FieldData> > EdgeData_data;
  vector<int> FaceEdgeSense;
  PetscErrorCode GetFaceIndicesAndData(EntityHandle face);
  PetscErrorCode GetTangentExt();
  PetscErrorCode GetFExt();

  PetscErrorCode OpComplexForLazyStart();

  private:
  vector<int> order_edges;
  vector<int> order_faces;
  int order_volume;

};

}

#endif //__MOABFEMETHOD_COMPLEXFORLAZY_HPP__
