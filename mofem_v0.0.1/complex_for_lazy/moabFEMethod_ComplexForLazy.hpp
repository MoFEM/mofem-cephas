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

  FEMethod_ComplexForLazy(Interface& _moab,
    double _lambda,double _mu,
    int _verbose = 0): 
    FEMethod_UpLevelStudent(_moab,_verbose),
    type_of_analysis(type_of_analysis),
    lambda(_lambda),mu(_mu),
    eps(1e-12) {
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
      order_edges.resize(6);
      order_faces.resize(4);
      edgeNinvJac.resize(6);
      faceNinvJac.resize(4);
      diff_edgeNinvJac.resize(6);
      diff_faceNinvJac.resize(4);
      Kedgeh_data.resize(6);
      Kfaceh_data.resize(4);
      Kedgeh.resize(6);
      Kfaceh.resize(4);
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      g_NTRI.resize(3*7);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X7,G_TRI_Y7,7); 
      g_TET_W = G_TET_W45;
  }

  vector<double> g_NTET,g_NTRI;
  const double *g_TET_W;
    
  ErrorCode rval;  
  PetscErrorCode ierr;
  ParallelComm* pcomm;

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;

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

  ublas::matrix<double> Khh,Kvolumeh;
  vector<ublas::matrix<double> > Kedgeh_data,Kfaceh_data;
  vector<double*> Kedgeh,Kfaceh;

  ublas::vector<double> dofs_x,dofs_x_volume;
  vector<ublas::vector<double> > dofs_x_edge_data,dofs_x_face_data;
  vector<double*> dofs_x_edge,dofs_x_face;
  
  PetscErrorCode GetIndices();
  PetscErrorCode GetTangent();
  PetscErrorCode GetFint();
  PetscErrorCode GetFext();

  PetscErrorCode OpComplexForLazyStart();

  private:
  vector<int> order_edges;
  vector<int> order_faces;
  int order_volume;

};

}

#endif //__MOABFEMETHOD_COMPLEXFORLAZY_HPP__
