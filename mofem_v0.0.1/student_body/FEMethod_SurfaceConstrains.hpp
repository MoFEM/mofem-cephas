/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __MOABSURFACECONSTRAINS_HPP__
#define __MOABSURFACECONSTRAINS_HPP__

#include "FieldInterface.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

struct C_SURFACE_FEMethod:public FieldInterface::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Mat C;
  Tag th_material_normal;
  vector<double> diffNTRI;
  vector<double> g_NTRI3;
  const double *G_TRI_W;
  
  void run_in_constructor();

  string lambda_field_name;
 
  C_SURFACE_FEMethod(Interface& _moab,Mat _C,string _lambda_field_name,int _verbose = 0);
  C_SURFACE_FEMethod(Interface& _moab,Mat _C,int _verbose = 0);

  PetscErrorCode preProcess();

  virtual PetscErrorCode SaveConstrainOnTags();

  ublas::bounded_matrix<double,3,9 > C_MAT_ELEM;
  ublas::vector<DofIdx,ublas::bounded_array<DofIdx,9> > ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_lambda_data;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map0;
  ublas::vector<double,ublas::bounded_array<double,9> > coords;

  virtual PetscErrorCode Integrate();
  PetscErrorCode operator()();
  PetscErrorCode postProcess();

};

struct g_SURFACE_FEMethod: public C_SURFACE_FEMethod {

  Vec g;
  g_SURFACE_FEMethod(Interface& _moab,Vec _g,string _lambda_field_name,int _verbose = 0); 
  g_SURFACE_FEMethod(Interface& _moab,Vec _g,int _verbose = 0); 

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;
  PetscErrorCode Integrate();

};

struct C_CORNER_FEMethod:public FieldInterface::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Mat C;
  C_CORNER_FEMethod(Interface& _moab,Mat _C,int _verbose = 0); 

  PetscErrorCode preProcess(); 
  
  ublas::bounded_matrix<double,3,3 > C_MAT_ELEM;
  ublas::vector<DofIdx,ublas::bounded_array<DofIdx,3> > ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_lambda_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
  ublas::vector<double,ublas::bounded_array<double,3> > coords; 

  virtual PetscErrorCode Integrate();

  PetscErrorCode operator()();

  PetscErrorCode postProcess();

};

struct g_CORNER_FEMethod: public C_CORNER_FEMethod {
  
  Vec g;
  g_CORNER_FEMethod(Interface& _moab,Vec _g,int _verbose = 0); 

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;

  PetscErrorCode Integrate();

};


}

#endif //__MOABSURFACECONSTRAINS_HPP__
