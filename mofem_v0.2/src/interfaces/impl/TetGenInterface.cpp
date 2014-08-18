/** \file TetGenInterface.cpp
 * \brief TetGen inteface for resmeshing and on the fly mesh craetion
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
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

#ifdef WITH_TETGEM

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#endif

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <moab/ParallelComm.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <version.h>
#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>
#include <LoopMethods.hpp>
#include <Core.hpp>

#include <FieldInterface.hpp>

#ifdef WITH_TETGEM

#include <TetGenInterface.hpp>

namespace MoFEM {

PetscErrorCode TetGenInterface::queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMTetGegInterface) {
    *iface = dynamic_cast<TetGenInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<FieldUnknownInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"unknown inteface");

  PetscFunctionReturn(0);
}

PetscErrorCode TetGenInterface::inData(
    Range& ents,tetgenio& in,
    map<EntityHandle,unsigned long>& moab_tetgen_map,
    map<unsigned long,EntityHandle>& tetgen_moab_map) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;

  //PetscErrorCode ierr;
  ErrorCode rval;

  //All indices start from 0
  in.firstnumber = 0;

  in.numberofpoints = ents.subset_by_dimension(0).size();
  in.pointlist = new double[in.numberofpoints * 3];
  rval = m_field.get_moab().get_coords(ents.subset_by_dimension(0),in.pointlist); CHKERR_PETSC(rval);

  Range::iterator it = ents.subset_by_dimension(0).begin();
  for(int ii = 0;it != ents.subset_by_dimension(0).begin(); it++,ii++) {
    tetgen_moab_map[MBVERTEX|(ii<<sizeof(int))] = *it;
    moab_tetgen_map[*it] = MBVERTEX|(ii<<sizeof(int));
  }

  in.numberoftetrahedra = ents.subset_by_type(MBTET).size();
  in.tetrahedronlist = new int[4*ents.subset_by_type(MBTET).size()];
  it = ents.subset_by_type(MBTET).begin();
  for(int ii = 0;it!=ents.subset_by_type(MBTET).end();it++,ii++) {
    int num_nodes;
    const EntityHandle* conn;
    rval = m_field.get_moab().get_connectivity(*it,conn,num_nodes,true); CHKERR_THROW(rval);
    tetgen_moab_map[MBTET|(ii<<sizeof(int))] = *it;
    moab_tetgen_map[*it] = MBTET|(ii<<sizeof(int));
    for(int nn = 0;nn<4;nn++) {
      in.tetrahedronlist[4*ii+nn] = moab_tetgen_map[conn[nn]]>>sizeof(int);
    }
  }

  in.numberoffacets = ents.subset_by_type(MBTRI).size();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  bzero(in.facetmarkerlist,in.numberoffacets*sizeof(int));
  it = ents.subset_by_dimension(2).begin();
  for(int ii = 0;it!=ents.subset_by_dimension(2).end();it++,ii++) {
    tetgenio::facet *f = &(in.facetlist[ii]);
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    tetgenio::polygon *p = &f->polygonlist[0];
    int num_nodes;
    const EntityHandle* conn;
    rval = m_field.get_moab().get_connectivity(*it,conn,num_nodes,true); CHKERR_THROW(rval);
    p->numberofvertices = num_nodes;
    p->vertexlist = new int[p->numberofvertices];
    for(int nn = 0;nn<num_nodes;nn++) {
      p->vertexlist[nn] = moab_tetgen_map[conn[nn]]>>sizeof(int);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TetGenInterface::outData(
  Range& ents,tetgenio& in,tetgenio& out,
  map<EntityHandle,unsigned long>& moab_tetgen_map,
  map<unsigned long,EntityHandle>& tetgen_moab_map) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;

  ErrorCode rval;

  int ii = 0;
  for(;ii<out.numberofpoints;ii++) {
    if(ii<in.numberofpoints) {
      if(memcmp(&in.pointlist[3*ii],&out.pointlist[3*ii],3*sizeof(double)) == 0) {
	unsigned long iii = MBVERTEX|(ii<<sizeof(int));
	if(tetgen_moab_map.find(iii)!=tetgen_moab_map.end()) {
	  continue;
	} else {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
	}
      }
    }
    EntityHandle node;
    rval = m_field.get_moab().create_vertex(&out.pointlist[3*ii],node); CHKERR_PETSC(rval);
    moab_tetgen_map[node] = MBVERTEX|(ii<<sizeof(int));
    tetgen_moab_map[MBVERTEX|(ii<<sizeof(int))] = node;
  }

  ii = 0;
  for(;ii<out.numberoftetrahedra;ii++) {
    if(ii<in.numberoftetrahedra) {
      if(memcmp(&in.tetrahedronlist[4*ii],&out.tetrahedronlist[4*ii],3*sizeof(int)) == 0) {
	unsigned long iii = MBTET|(ii<<sizeof(int));
	if(tetgen_moab_map.find(iii)!=tetgen_moab_map.end()) {
	  continue;
	} else {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
	}
      }
    }
    EntityHandle conn[4];
    for(int nn = 0;nn<4;nn++) {
      int iii = out.tetrahedronlist[4*ii+nn];
      conn[nn] = tetgen_moab_map.at(MBVERTEX|(iii<<sizeof(int)));
    }
    EntityHandle tet;
    rval = m_field.get_moab().create_element(MBTET,conn,4,tet); CHKERR_PETSC(rval);
    moab_tetgen_map[tet] = MBTET|(ii<<sizeof(int));
    tetgen_moab_map[MBTET|(ii<<sizeof(int))] = tet;
  }


  PetscFunctionReturn(0);
}

PetscErrorCode TetGenInterface::setFacetMarkers(Range& ents,int marker,tetgenio& in) {
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

}

#endif //TETGEN
