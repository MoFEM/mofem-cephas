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
//#ifdef REAL
  //#undef REAL
//#endif

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
  Range::iterator it;

  //All indices start from 0
  in.firstnumber = 0;
  in.numberofpoints = ents.subset_by_dimension(0).size();
  if(ents.subset_by_dimension(0).size()>0) {
    in.pointlist = new double[in.numberofpoints * 3];
    Range points = ents.subset_by_dimension(0);
    rval = m_field.get_moab().get_coords(points,in.pointlist); CHKERR_PETSC(rval);
    it = points.begin();
    for(int ii = 0;it != points.end(); it++,ii++) {
      unsigned long iii = MBVERTEX|(ii<<sizeof(int));
      tetgen_moab_map[iii] = *it;
      moab_tetgen_map[*it] = iii;
    }
  }

  in.numberoftetrahedra = ents.subset_by_type(MBTET).size();
  if(in.numberoftetrahedra>0) {
    in.tetrahedronlist = new int[4*ents.subset_by_type(MBTET).size()];
    Range tets = ents.subset_by_type(MBTET);
    it = tets.begin();
    for(int ii = 0;it!=tets.end();it++,ii++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = m_field.get_moab().get_connectivity(*it,conn,num_nodes,true); CHKERR_THROW(rval);
      tetgen_moab_map[MBTET|(ii<<sizeof(int))] = *it;
      moab_tetgen_map[*it] = MBTET|(ii<<sizeof(int));
      for(int nn = 0;nn<4;nn++) {
	in.tetrahedronlist[4*ii+nn] = moab_tetgen_map[conn[nn]]>>sizeof(int);
      }
    }
  }

  in.numberoffacets = ents.subset_by_dimension(2).size();
  if(in.numberoffacets>0) {
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    bzero(in.facetmarkerlist,in.numberoffacets*sizeof(int));
    Range faces = ents.subset_by_dimension(2);
    it = faces.begin();
    for(int ii = 0;it!=faces.end();it++,ii++) {
      int iii = m_field.get_moab().type_from_handle(*it)|ii<<sizeof(int);
      tetgen_moab_map[iii] = *it;
      moab_tetgen_map[*it] = iii;
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
    ents.insert(node);
  }

  ii = 0;
  for(;ii<out.numberoftetrahedra;ii++) {
    if(ii<in.numberoftetrahedra) {
      if(memcmp(&in.tetrahedronlist[4*ii],&out.tetrahedronlist[4*ii],4*sizeof(int)) == 0) {
	unsigned long iii = MBTET|(ii<<sizeof(int));
	if(tetgen_moab_map.find(iii)!=tetgen_moab_map.end()) {
          ents.insert(tetgen_moab_map[iii]);
	  continue;
	} else {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
	}
      }
    }
    EntityHandle conn[4];
    for(int nn = 0;nn<4;nn++) {
      int iii = out.tetrahedronlist[4*ii+nn];
      if(tetgen_moab_map.find(MBVERTEX|(iii<<sizeof(int)))==tetgen_moab_map.end()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
      }
      conn[nn] = tetgen_moab_map.at(MBVERTEX|(iii<<sizeof(int)));
    }
    EntityHandle tet;
    rval = m_field.get_moab().create_element(MBTET,conn,4,tet); CHKERR_PETSC(rval);
    moab_tetgen_map[tet] = MBTET|(ii<<sizeof(int));
    tetgen_moab_map[MBTET|(ii<<sizeof(int))] = tet;
    ents.insert(tet);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TetGenInterface::outData(
  BitRefLevel bit,tetgenio& in,tetgenio& out,
  map<EntityHandle,unsigned long>& moab_tetgen_map,
  map<unsigned long,EntityHandle>& tetgen_moab_map) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Range ents;
  ierr = outData(ents,in,out,moab_tetgen_map,tetgen_moab_map); CHKERRQ(ierr);
  //cerr << ents.size() << endl;
  FieldInterface& m_field = cOre;
  ierr = m_field.seed_ref_level_3D(ents.subset_by_type(MBTET),bit); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TetGenInterface::setFacetMarkers(Range& ents,map<EntityHandle,unsigned long>& moab_tetgen_map,int marker,tetgenio& in) {
  PetscFunctionBegin;

  Range::iterator it = ents.subset_by_dimension(2).begin();
  for(;it!=ents.subset_by_dimension(2).end();it++) {
    int ii = moab_tetgen_map[*it]>>sizeof(int);
    in.facetmarkerlist[ii] = marker;
  }

  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::setFacetMarkers(EntityHandle ent[],int nb_ents,map<EntityHandle,unsigned long>& moab_tetgen_map,int marker,tetgenio& in) {
  PetscFunctionBegin;

  int jj = 0;
  for(;jj<nb_ents;jj++) {
    int ii = moab_tetgen_map[ent[jj]]>>sizeof(int);
    in.facetmarkerlist[ii] = marker;
  }

  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::setFaceCubitSideSetMarkers(map<EntityHandle,unsigned long>& moab_tetgen_map,tetgenio& in) {
  PetscFunctionBegin;

  ErrorCode rval;
  FieldInterface& m_field = cOre;

  int shift = 0;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    Range faces;
    rval = m_field.get_moab().get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
    Range::iterator it = faces.begin();
    for(;it!=faces.end();it++) {
      if(moab_tetgen_map.find(*it)!=moab_tetgen_map.end()) {
	int ii = moab_tetgen_map[*it]>>sizeof(int);
	in.facetmarkerlist[ii] |= 1<<shift;
      }
    }
    shift++;
  }

  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::getFaceCubitSideSetMarkers(map<EntityHandle,unsigned long>& tetgen_moab_map,tetgenio& out) {
  PetscFunctionBegin;

  ErrorCode rval;
  FieldInterface& m_field = cOre;

  int shift = 0;
  map<int,EntityHandle> shift_meshset_map;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,sit)) {
    shift_meshset_map[shift] = sit->get_meshset();
    shift++;
  }

  int ii = 0;
  for(;ii<out.numberoftrifaces;ii++) {
    for(int sh = 0;sh<shift;sh++) {
      if(out.trifacemarkerlist[ii]&(1<<sh)) {
	EntityHandle conn[3];
	for(int nn = 0;nn<3;nn++) {
	  int iii = MBVERTEX|(out.trifacelist[3*ii+nn]<<sizeof(int));
	  if(tetgen_moab_map.find(iii) == tetgen_moab_map.end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
	  } else {
	    conn[nn] = tetgen_moab_map[iii];
	  }
	}
	Range face;
	rval = m_field.get_moab().get_adjacencies(conn,3,2,true,face); CHKERR_PETSC(rval);
	if(face.size()!=1) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
	}
	if(shift_meshset_map.find(sh) == shift_meshset_map.end()) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"unknown boundary marker and SideSet");
	}
	rval = m_field.get_moab().add_entities(shift_meshset_map[sh],face); CHKERR(rval);
      }
    }
  }

  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::tetRahedralize(char switches[],tetgenio& in,tetgenio& out) {
  PetscFunctionBegin;

  tetgenbehavior a;
  a.parse_commandline(switches);
  tetrahedralize(&a,&in,&out);

  PetscFunctionReturn(0);
}



}

#endif //TETGEN
