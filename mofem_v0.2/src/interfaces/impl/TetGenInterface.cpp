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

#include <moab/Skinner.hpp>

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
      rval = m_field.get_moab().get_connectivity(*it,conn,num_nodes,true); CHKERR_PETSC(rval);
      tetgen_moab_map[MBTET|(ii<<sizeof(int))] = *it;
      moab_tetgen_map[*it] = MBTET|(ii<<sizeof(int));
      for(int nn = 0;nn<4;nn++) {
	in.tetrahedronlist[4*ii+nn] = moab_tetgen_map[conn[nn]]>>sizeof(int);
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
	} /*else {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
	}*/
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
    unsigned long iii = MBTET|(ii<<sizeof(int));
    if(ii<in.numberoftetrahedra) {
      if(memcmp(&in.tetrahedronlist[4*ii],&out.tetrahedronlist[4*ii],4*sizeof(int)) == 0) {
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
      int nnn = out.tetrahedronlist[4*ii+nn];
      if(tetgen_moab_map.find(MBVERTEX|(nnn<<sizeof(int)))==tetgen_moab_map.end()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency between TetGen and MoAB");
      }
      conn[nn] = tetgen_moab_map.at(MBVERTEX|(nnn<<sizeof(int)));
    }
    EntityHandle tet;
    rval = m_field.get_moab().create_element(MBTET,conn,4,tet); CHKERR_PETSC(rval);
    moab_tetgen_map[tet] = iii;
    tetgen_moab_map[iii] = tet;
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
PetscErrorCode TetGenInterface::setFaceData(
  vector<pair<Range,int> >& markers,
  tetgenio& in,
  map<EntityHandle,unsigned long>& moab_tetgen_map,
  map<unsigned long,EntityHandle>& tetgen_moab_map) {
  PetscFunctionBegin;
  ErrorCode rval;
  FieldInterface& m_field = cOre;
  in.numberoffacets = markers.size();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  vector<pair<Range,int> >::iterator mit = markers.begin();
  for(int ii = 0;mit!=markers.end();mit++,ii++) {
    //cerr << "AAAAAAAAA\n";
    in.facetmarkerlist[ii] = mit->second;
    Range& faces = mit->first;
    tetgenio::facet *f = &(in.facetlist[ii]);
    f->numberofpolygons = faces.size();
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    Range::iterator it = faces.begin();
    for(int jj = 0;it!=faces.end();it++,jj++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = m_field.get_moab().get_connectivity(*it,conn,num_nodes,true); CHKERR_PETSC(rval);
      tetgenio::polygon *p = &(f->polygonlist[jj]);
      p->numberofvertices = num_nodes;
      p->vertexlist = new int[p->numberofvertices];
      for(int nn = 0;nn<num_nodes;nn++) {
	p->vertexlist[nn] = moab_tetgen_map[conn[nn]]>>sizeof(int);
	//cerr << p->vertexlist[nn] << "( " << conn[nn] << ") ";
      }
      //cerr << endl;
    }
    //holes
    f->numberofholes = 0;
    f->holelist = NULL;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::getTiangleAttributes(map<EntityHandle,unsigned long>& tetgen_moab_map,tetgenio& out) {
  PetscFunctionBegin;
  ErrorCode rval;
  FieldInterface& m_field = cOre;
  Tag th_attribute;
  int def_marker = -1;
  rval = m_field.get_moab().tag_get_handle(
    "TETGEN_TRIANGLEATTRIBUTE",1,MB_TYPE_INTEGER,th_attribute,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_PETSC(rval); 
  int ii = 0;
  for(;ii<out.numberoftrifaces;ii++) {
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
    rval = m_field.get_moab().tag_set_data(th_attribute,&*face.begin(),1,&out.trifacemarkerlist[ii]); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::setReginData(vector<pair<Range,int> >& regions,tetgenio& in) {
  PetscFunctionBegin;
  ErrorCode rval;
  FieldInterface& m_field = cOre;
  int nb = 0;
  vector<pair<Range,int> >::iterator it = regions.begin();
  for(int ii = 0;it!=regions.end();it++,ii++) {
    nb += it->first.size();
  }
  in.numberofregions = nb;
  in.regionlist = new double[5*nb];
  int kk = 0;
  it = regions.begin();
  for(int ii = 0;it!=regions.end();it++,ii++) {
    Range& ents = it->first;
    Range::iterator iit = ents.begin();
    for(int jj = 0;iit!=ents.end();iit++,jj++) {
      double coords[3];
      switch(m_field.get_moab().type_from_handle(*iit)) {
	case MBTET: {
	    int num_nodes;
	    const EntityHandle* conn;
	    rval = m_field.get_moab().get_connectivity(*iit,conn,num_nodes,true); CHKERR_PETSC(rval);
	    double _coords[12];
	    rval = m_field.get_moab().get_coords(conn,num_nodes,_coords); CHKERR_PETSC(rval);
	    coords[0] = (_coords[0] + _coords[3] + _coords[6] + _coords[9 ])/4.;
	    coords[1] = (_coords[1] + _coords[4] + _coords[7] + _coords[10])/4.;
	    coords[2] = (_coords[2] + _coords[5] + _coords[8] + _coords[11])/4.;
	  }
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      for(int nn = 0;nn<3;nn++) {
	in.regionlist[kk++] = coords[nn];
      }
      in.regionlist[kk++] = ii;
      in.regionlist[kk++] = ii;
    }
  }

  PetscFunctionReturn(0);
}
PetscErrorCode TetGenInterface::getReginData(map<EntityHandle,unsigned long>& tetgen_moab_map,tetgenio& out) {
  PetscFunctionBegin;
  ErrorCode rval;
  FieldInterface& m_field = cOre;
  int nbattributes = out.numberoftetrahedronattributes;
  Tag th_region;
  double def_marker = 0;
  rval = m_field.get_moab().tag_get_handle(
    "TETGEN_REGION",nbattributes,MB_TYPE_DOUBLE,
    th_region,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_PETSC(rval); 
  int ii = 0;
  for(;ii<out.numberoftetrahedra;ii++) {
    int jj = 0;
    for(;jj<nbattributes;jj++) {
      double id = out.tetrahedronattributelist[ii*nbattributes+jj];
      //cerr << id << endl;
      int iii = MBTET|(ii<<sizeof(int));
      if(tetgen_moab_map.find(iii)==tetgen_moab_map.end()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
	  "data inconsistency between TetGen and MoAB");
      }
      EntityHandle ent = tetgen_moab_map[iii];
      rval = m_field.get_moab().tag_set_data(th_region,&ent,1,&id); CHKERR_PETSC(rval);
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
PetscErrorCode TetGenInterface::load_poly(char file_name[],tetgenio& in) {
  PetscFunctionBegin;
  if(!in.load_poly(file_name)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
      "can not read TetGen poly file");
  }
  PetscFunctionReturn(0);
}


}

#endif //TETGEN
