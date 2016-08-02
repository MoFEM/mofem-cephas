/** \file NetGenInterface.cpp
 * \brief NetGen inteface for resmeshing and on the fly mesh craetion
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

#include <Includes.hpp>

//SRC APPROXIMATION
#include <config.h>
#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

//SRC/MULTI-INDICES
#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <SeriesMultiIndices.hpp>

//SRC/INTERFACES
#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#ifdef WITH_NETGEN

/*

  Interface to the netgen meshing kernel

*/
#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>
#include <../visualization/soldata.hpp>

#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif

//#include <nginterface.h>

namespace nglib {
#include <nglib.h>
}
using namespace nglib;
using namespace netgen;

#include <NetGenInterface.hpp>

#include <moab/Skinner.hpp>

namespace MoFEM {

PetscErrorCode NetGenInterface::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMNetGegInterface) {
    *iface = dynamic_cast<NetGenInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<UnknownInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown inteface");
  PetscFunctionReturn(0);
}

PetscErrorCode NetGenInterface::stlSetSurfaceTriangles(Ng_STL_Geometry *stl_geom,Range &ents,double *nv,int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  //PetscErroCode ierr;

  MoFEM::Interface& m_field = cOre;
  Range tets = ents.subset_by_type(MBTET);
  Range tris = ents.subset_by_type(MBTRI);
  Range::iterator tit = tris.begin();
  for(;tit!=tris.end();tit++) {

    int num_nodes;
    const EntityHandle* conn;
    rval = m_field.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    double coords[9];
    rval = m_field.get_moab().get_coords(conn,3,coords); CHKERRQ_MOAB(rval);

    int first = 0,second = 1;
    Range adj_tet;
    rval = m_field.get_moab().get_adjacencies(&*tit,1,3,false,adj_tet); CHKERRQ_MOAB(rval);
    adj_tet = intersect(adj_tet,tets);
    if(!adj_tet.empty()) {
      int side_number;
      int sense;
      int offset;
      rval = m_field.get_moab().side_number(adj_tet[0],*tit,side_number,sense,offset); CHKERRQ_MOAB(rval);
      if(sense == -1) {
	first = 1;
	second = 0;
      }
    }

    // fils STL Geometry
    // positive orientation
    // normal vector may be null-pointer
    Ng_STL_AddTriangle(stl_geom,&coords[3*first],&coords[3*second],&coords[3*2],nv);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode NetGenInterface::stlSetSurfaceEdges(Ng_STL_Geometry *stl_geom,Range &ents,int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  //PetscErroCode ierr;

  MoFEM::Interface& m_field = cOre;
  Range edges = ents.subset_by_type(MBEDGE);
  Range::iterator eit = edges.begin();
  for(;eit!=edges.end();eit++) {

    int num_nodes;
    const EntityHandle* conn;
    rval = m_field.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    double coords[6];
    rval = m_field.get_moab().get_coords(conn,2,coords); CHKERRQ_MOAB(rval);

    Ng_STL_AddEdge(stl_geom,&coords[3*0],&coords[3*1]);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode NetGenInterface::setPoints(Ng_Mesh *mesh,std::vector<EntityHandle> &pts,int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  MoFEM::Interface& m_field = cOre;

  std::vector<double> coords(pts.size()*3);
  rval = m_field.get_moab().get_coords(&*pts.begin(),pts.size(),&*coords.begin()); CHKERRQ_MOAB(rval);

  std::vector<EntityHandle>::iterator pit = pts.begin();
  for(int nn = 0;pit!=pts.end();pit++,nn++) {
    Ng_AddPoint (mesh,&coords[3*nn]);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NetGenInterface::setSurfaceElements(Ng_Mesh *mesh,std::vector<EntityHandle> &pts,std::vector<EntityHandle> &elms,Range *tets,int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  //PetscErroCode ierr;

  std::map<EntityHandle,int> map_pts;
  std::vector<EntityHandle>::iterator pit = pts.begin();
  for(int nn = 1;pit!=pts.end();pit++,nn++) {
    map_pts[*pit] = nn;
  }

  MoFEM::Interface& m_field = cOre;
  std::vector<EntityHandle>::iterator eit = elms.begin();
  for(;eit!=elms.end();eit++) {

    int num_nodes;
    const EntityHandle* conn;
    rval = m_field.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);

    int order[] = { 0,1,2 };
    Range adj_tet;
    rval = m_field.get_moab().get_adjacencies(&*eit,1,3,false,adj_tet); CHKERRQ_MOAB(rval);
    adj_tet = intersect(adj_tet,*tets);
    if(!adj_tet.empty()) {
      int side_number;
      int sense;
      int offset;
      rval = m_field.get_moab().side_number(adj_tet[0],*eit,side_number,sense,offset); CHKERRQ_MOAB(rval);
      if(sense == -1) {
	order[0] = 1;
	order[1] = 0;
      }
    }

    int pi[NG_SURFACE_ELEMENT_MAXPOINTS];
    for(int nn = 0;nn<num_nodes;nn++) {
      pi[nn] = map_pts[conn[order[nn]]];
    }

    switch(m_field.get_moab().type_from_handle(*eit)) {
      case MBTRI:
	Ng_AddSurfaceElement (mesh,NG_TRIG,pi);
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this element is not yet implemented or not surface element");
    }


  }

  PetscFunctionReturn(0);
}
PetscErrorCode NetGenInterface::getPoints(Ng_Mesh *mesh,std::vector<EntityHandle> &pts) {
  PetscFunctionBegin;
  MoFEM::Interface& m_field = cOre;
  ErrorCode rval;

  int nb_pts = Ng_GetNP(mesh);
  pts.resize(nb_pts,0);
  for(int nn = 1;nn<=nb_pts;nn++) {
    double x[3];
    Ng_GetPoint(mesh,nn,x);
    if(pts[nn-1]!=0) {
      double coords[3];
      rval = m_field.get_moab().get_coords(&pts[nn-1],1,coords); CHKERRQ_MOAB(rval);
      cblas_daxpy(3,-1,x,1,coords,1);
      double d = cblas_dnrm2(3,coords,1);
      const double eps = 1e-12;
      if(d < eps) {
	continue;
      }
    }
    EntityHandle node;
    rval = m_field.get_moab().create_vertex(x,node); CHKERRQ_MOAB(rval);
    pts[nn-1] = node;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NetGenInterface::getSurfaceElements(Ng_Mesh *mesh,std::vector<EntityHandle> &pts,std::vector<EntityHandle> &elms) {
  PetscFunctionBegin;
  MoFEM::Interface& m_field = cOre;
  ErrorCode rval;
  //PetscErrorCode ierr;

  Tag th_geom_info;
  int def_marker = 0;
  rval = m_field.get_moab().tag_get_handle(
    "NETGEN_GEOMINFO",1,MB_TYPE_INTEGER,th_geom_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERRQ_MOAB(rval);

  int ne;
  ne = Ng_GetNSE(mesh);
  for(int ee = 1;ee<=ne;ee++) {
    int pi[NG_SURFACE_ELEMENT_MAXPOINTS];
    int matnum = 0;
    Ng_Surface_Element_Type type = Ng_GetElement_2D(mesh,ee,pi,&matnum);
    //by pass nglib interface
    const Element2d &el = ((Mesh*)mesh)->SurfaceElement(ee);
    EntityHandle elem;
    EntityHandle conn[NG_SURFACE_ELEMENT_MAXPOINTS];
    switch(type) {
      case NG_TRIG:
	for(int nn = 0;nn<3;nn++) {
	  conn[nn] = pts[pi[nn]-1];
	}
	rval = m_field.get_moab().create_element(MBTRI,conn,3,elem); CHKERRQ_MOAB(rval);
	elms.push_back(elem);
	for(int nn = 0;nn<3;nn++) {
	  //std::cerr << "GeomInfo " << el.GeomInfoPi(nn+1) << std::endl;
	  int stl_tri_num = el.GeomInfoPi(nn+1).trignum;
	  rval = m_field.get_moab().tag_set_data(th_geom_info,&conn[nn],1,&stl_tri_num); CHKERRQ_MOAB(rval);
	}
	break;
      case NG_QUAD:
      case NG_TRIG6:
      case NG_QUAD6:
      case NG_QUAD8:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this element is not yet implemented");
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NetGenInterface::getVolumeElements(Ng_Mesh *mesh,std::vector<EntityHandle> &pts,std::vector<EntityHandle> &elms) {
  PetscFunctionBegin;

  MoFEM::Interface& m_field = cOre;
  ErrorCode rval;

  int ne;
  ne = Ng_GetNE(mesh);
  for(int ee = 1;ee<=ne;ee++) {
    int pi[NG_SURFACE_ELEMENT_MAXPOINTS];
    Ng_Volume_Element_Type type = Ng_GetVolumeElement(mesh,ee,pi);
    EntityHandle elem;
    EntityHandle conn[NG_VOLUME_ELEMENT_MAXPOINTS];
    EntityHandle conn0;
    switch(type) {
      case NG_TET:
	for(int nn = 0;nn<4;nn++) {
	  conn[nn] = pts[pi[nn]-1];
	}
	conn0 = conn[0];
	conn[0] = conn[1];
	conn[1] = conn0;
	rval = m_field.get_moab().create_element(MBTET,conn,4,elem); CHKERRQ_MOAB(rval);
	elms.push_back(elem);
	break;
      case NG_PYRAMID:
      case NG_PRISM:
      case NG_TET10:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this element is not yet implemented");
    }
  }

  PetscFunctionReturn(0);
}



}

#endif //NETGEN
