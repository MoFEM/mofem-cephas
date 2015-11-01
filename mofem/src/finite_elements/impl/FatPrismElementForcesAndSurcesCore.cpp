/** \file FatPrismElementForcesAndSurcesCore.cpp

\brief Implementation of fat prism element

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

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <FatPrismElementForcesAndSurcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

PetscErrorCode FatPrismElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBPRISM) PetscFunctionReturn(0);

  try {

    ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);
    ierr = getSpacesOnEntities(dataH1TrianglesOnly); CHKERRQ(ierr);
    ierr = getSpacesOnEntities(dataH1TroughThickness); CHKERRQ(ierr);
    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getQuadSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getQuadOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getPrismOrder(dataH1,H1); CHKERRQ(ierr);
      // Triangles only
      ierr = getEdgesSense(dataH1TrianglesOnly); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1TrianglesOnly); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1TrianglesOnly,H1); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1TrianglesOnly,H1); CHKERRQ(ierr);
      // Through thickness
      ierr = getEdgesSense(dataH1TroughThickness); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1TroughThickness,H1); CHKERRQ(ierr);
    }
    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }
    //Hcurl
    if((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }
    //L2
    if((dataH1.spacesOnEntities[MBPRISM]).test(L2)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }

    // get approx. on triangles, i.e. faces 3 and 4
    int nb_gauss_pts_on_faces;
    {
      int order_triangles_only = 1;
      for(unsigned int ee = 0;ee<3;ee++) {
        order_triangles_only = max(
          order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getOrder()
        );
      }
      for(unsigned int ee = 6;ee<dataH1TrianglesOnly.dataOnEntities[MBEDGE].size();ee++) {
        order_triangles_only = max(
          order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getOrder()
        );
      }
      for(unsigned int ff = 0;ff<dataH1TrianglesOnly.dataOnEntities[MBTRI].size();ff++) {
        order_triangles_only = max(
          order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getOrder()
        );
      }
      // integration pts on the triangles surfaces
      int rule = getRuleTrianglesOnly(order_triangles_only);
      if(rule >= 0) {
        nb_gauss_pts_on_faces = gm_rule_size(rule,2);
        gaussPtsTrianglesOnly.resize(3,nb_gauss_pts_on_faces,false);
        ierr = Grundmann_Moeller_integration_points_2D_TRI(
          rule,
          &gaussPtsTrianglesOnly(0,0),
          &gaussPtsTrianglesOnly(1,0),
          &gaussPtsTrianglesOnly(2,0)
        ); CHKERRQ(ierr);
      } else {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        // ierr = setGaussPtsTrianglesOnly(order_triangles_only); CHKERRQ(ierr);
        // nb_gauss_pts_on_faces = gaussPtsTrianglesOnly.size2();
      }
      if(nb_gauss_pts_on_faces == 0) PetscFunctionReturn(0);
      // calculate shape functions
      ierr = shapeFlatPRISMFunctions_H1(
        dataH1TrianglesOnly,
        &gaussPtsTrianglesOnly(0,0),
        &gaussPtsTrianglesOnly(1,0),
        nb_gauss_pts_on_faces
      ); CHKERRQ(ierr);
    }

    // approx. trough prism thickness
    int nb_gauss_pts_through_thickness;
    {
      int order_thickness = 1;
      for(unsigned int ee = 3;ee<6;ee++) {
        order_thickness = max(
          order_thickness,dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getOrder()
        );
      }
      for(unsigned int qq = 0;qq<dataH1TroughThickness.dataOnEntities[MBQUAD].size();qq++) {
        order_thickness = max(
          order_thickness,dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getOrder()
        );
      }
      // integration points
      int rule = getRuleThroughThickness(order_thickness);
      if(rule >= 0) {
        nb_gauss_pts_through_thickness = gm_rule_size(rule,1);
        gaussPtsThroughThickness.resize(2,nb_gauss_pts_through_thickness,false);
        ierr = Grundmann_Moeller_integration_points_1D_EDGE(
          rule,
          &gaussPtsThroughThickness(0,0),
          &gaussPtsThroughThickness(1,0)
        ); CHKERRQ(ierr);
      } else {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        // ierr = setGaussPtsThroughThickness(order_thickness); CHKERRQ(ierr);
        // nb_gauss_pts_through_thickness = gaussPtsThroughThickness.size2();
      }
      if(nb_gauss_pts_through_thickness == 0) PetscFunctionReturn(0);
      // calculate Legendre approx. on edges
      for(unsigned int ee = 0;ee<9;ee++) {
        int sense = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getSense();
        int order = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getOrder()-2;
        dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN().resize(
          nb_gauss_pts_through_thickness,order<0?0:1+order,false
        );
        dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDiffN().resize(
          nb_gauss_pts_through_thickness,order<0?0:1+order,false
        );
        if(order<0) continue;
        double diff_s = 0.5; // s = s(xi), ds/dxi = 0.5, because change of basis
        for(int gg = 0;gg<nb_gauss_pts_through_thickness;gg++) {
          double s = 2*gaussPtsThroughThickness(0,gg)-1; // makes form -1..1
          if(!sense) s *= -1;
          // calculate Legendre polynomials at integration points
          ierr = Legendre_polynomials(
            order,s,&diff_s,
            &dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN()(gg,0),
            &dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDiffN()(gg,0),
            1
          ); CHKERRQ(ierr);
        }
      }
    }

    // Build prism approx.
    int nb_gauss_pts = nb_gauss_pts_on_faces*nb_gauss_pts_through_thickness;
    {
      gaussPts.resize(4,nb_gauss_pts,false);
      {
        int gg = 0;
        for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
          for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
            gaussPts(0,gg) = gaussPtsTrianglesOnly(0,ggf);
            gaussPts(1,gg) = gaussPtsTrianglesOnly(1,ggf);
            gaussPts(2,gg) = gaussPtsThroughThickness(0,ggt);
            gaussPts(3,gg) = gaussPtsTrianglesOnly(2,ggf)*gaussPtsThroughThickness(1,ggt);
          }
        }
      }

      {
        // nodes
        dataH1.dataOnEntities[MBVERTEX][0].getN().resize(nb_gauss_pts,6,false);
        dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(nb_gauss_pts,18);
        for(int dd = 0;dd<6;dd++) {
          int gg = 0;
          for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
            double tri_n = dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,dd);
            double dksi_tri_n = dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*dd+0);
            double deta_tri_n = dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*dd+1);
            for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
              double zeta = gaussPtsThroughThickness(0,ggt);
              if(dd<3) {
                double dzeta = +1;
                dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,dd) =
                tri_n*N_MBEDGE0(zeta);
                dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+0) =
                dksi_tri_n*N_MBEDGE0(zeta);
                dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+1) =
                deta_tri_n*N_MBEDGE0(zeta);
                dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+2) =
                tri_n*dzeta;
              } else {
                double dzeta = -1;
                dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,dd) =
                tri_n*N_MBEDGE1(zeta);
                dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+0) =
                dksi_tri_n*N_MBEDGE1(zeta);
                dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+1) =
                deta_tri_n*N_MBEDGE1(zeta);
                dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+2) =
                tri_n*dzeta;
              }
            }
          }
        }
        const int prism_edge_map[9][2] = {
          {0,1}, {1,2}, {2,0}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {5,3}
        };
        // edges on triangles
        for(int ee = 0;ee<9;ee++) {
          if(ee>2&&ee<6) {
            // through thickness ho approx.
            // cerr << ee << " "
            // << dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getOrder() << " "
            // << dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN() << endl;
            if(dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getOrder()<2) {
              dataH1.dataOnEntities[MBEDGE][ee].getN().resize(nb_gauss_pts,0,false);
              dataH1.dataOnEntities[MBEDGE][ee].getDiffN().resize(nb_gauss_pts,0,false);
              continue;
            }
            int nb_dofs = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN().size2();
            dataH1.dataOnEntities[MBEDGE][ee].getN().resize(nb_gauss_pts,nb_dofs,false);
            dataH1.dataOnEntities[MBEDGE][ee].getDiffN().resize(nb_gauss_pts,3*nb_dofs,false);
            for(int dd = 0;dd<nb_dofs;dd++) {
              int gg = 0;
              for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
                double tri_n =
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][0])*
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][1]);
                double dksi_tri_n =
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*prism_edge_map[ee][0]+0)*
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][1])
                +
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][0])*
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*prism_edge_map[ee][1]+0);
                double deta_tri_n =
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*prism_edge_map[ee][0]+1)*
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][1])
                +
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][0])*
                dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*prism_edge_map[ee][1]+1);
                for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
                  double tri_m = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN()(ggt,dd);
                  double dzeta_tri_m = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDiffN()(ggt,dd);
                  dataH1.dataOnEntities[MBEDGE][ee].getN()(gg,dd) = tri_n*tri_m;
                  dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+0) = dksi_tri_n*tri_m;
                  dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+1) = deta_tri_n*tri_m;
                  dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+2) = tri_n*dzeta_tri_m;
                }
              }
            }
          } else {
            // on triangles ho approx.
            int nb_dofs = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getN().size2();
            dataH1.dataOnEntities[MBEDGE][ee].getN().resize(nb_gauss_pts,nb_dofs,false);
            dataH1.dataOnEntities[MBEDGE][ee].getDiffN().resize(nb_gauss_pts,3*nb_dofs);
            for(int dd = 0;dd<nb_dofs;dd++) {
              int gg = 0;
              for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
                double tri_n = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getN()(ggf,dd);
                double dksi_tri_n = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN()(ggf,2*dd+0);
                double deta_tri_n = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN()(ggf,2*dd+1);
                for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
                  double zeta = gaussPtsThroughThickness(0,ggt);
                  if(ee<3) {
                    double dzeta = 1;
                    dataH1.dataOnEntities[MBEDGE][ee].getN()(gg,dd) =
                    tri_n*N_MBEDGE0(zeta);
                    dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+0) =
                    dksi_tri_n*N_MBEDGE0(zeta);
                    dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+1) =
                    deta_tri_n*N_MBEDGE0(zeta);
                    dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+0) =
                    tri_n*dzeta;
                  } else {
                    double dzeta = -1;
                    dataH1.dataOnEntities[MBEDGE][ee].getN()(gg,dd) =
                    tri_n*N_MBEDGE1(zeta);
                    dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+0) =
                    dksi_tri_n*N_MBEDGE1(zeta);
                    dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+1) =
                    deta_tri_n*N_MBEDGE1(zeta);
                    dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+0) =
                    tri_n*dzeta;
                  }
                }
              }
            }
          }
        }
        // triangles
        for(int ff = 3;ff<5;ff++) {
          int nb_dofs;
          nb_dofs = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getN().size2();
          dataH1.dataOnEntities[MBTRI][ff].getN().resize(nb_gauss_pts,nb_dofs);
          dataH1.dataOnEntities[MBTRI][ff].getDiffN().resize(nb_gauss_pts,3*nb_dofs);
          for(int dd = 0;dd<nb_dofs;dd++) {
            int gg = 0;
            for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
              double tri_n = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getN()(ggf,dd);
              double dksi_tri_n = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getDiffN()(ggf,3*dd+0);
              double deta_tri_n = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getDiffN()(ggf,3*dd+1);
              for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
                double zeta = gaussPtsThroughThickness(0,ggt);
                if(ff == 3) {
                  double dzeta = +1;
                  dataH1.dataOnEntities[MBTRI][ff].getN()(gg,dd) =
                  tri_n*N_MBEDGE0(zeta);
                  dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+0) =
                  dksi_tri_n*N_MBEDGE0(zeta);
                  dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+1) =
                  deta_tri_n*N_MBEDGE0(zeta);
                  dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+0) =
                  tri_n*dzeta;
                } else {
                  double dzeta = -1;
                  dataH1.dataOnEntities[MBTRI][ff].getN()(gg,dd) =
                  tri_n*N_MBEDGE1(zeta);
                  dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+0) =
                  dksi_tri_n*N_MBEDGE1(zeta);
                  dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+1) =
                  deta_tri_n*N_MBEDGE1(zeta);
                  dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+0) =
                  tri_n*dzeta;
                }
              }
            }
          }
        }
        // quads
        {
          int quads_nodes[3*4];
          int quad_order[3];
          double *quad_n[3],*diff_quad_n[3];
          SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
          SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBQUAD,0));
          SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBQUAD,3));
          const EntityHandle *conn_prism;
          int num_nodes_prism;
          EntityHandle ent = fePtr->get_ent();
          rval = mField.get_moab().get_connectivity(ent,conn_prism,num_nodes_prism,true); CHKERR_PETSC(rval);
          for(;siit!=hi_siit;siit++) {
            EntityHandle quad = siit->ent;
            int num_nodes_quad;
            const EntityHandle *conn_quad;
            rval = mField.get_moab().get_connectivity(
              quad,conn_quad,num_nodes_quad,true
            ); CHKERR_PETSC(rval);
            for(int nn = 0;nn<num_nodes_quad;nn++) {
              quads_nodes[4*siit->side_number+nn] = distance(conn_prism,find(conn_prism,&conn_prism[7],conn_quad[nn]));
            }
            int order = dataH1.dataOnEntities[MBQUAD][siit->side_number].getOrder();
            dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().resize(
              nb_gauss_pts_through_thickness,order<0?0:1+order,false
            );
            dataH1.dataOnEntities[MBQUAD][siit->side_number].getDiffN().resize(
              nb_gauss_pts_through_thickness,order<0?0:1+order,false
            );
            quad_order[siit->side_number] = order;
            dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().resize(nb_gauss_pts,NBFACEQUAD_H1(order),false);
            dataH1.dataOnEntities[MBQUAD][siit->side_number].getDiffN().resize(nb_gauss_pts,3*NBFACEQUAD_H1(order),false);
            if(dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().size2()>0) {
              quad_n[siit->side_number] = &*dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().data().begin();
              diff_quad_n[siit->side_number] = &*dataH1.dataOnEntities[MBQUAD][siit->side_number].getDiffN().data().begin();
            }
          }
          ierr = H1_QuadShapeFunctions_MBPRISM(
            quads_nodes,
            quad_order,
            &dataH1.dataOnEntities[MBVERTEX][0].getN()(0,0),
            &dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(0,0),
            quad_n,
            diff_quad_n,
            nb_gauss_pts
          ); CHKERRQ(ierr);
        }
      }
      // prism
      {
        int order = dataH1.dataOnEntities[MBPRISM][0].getOrder();
        double *n  = &dataH1.dataOnEntities[MBVERTEX][0].getN()(0,0);
        double *diff_n = &dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(0,0);
        dataH1.dataOnEntities[MBPRISM][0].getN().resize(nb_gauss_pts,NBVOLUMEPRISM_H1(order),false);
        dataH1.dataOnEntities[MBPRISM][0].getDiffN().resize(nb_gauss_pts,3*NBVOLUMEPRISM_H1(order),false);
        if(NBVOLUMEPRISM_H1(order)>0) {
          ierr = H1_VolumeShapeFunctions_MBPRISM(
            order,
            n,
            diff_n,
            &dataH1.dataOnEntities[MBPRISM][0].getN()(0,0),
            &dataH1.dataOnEntities[MBPRISM][0].getDiffN()(0,0),
            nb_gauss_pts
          ); CHKERRQ(ierr);
        }
      }
    }

    EntityHandle ent = fePtr->get_ent();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
    coords.resize(num_nodes*3,false);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

    normal.resize(6,false);
    ierr = ShapeFaceNormalMBTRI(
      &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &coords[0],
      &normal[0]
    ); CHKERRQ(ierr);
    ierr = ShapeFaceNormalMBTRI(
      &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &coords[9],
      &normal[3]
    ); CHKERRQ(ierr);
    aRea[0] = cblas_dnrm2(3,&normal[0],1)*0.5;
    aRea[1] = cblas_dnrm2(3,&normal[3],1)*0.5;

    coordsAtGaussPtsTrianglesOnly.resize(nb_gauss_pts_on_faces,6,false);
    for(int gg = 0;gg<nb_gauss_pts_on_faces;gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPtsTrianglesOnly(gg,dd) = cblas_ddot(
          3,&dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3
        );
        coordsAtGaussPtsTrianglesOnly(gg,3+dd) = cblas_ddot(
          3,&dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[9+dd],3
        );
      }
    }

    try {

      if(
        dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
        dataPtr->get<FieldName_mi_tag>().end()
      ) {
        hoCoordsAtGaussPtsF3.resize(nb_gauss_pts_on_faces,3,false);
        nOrmals_at_GaussPtF3.resize(nb_gauss_pts_on_faces,3,false);
        tAngent1_at_GaussPtF3.resize(nb_gauss_pts_on_faces,3,false);
        tAngent2_at_GaussPtF3.resize(nb_gauss_pts_on_faces,3,false);
        hoCoordsAtGaussPtsF4.resize(nb_gauss_pts_on_faces,3,false);
        nOrmals_at_GaussPtF4.resize(nb_gauss_pts_on_faces,3,false);
        tAngent1_at_GaussPtF4.resize(nb_gauss_pts_on_faces,3,false);
        tAngent2_at_GaussPtF4.resize(nb_gauss_pts_on_faces,3,false);
        ierr = getEdgesOrder(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisOrder(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldDofs(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        try {
          ierr = opHOCoordsAndNormals.opRhs(dataH1TrianglesOnly); CHKERRQ(ierr);
          ierr = opHOCoordsAndNormals.calculateNormals(); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
      } else {
        hoCoordsAtGaussPtsF3.resize(0,0,false);
        nOrmals_at_GaussPtF3.resize(0,0,false);
        tAngent1_at_GaussPtF3.resize(0,0,false);
        tAngent2_at_GaussPtF3.resize(0,0,false);
        hoCoordsAtGaussPtsF4.resize(0,0,false);
        nOrmals_at_GaussPtF4.resize(0,0,false);
        tAngent1_at_GaussPtF4.resize(0,0,false);
        tAngent2_at_GaussPtF4.resize(0,0,false);
      }
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    const UserDataOperator::OpType types[2] = {
      UserDataOperator::OPROW, UserDataOperator::OPCOL
    };
    vector<string> last_eval_field_name(2);
    DataForcesAndSurcesCore *op_data[2];
    FieldSpace space[2];

    boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
    oit = opPtrVector.begin();
    hi_oit = opPtrVector.end();

    for(;oit!=hi_oit;oit++) {

      oit->setPtrFE(this);

      for(int ss = 0;ss!=2;ss++) {

        string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        BitFieldId data_id = mField.get_field_structure(field_name)->get_id();
        if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data field < %s > on finite element < %s >",
            field_name.c_str(),feName.c_str()
          );
        }

        if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

          space[ss] = mField.get_field_structure(field_name)->get_space();

          switch(space[ss]) {
            case H1:
            op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
            break;
            case HCURL:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
            break;
            case HDIV:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
            break;
            case L2:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
            break;
            case NOFIELD:
            op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
            break;
            case LASTSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
            break;
          }

          if(last_eval_field_name[ss]!=field_name) {

            switch(space[ss]) {
              case H1:
              if(!ss) {
                ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getNodesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case HCURL:
              if(!ss) {
                ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getEdgesOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case HDIV:
              if(!ss) {
                ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTrisOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              if(!ss) {
                ierr = getQuadRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getQuadColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getQuadOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getQuadFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getQuadFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case L2:
              if(!ss) {
                ierr = getPrismRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getPrismColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getPrismOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getPrismFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getPrismFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              break;
              case NOFIELD:
              if(!getNinTheLoop()) {
                // NOFIELD data are the same for each element, can be retrieved only once
                if(!ss) {
                  ierr = getNoFieldRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
                } else {
                  ierr = getNoFieldColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
                }
                ierr = getNoFieldFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              break;
              case LASTSPACE:
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
              break;
            }
            last_eval_field_name[ss]=field_name;

          }
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROW) {
        try {
          ierr = oit->opRhs(*op_data[0]); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPCOL) {
        try {
          ierr = oit->opRhs(*op_data[1]); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROWCOL) {
        try {
          ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

}
