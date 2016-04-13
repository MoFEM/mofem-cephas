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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FatPrismElementForcesAndSurcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  // #include <gm_rule.h>
  #include <quad.h>
#ifdef __cplusplus
}
#endif


namespace MoFEM {

PetscErrorCode FatPrismElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBPRISM) PetscFunctionReturn(0);

  try {

    ierr = getSpacesAndBaseOnEntities(dataH1); CHKERRQ(ierr);
    ierr = getSpacesAndBaseOnEntities(dataH1TrianglesOnly); CHKERRQ(ierr);
    ierr = getSpacesAndBaseOnEntities(dataH1TroughThickness); CHKERRQ(ierr);
    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getQuadSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getQuadDataOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getPrismDataOrder(dataH1,H1); CHKERRQ(ierr);
      // Triangles only
      ierr = getEdgesSense(dataH1TrianglesOnly); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1TrianglesOnly); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataH1TrianglesOnly,H1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataH1TrianglesOnly,H1); CHKERRQ(ierr);
      // Through thickness
      ierr = getEdgesSense(dataH1TroughThickness); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataH1TroughThickness,H1); CHKERRQ(ierr);
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

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  // get approx. on triangles, i.e. faces 3 and 4
  int nb_gauss_pts_on_faces;
  try {
    int order_triangles_only = 1;
    int valid_edges[] = { 1,1,1, 0,0,0, 1,1,1 };
    for(unsigned int ee = 0;ee<9;ee++) {
      if(!valid_edges[ee]) continue;
      order_triangles_only = max(
        order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getDataOrder()
      );
    }
    for(unsigned int ff = 3;ff<=4;ff++) {
      order_triangles_only = max(
        order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getDataOrder()
      );
    }
    for(unsigned int qq = 0;qq<3;qq++) {
      order_triangles_only = max(
        order_triangles_only,dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getDataOrder()
      );
    }
    order_triangles_only = max(
      order_triangles_only,dataH1TroughThickness.dataOnEntities[MBPRISM][0].getDataOrder()
    );
    // integration pts on the triangles surfaces
    int rule = getRuleTrianglesOnly(order_triangles_only);
    if(rule >= 0) {
      if(rule<QUAD_2D_TABLE_SIZE) {
        if(QUAD_2D_TABLE[rule]->dim!=2) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong dimension");
        }
        if(QUAD_2D_TABLE[rule]->order<rule) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",
            QUAD_2D_TABLE[rule]->order,rule
          );
        }
        nb_gauss_pts_on_faces = QUAD_2D_TABLE[rule]->npoints;
        gaussPtsTrianglesOnly.resize(3,nb_gauss_pts_on_faces,false);
        cblas_dcopy(
          nb_gauss_pts_on_faces,&QUAD_2D_TABLE[rule]->points[1],3,&gaussPtsTrianglesOnly(0,0),1
        );
        cblas_dcopy(
          nb_gauss_pts_on_faces,&QUAD_2D_TABLE[rule]->points[2],3,&gaussPtsTrianglesOnly(1,0),1
        );
        cblas_dcopy(
          nb_gauss_pts_on_faces,QUAD_2D_TABLE[rule]->weights,1,&gaussPtsTrianglesOnly(2,0),1
        );
      } else {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
          rule,QUAD_2D_TABLE_SIZE
        );
        nb_gauss_pts_on_faces = 0;
      }
      // nb_gauss_pts_on_faces = gm_rule_size(rule,2);
      // gaussPtsTrianglesOnly.resize(3,nb_gauss_pts_on_faces,false);
      // ierr = Grundmann_Moeller_integration_points_2D_TRI(
      //   rule,
      //   &gaussPtsTrianglesOnly(0,0),
      //   &gaussPtsTrianglesOnly(1,0),
      //   &gaussPtsTrianglesOnly(2,0)
      // ); CHKERRQ(ierr);
    } else {
      ierr = setGaussPtsTrianglesOnly(order_triangles_only); CHKERRQ(ierr);
      nb_gauss_pts_on_faces = gaussPtsTrianglesOnly.size2();
    }
    if(nb_gauss_pts_on_faces == 0) PetscFunctionReturn(0);
    // calculate shape functions
    ierr = shapeFlatPRISMFunctions_H1(
      dataH1TrianglesOnly,
      &gaussPtsTrianglesOnly(0,0),
      &gaussPtsTrianglesOnly(1,0),
      nb_gauss_pts_on_faces
    ); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  // approx. trough prism thickness
  int nb_gauss_pts_through_thickness;
  try {
    int order_thickness = 1;
    for(unsigned int ee = 3;ee<=5;ee++) {
      order_thickness = max(
        order_thickness,dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDataOrder()
      );
    }
    for(unsigned int qq = 0;qq<3;qq++) {
      order_thickness = max(
        order_thickness,dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getDataOrder()
      );
    }
    order_thickness = max(
      order_thickness,dataH1TroughThickness.dataOnEntities[MBPRISM][0].getDataOrder()
    );
    // integration points
    int rule = getRuleThroughThickness(order_thickness);
    if(rule >= 0) {
      if(rule<QUAD_1D_TABLE_SIZE) {
        if(QUAD_1D_TABLE[rule]->dim!=1) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong dimension");
        }
        if(QUAD_1D_TABLE[rule]->order<rule) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",
            QUAD_1D_TABLE[rule]->order,rule
          );
        }
        nb_gauss_pts_through_thickness = QUAD_1D_TABLE[rule]->npoints;
        gaussPtsThroughThickness.resize(2,nb_gauss_pts_through_thickness,false);
        cblas_dcopy(
          nb_gauss_pts_through_thickness,&QUAD_1D_TABLE[rule]->points[1],2,&gaussPtsThroughThickness(0,0),1
        );
        cblas_dcopy(
          nb_gauss_pts_through_thickness,QUAD_1D_TABLE[rule]->weights,1,&gaussPtsThroughThickness(1,0),1
        );
      } else {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
          rule,QUAD_1D_TABLE_SIZE
        );
        nb_gauss_pts_through_thickness = 0;
      }
      // nb_gauss_pts_through_thickness = gm_rule_size(rule,1);
      // gaussPtsThroughThickness.resize(2,nb_gauss_pts_through_thickness,false);
      // ierr = Grundmann_Moeller_integration_points_1D_EDGE(
      //   rule,
      //   &gaussPtsThroughThickness(0,0),
      //   &gaussPtsThroughThickness(1,0)
      // ); CHKERRQ(ierr);
      // cerr << gaussPtsThroughThickness << endl;
    } else {
      ierr = setGaussPtsThroughThickness(order_thickness); CHKERRQ(ierr);
      nb_gauss_pts_through_thickness = gaussPtsThroughThickness.size2();
    }
    if(nb_gauss_pts_through_thickness == 0) PetscFunctionReturn(0);
    // calculate Legendre approx. on edges
    for(unsigned int ee = 3;ee<=5;ee++) {
      int sense = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getSense();
      // cerr << "sense " << sense << endl;
      int order = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDataOrder()-2;
      dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN().resize(
        nb_gauss_pts_through_thickness,order<0?0:1+order,false
      );
      dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDiffN().resize(
        nb_gauss_pts_through_thickness,order<0?0:1+order,false
      );
      if(order<0) continue;
      double diff_s = 2.; // s = s(xi), ds/dxi = 2., because change of basis
      for(int gg = 0;gg<nb_gauss_pts_through_thickness;gg++) {
        double s = 2*gaussPtsThroughThickness(0,gg)-1; // makes form -1..1
        if(!sense) {
          s *= -1;
          diff_s *= -1;
        }
        // calculate Legendre polynomials at integration points
        ierr = Legendre_polynomials(
          order,s,&diff_s,
          &dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN()(gg,0),
          &dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDiffN()(gg,0),
          1
        ); CHKERRQ(ierr);
      }
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  // Build prism approx.
  int nb_gauss_pts = nb_gauss_pts_on_faces*nb_gauss_pts_through_thickness;
  try {
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
    // nodes
    // linear for xi,eta and zeta
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
          double dzeta,edge_shape;
          if(dd<3) {
            dzeta = diffN_MBEDGE0;
            edge_shape = N_MBEDGE0(zeta);
          } else {
            dzeta = diffN_MBEDGE1;
            edge_shape = N_MBEDGE1(zeta);
          }
          dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,dd) = tri_n*edge_shape;
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+0) = dksi_tri_n*edge_shape;
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+1) = deta_tri_n*edge_shape;
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(gg,3*dd+2) = tri_n*dzeta;
        }
      }
    }
    // edges on triangles
    for(int ee = 0;ee<9;ee++) {
      if(ee>=3&&ee<=5) {
        // through thickness ho approximation
        // linear xi,eta, ho terms for zeta
        int order = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDataOrder();
        int nb_dofs = NBEDGE_H1_AINSWORTH_COLE(order);
        if((unsigned int)nb_dofs!=dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN().size2()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"nb_dofs != nb_dofs",
            nb_dofs,dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN().size2()
          );
        }
        dataH1.dataOnEntities[MBEDGE][ee].getN().resize(nb_gauss_pts,nb_dofs,false);
        dataH1.dataOnEntities[MBEDGE][ee].getDiffN().resize(nb_gauss_pts,3*nb_dofs,false);
        if(nb_dofs == 0) continue;
        const int prism_edge_map[9][2] = {
          {0,1}, {1,2}, {2,0}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {5,3}
        };
        int gg = 0;
        for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
          double tri_n = dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(ggf,prism_edge_map[ee][0]);
          double dksi_tri_n[2];
          for(int kk = 0;kk<2;kk++) {
            dksi_tri_n[kk] = dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN()(ggf,2*prism_edge_map[ee][0]+kk);
          }
          for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
            double zeta = gaussPtsThroughThickness(0,ggt);
            double n0 = N_MBEDGE0(zeta);
            double n1 = N_MBEDGE1(zeta);
            double n0n1 = n0*n1;
            for(int dd = 0;dd<nb_dofs;dd++) {
              double l = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getN()(ggt,dd);
              double diff_l = dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDiffN()(ggt,dd);
              double edge_m = n0n1*l;
              double dzeta_edge_m = (diffN_MBEDGE0*n1+n0*diffN_MBEDGE1)*l + n0n1*diff_l;
              dataH1.dataOnEntities[MBEDGE][ee].getN()(gg,dd) = tri_n*edge_m;
              for(int kk = 0;kk<2;kk++) {
                dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+kk) = dksi_tri_n[kk]*edge_m;
              }
              dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+2) = tri_n*dzeta_edge_m;
            }
          }
        }
      } else {
        // on triangles ho approximation
        // ho terms on edges, linear zeta
        int nb_dofs = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getN().size2();
        dataH1.dataOnEntities[MBEDGE][ee].getN().resize(nb_gauss_pts,nb_dofs,false);
        dataH1.dataOnEntities[MBEDGE][ee].getDiffN().resize(nb_gauss_pts,3*nb_dofs,false);
        for(int dd = 0;dd<nb_dofs;dd++) {
          int gg = 0;
          for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
            double tri_n = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getN()(ggf,dd);
            double dksi_tri_n = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN()(ggf,2*dd+0);
            double deta_tri_n = dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN()(ggf,2*dd+1);
            for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
              double zeta = gaussPtsThroughThickness(0,ggt);
              double dzeta,edge_shape;
              if(ee<3) {
                dzeta = diffN_MBEDGE0;
                edge_shape = N_MBEDGE0(zeta);
              } else {
                dzeta = diffN_MBEDGE1;
                edge_shape = N_MBEDGE1(zeta);
              }
              dataH1.dataOnEntities[MBEDGE][ee].getN()(gg,dd) = tri_n*edge_shape;
              dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+0) = dksi_tri_n*edge_shape;
              dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+1) = deta_tri_n*edge_shape;
              dataH1.dataOnEntities[MBEDGE][ee].getDiffN()(gg,3*dd+2) = tri_n*dzeta;
            }
          }
        }
      }
    }
    // triangles
    // ho on triangles, linear zeta
    for(int ff = 3;ff<=4;ff++) {
      int nb_dofs;
      nb_dofs = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getN().size2();
      dataH1.dataOnEntities[MBTRI][ff].getN().resize(nb_gauss_pts,nb_dofs,false);
      dataH1.dataOnEntities[MBTRI][ff].getDiffN().resize(nb_gauss_pts,3*nb_dofs,false);
      for(int dd = 0;dd<nb_dofs;dd++) {
        int gg = 0;
        for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
          double tri_n = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getN()(ggf,dd);
          double dksi_tri_n[2];
          for(int kk = 0;kk<2;kk++) {
            dksi_tri_n[kk] = dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getDiffN()(ggf,2*dd+kk);
          }
          for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
            double zeta = gaussPtsThroughThickness(0,ggt);
            double dzeta,edge_shape;
            if(ff == 3) {
              dzeta = diffN_MBEDGE0;
              edge_shape = N_MBEDGE0(zeta);
            } else {
              dzeta = diffN_MBEDGE1;
              edge_shape = N_MBEDGE1(zeta);
            }
            dataH1.dataOnEntities[MBTRI][ff].getN()(gg,dd) = tri_n*edge_shape;
            for(int kk = 0;kk<2;kk++) {
              dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+kk) = dksi_tri_n[kk]*edge_shape;
            }
            dataH1.dataOnEntities[MBTRI][ff].getDiffN()(gg,3*dd+2) = tri_n*dzeta;
          }
        }
      }
    }
    // quads
    // higher order edges and zeta
    {
      int quads_nodes[3*4];
      int quad_order[3];
      double *quad_n[3],*diff_quad_n[3];
      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit;
      siit = side_table.get<1>().lower_bound(boost::make_tuple(MBQUAD,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit;
      hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBQUAD,3));
      EntityHandle ent = fePtr->get_ent();
      rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      // cerr << "\n\n" << endl;
      // const int quad_nodes[3][4] = { {0,1,4,3}, {1,2,5,4}, {0,2,5,3} };
      for(;siit!=hi_siit;siit++) {
        //  cerr << "sn " << siit->side_number << endl;
        int num_nodes_quad;
        const EntityHandle *conn_quad;
        EntityHandle quad = siit->ent;
        rval = mField.get_moab().get_connectivity(
          quad,conn_quad,num_nodes_quad,true
        ); CHKERRQ_MOAB(rval);
        for(int nn = 0;nn<num_nodes_quad;nn++) {
          quads_nodes[4*siit->side_number+nn] = distance(conn,find(conn,conn+6,conn_quad[nn]));
          // cerr
          // << "quad " << quad
          // << " side number " << siit->side_number
          // << " " << quads_nodes[4*siit->side_number+nn]
          // << " " << conn[quads_nodes[4*siit->side_number+nn]]
          // << " " << conn_quad[nn]
          // << endl;
        }
        int order = dataH1.dataOnEntities[MBQUAD][siit->side_number].getDataOrder();
        quad_order[siit->side_number] = order;
        dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().resize(nb_gauss_pts,NBFACEQUAD_H1_AINSWORTH_COLE(order),false);
        dataH1.dataOnEntities[MBQUAD][siit->side_number].getDiffN().resize(nb_gauss_pts,3*NBFACEQUAD_H1_AINSWORTH_COLE(order),false);
        if(dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().size2()>0) {
          quad_n[siit->side_number] = &*dataH1.dataOnEntities[MBQUAD][siit->side_number].getN().data().begin();
          diff_quad_n[siit->side_number] = &*dataH1.dataOnEntities[MBQUAD][siit->side_number].getDiffN().data().begin();
        } else {
          quad_n[siit->side_number] = NULL;
          diff_quad_n[siit->side_number] = NULL;
        }
      }
      {
        double *vertex_n = &*dataH1.dataOnEntities[MBVERTEX][0].getN().data().begin();
        double *diff_vertex_n = &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin();
        ierr = H1_QuadShapeFunctions_MBPRISM(
          quads_nodes,quad_order,vertex_n,diff_vertex_n,quad_n,diff_quad_n,nb_gauss_pts,Legendre_polynomials
        ); CHKERRQ(ierr);
      }
    }
    // prism
    {
      int order = dataH1.dataOnEntities[MBPRISM][0].getDataOrder();
      double *vertex_n  = &dataH1.dataOnEntities[MBVERTEX][0].getN()(0,0);
      double *diff_vertex_n = &dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(0,0);
      dataH1.dataOnEntities[MBPRISM][0].getN().resize(nb_gauss_pts,NBVOLUMEPRISM_H1_AINSWORTH_COLE(order),false);
      dataH1.dataOnEntities[MBPRISM][0].getDiffN().resize(nb_gauss_pts,3*NBVOLUMEPRISM_H1_AINSWORTH_COLE(order),false);
      if(NBVOLUMEPRISM_H1_AINSWORTH_COLE(order)>0) {
        ierr = H1_VolumeShapeFunctions_MBPRISM(
          order,
          vertex_n,
          diff_vertex_n,
          &dataH1.dataOnEntities[MBPRISM][0].getN()(0,0),
          &dataH1.dataOnEntities[MBPRISM][0].getDiffN()(0,0),
          nb_gauss_pts,
          Legendre_polynomials
        ); CHKERRQ(ierr);
      }
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  try {

    EntityHandle ent = fePtr->get_ent();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    coords.resize(num_nodes*3,false);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERRQ_MOAB(rval);

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

    coordsAtGaussPts.resize(nb_gauss_pts,3,false);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPts(gg,dd) = cblas_ddot(
          6,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3
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
        ierr = getEdgesDataOrderSpaceAndBase(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisDataOrderSpaceAndBase(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
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
    FieldApproximationBase base[2];

    boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
    oit = opPtrVector.begin();
    hi_oit = opPtrVector.end();

    for(;oit!=hi_oit;oit++) {

      oit->setPtrFE(this);

      for(int ss = 0;ss!=2;ss++) {

        string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        const MoFEMField* field_struture = mField.get_field_structure(field_name);
        BitFieldId data_id = field_struture->get_id();

        if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data field < %s > on finite element < %s >",
            field_name.c_str(),feName.c_str()
          );
        }

        if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

          space[ss] = field_struture->get_space();
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

          base[ss] = field_struture->get_approx_base();
          switch(base[ss]) {
            case AINSWORTH_COLE_BASE:
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown or not implemented base");
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
              case HCURL:
              if(!ss) {
                ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getEdgesDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              // ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case HDIV:
              if(!ss) {
                ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTrisDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              if(!ss) {
                ierr = getQuadRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getQuadColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getQuadDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getQuadFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              case L2:
              if(!ss) {
                ierr = getPrismRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getPrismColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getPrismDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getPrismFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
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
          ierr = oit->opRhs(
            *op_data[0],
            oit->doVerticesRow,
            oit->doEdgesRow,
            oit->doQuadsRow,
            oit->doTrisRow,
            false,
            oit->doPrismsRow
          ); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPCOL) {
        try {
          ierr = oit->opRhs(
            *op_data[1],
            oit->doVerticesCol,
            oit->doEdgesCol,
            oit->doQuadsCol,
            oit->doTrisCol,
            false,
            oit->doPrismsCol
          ); CHKERRQ(ierr);
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
