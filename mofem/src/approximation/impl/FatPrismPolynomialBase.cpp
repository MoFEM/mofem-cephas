/** \file FatPrismPolynomialBase.cpp
\brief Implementation of Ainsworth-Cole H1 base on edge
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

#include <version.h>
#include <config.h>
#include <definitions.h>
#include <Includes.hpp>

#include <base_functions.h>
#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>
#include <Common.hpp>
#include <UnknownInterface.hpp>
using namespace MoFEM;

#include <FTensor.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <DataStructures.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <LoopMethods.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <FlatPrismPolynomialBase.hpp>
#include <FatPrismPolynomialBase.hpp>

PetscErrorCode FatPrismPolynomialBaseCtx::queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(
    uuid == IDD_FATPRISM_BASE_FUNCTION
  ) {
    *iface = static_cast<FatPrismPolynomialBaseCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = EntPolynomialBaseCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

FatPrismPolynomialBaseCtx::FatPrismPolynomialBaseCtx(
  DataForcesAndSurcesCore &data,
  DataForcesAndSurcesCore &data_triangles_only,
  DataForcesAndSurcesCore &data_trough_thickness,
  ublas::matrix<double>& gauss_pts_triangles_only,
  ublas::matrix<double>& gauss_pts_through_thickness,
  moab::Interface &moab,
  const NumeredEntFiniteElement *fe_ptr,
  const FieldSpace space,
  const FieldApproximationBase base,
  const FieldApproximationBase copy_node_base
):
EntPolynomialBaseCtx(data,space,base,copy_node_base),
dataTrianglesOnly(data_triangles_only),
dataTroughThickness(data_trough_thickness),
gaussPtsTrianglesOnly(gauss_pts_triangles_only),
gaussPtsThroughThickness(gauss_pts_through_thickness),
mOab(moab),
fePtr(fe_ptr) {
  PetscErrorCode ierr;
  ierr = setBase(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}
FatPrismPolynomialBaseCtx::~FatPrismPolynomialBaseCtx() {
}

PetscErrorCode FatPrismPolynomialBase::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_FATPRISM_BASE_FUNCTION) {
    *iface = static_cast<FatPrismPolynomialBase*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

FatPrismPolynomialBase::~FatPrismPolynomialBase() {}
FatPrismPolynomialBase::FatPrismPolynomialBase() {}

PetscErrorCode FatPrismPolynomialBase::getValue(
  ublas::matrix<double> &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_FATPRISM_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  cTx = reinterpret_cast<FatPrismPolynomialBaseCtx*>(iface);
  if(!cTx->fePtr) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
      "Pointer to element should be given "
      "when EntPolynomialBaseCtx is constructed "
      "(use different constructor)"
    );
  }

  int nb_gauss_pts = pts.size2();
  if(!nb_gauss_pts) {
    PetscFunctionReturn(0);
  }

  if(pts.size1()<1) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Wrong dimension of pts, should be at least 3 rows with coordinates"
    );
  }

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSurcesCore& data = cTx->dAta;

  if(cTx->copyNodeBase==LASTBASE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is assumed that base for vertices is calculated");
  } else {
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) = data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
  }
  data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,6,false);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts,12,false);
  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1()!=(unsigned int)nb_gauss_pts) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Base functions or nodes has wrong number of integration points for base %s",
      ApproximationBaseNames[base]
    );
  }
  if(
    cTx->gaussPtsTrianglesOnly.size2()*cTx->gaussPtsThroughThickness.size2()!=
    pts.size2()
  ) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }

  switch (cTx->sPace) {
    case H1:
    ierr = getValueH1TrianglesOnly(); CHKERRQ(ierr);
    ierr = getValueH1ThroughThickness(); CHKERRQ(ierr);
    ierr = getValueH1(pts); CHKERRQ(ierr);
    break;
    case HDIV:
    ierr = getValueHdiv(pts); CHKERRQ(ierr);
    break;
    case HCURL:
    ierr = getValueHCurl(pts); CHKERRQ(ierr);
    break;
    case L2:
    ierr = getValueL2(pts); CHKERRQ(ierr);
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FatPrismPolynomialBase::getValueH1TrianglesOnly() {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  const FieldApproximationBase base = cTx->bAse;
  // PetscErrorCode (*base_polynomials)(
  //   int p,double s,double *diff_s,double *L,double *diffL,const int dim
  // ) = cTx->basePolynomialsType0;

  ierr = FlatPrismPolynomialBase().getValue(
    cTx->gaussPtsTrianglesOnly,
    boost::shared_ptr<BaseFunctionCtx>(
      new FlatPrismPolynomialBaseCtx(
        cTx->dataTrianglesOnly,cTx->mOab,cTx->fePtr,H1,base,NOBASE
      )
    )
  ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode FatPrismPolynomialBase::getValueH1ThroughThickness() {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(
    int p,double s,double *diff_s,double *L,double *diffL,const int dim
  ) = cTx->basePolynomialsType0;

  int nb_gauss_pts_through_thickness = cTx->gaussPtsThroughThickness.size2();

  // Calculate Legendre approx. on edges
  for(unsigned int ee = 3;ee<=5;ee++) {
    int sense = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getSense();
    int order = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getDataOrder()-2;
    cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getN(base).resize(
      nb_gauss_pts_through_thickness,order<0?0:1+order,false
    );
    cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(
      nb_gauss_pts_through_thickness,order<0?0:1+order,false
    );
    if(order<0) continue;
    double diff_s = 2.; // s = s(xi), ds/dxi = 2., because change of basis
    for(int gg = 0;gg<nb_gauss_pts_through_thickness;gg++) {
      double s = 2*cTx->gaussPtsThroughThickness(0,gg)-1; // makes form -1..1
      if(!sense) {
        s *= -1;
        diff_s *= -1;
      }
      // calculate Legendre polynomials at integration points on edges thorough thickness
      ierr = base_polynomials(
        order,s,&diff_s,
        &cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getN(base)(gg,0),
        &cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg,0),
        1
      ); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FatPrismPolynomialBase::getValueH1(ublas::matrix<double> &pts) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  DataForcesAndSurcesCore& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  // PetscErrorCode (*base_polynomials)(
  //   int p,double s,double *diff_s,double *L,double *diffL,const int dim
  // ) = cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();
  int nb_gauss_pts_on_faces = cTx->gaussPtsTrianglesOnly.size2();
  int nb_gauss_pts_through_thickness = cTx->gaussPtsThroughThickness.size2();

  try {
    // nodes
    // linear for xi,eta and zeta
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,6,false);
    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts,18);
    noalias(data.dataOnEntities[MBVERTEX][0].getN(base)) = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);
    noalias(data.dataOnEntities[MBVERTEX][0].getDiffN(base)) = data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);

    // edges on triangles
    for(int ee = 0;ee<9;ee++) {
      if(ee>=3&&ee<=5) {
        // through thickness ho approximation
        // linear xi,eta, ho terms for zeta
        int order = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getDataOrder();
        int nb_dofs = NBEDGE_H1(order);
        if((unsigned int)nb_dofs!=cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getN(base).size2()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"nb_dofs != nb_dofs %d != %d",
            nb_dofs,cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getN(base).size2()
          );
        }
        data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,nb_dofs,false);
        data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
        if(nb_dofs == 0) continue;
        const int prism_edge_map[9][2] = {
          {0,1}, {1,2}, {2,0}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {5,3}
        };
        int gg = 0;
        for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
          double tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBVERTEX][0].getN(base)(ggf,prism_edge_map[ee][0]);
          double dksi_tri_n[2];
          for(int kk = 0;kk<2;kk++) {
            dksi_tri_n[kk] = cTx->dataTrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN(base)(
              ggf,2*prism_edge_map[ee][0]+kk
            );
          }
          for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {

            double zeta = cTx->gaussPtsThroughThickness(0,ggt);
            double n0 = N_MBEDGE0(zeta);
            double n1 = N_MBEDGE1(zeta);
            double n0n1 = n0*n1;

            for(int dd = 0;dd<nb_dofs;dd++) {

              double l = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getN(base)(ggt,dd);
              double diff_l = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee].getDiffN(base)(ggt,dd);

              double edge_m = n0n1*l;
              double dzeta_edge_m = (diffN_MBEDGE0*n1+n0*diffN_MBEDGE1)*l + n0n1*diff_l;
              data.dataOnEntities[MBEDGE][ee].getN(base)(gg,dd) = tri_n*edge_m;
              for(int kk = 0;kk<2;kk++) {
                data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg,3*dd+kk) = dksi_tri_n[kk]*edge_m;
              }
              data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg,3*dd+2) = tri_n*dzeta_edge_m;

            }
          }
        }
      } else {
        // on triangles ho approximation
        // ho terms on edges, linear zeta
        int nb_dofs = cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getN(base).size2();
        data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,nb_dofs,false);
        data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
        for(int dd = 0;dd<nb_dofs;dd++) {
          int gg = 0;
          for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
            double tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getN(base)(ggf,dd);
            double dksi_tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN(base)(ggf,2*dd+0);
            double deta_tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN(base)(ggf,2*dd+1);
            for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
              double zeta = cTx->gaussPtsThroughThickness(0,ggt);
              double dzeta,edge_shape;
              if(ee<3) {
                dzeta = diffN_MBEDGE0;
                edge_shape = N_MBEDGE0(zeta);
              } else {
                dzeta = diffN_MBEDGE1;
                edge_shape = N_MBEDGE1(zeta);
              }
              data.dataOnEntities[MBEDGE][ee].getN(base)(gg,dd) = tri_n*edge_shape;
              data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg,3*dd+0) = dksi_tri_n*edge_shape;
              data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg,3*dd+1) = deta_tri_n*edge_shape;
              data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg,3*dd+2) = tri_n*dzeta;
            }
          }
        }
      }
    }
    // triangles
    // ho on triangles, linear zeta
    for(int ff = 3;ff<=4;ff++) {
      int nb_dofs;
      nb_dofs = cTx->dataTrianglesOnly.dataOnEntities[MBTRI][ff].getN(base).size2();
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,nb_dofs,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,3*nb_dofs,false);
      for(int dd = 0;dd<nb_dofs;dd++) {
        int gg = 0;
        for(int ggf = 0;ggf<nb_gauss_pts_on_faces;ggf++) {
          double tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBTRI][ff].getN(base)(ggf,dd);
          double dksi_tri_n[2];
          for(int kk = 0;kk<2;kk++) {
            dksi_tri_n[kk] = cTx->dataTrianglesOnly.dataOnEntities[MBTRI][ff].getDiffN(base)(ggf,2*dd+kk);
          }
          for(int ggt = 0;ggt<nb_gauss_pts_through_thickness;ggt++,gg++) {
            double zeta = cTx->gaussPtsThroughThickness(0,ggt);
            double dzeta,edge_shape;
            if(ff == 3) {
              dzeta = diffN_MBEDGE0;
              edge_shape = N_MBEDGE0(zeta);
            } else {
              dzeta = diffN_MBEDGE1;
              edge_shape = N_MBEDGE1(zeta);
            }
            data.dataOnEntities[MBTRI][ff].getN(base)(gg,dd) = tri_n*edge_shape;
            for(int kk = 0;kk<2;kk++) {
              data.dataOnEntities[MBTRI][ff].getDiffN(base)(gg,3*dd+kk) = dksi_tri_n[kk]*edge_shape;
            }
            data.dataOnEntities[MBTRI][ff].getDiffN(base)(gg,3*dd+2) = tri_n*dzeta;
          }
        }
      }
    }
    // quads
    // higher order edges and zeta
    {
      MoABErrorCode rval;
      int quads_nodes[3*4];
      int quad_order[3] = { 0, 0, 0};
      double *quad_n[3],*diff_quad_n[3];
      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(cTx->fePtr->getSideNumberTable());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit;
      siit = side_table.get<1>().lower_bound(boost::make_tuple(MBQUAD,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit;
      hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBQUAD,3));
      EntityHandle ent = cTx->fePtr->getEnt();
      int num_nodes;
      const EntityHandle *conn;
      rval = cTx->mOab.get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      // std::cerr << "\n\n" << std::endl;
      // const int quad_nodes[3][4] = { {0,1,4,3}, {1,2,5,4}, {0,2,5,3} };
      for(;siit!=hi_siit;siit++) {
        //  std::cerr << "sn " << siit->side_number << std::endl;
        int num_nodes_quad;
        const EntityHandle *conn_quad;
        EntityHandle quad = siit->get()->ent;
        rval = cTx->mOab.get_connectivity(
          quad,conn_quad,num_nodes_quad,true
        ); CHKERRQ_MOAB(rval);
        for(int nn = 0;nn<num_nodes_quad;nn++) {
          quads_nodes[4*siit->get()->side_number+nn] = std::distance(conn,std::find(conn,conn+6,conn_quad[nn]));
          // std::cerr
          // << "quad " << quad
          // << " side number " << siit->side_number
          // << " " << quads_nodes[4*siit->side_number+nn]
          // << " " << conn[quads_nodes[4*siit->side_number+nn]]
          // << " " << conn_quad[nn]
          // << std::endl;
        }
        int order = data.dataOnEntities[MBQUAD][siit->get()->side_number].getDataOrder();
        quad_order[siit->get()->side_number] = order;
        data.dataOnEntities[MBQUAD][siit->get()->side_number].getN(base).resize(nb_gauss_pts,NBFACEQUAD_H1(order),false);
        data.dataOnEntities[MBQUAD][siit->get()->side_number].getDiffN(base).resize(nb_gauss_pts,3*NBFACEQUAD_H1(order),false);
        if(data.dataOnEntities[MBQUAD][siit->get()->side_number].getN(base).size2()>0) {
          quad_n[siit->get()->side_number] = &*data.dataOnEntities[MBQUAD][siit->get()->side_number].getN(base).data().begin();
          diff_quad_n[siit->get()->side_number] = &*data.dataOnEntities[MBQUAD][siit->get()->side_number].getDiffN(base).data().begin();
        } else {
          quad_n[siit->get()->side_number] = NULL;
          diff_quad_n[siit->get()->side_number] = NULL;
        }
      }
      if(quad_order[0]>0||quad_order[1]>0||quad_order[2]>0) {
        double *vertex_n = &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin();
        double *diff_vertex_n = &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin();
        ierr = H1_QuadShapeFunctions_MBPRISM(
          quads_nodes,quad_order,vertex_n,diff_vertex_n,quad_n,diff_quad_n,nb_gauss_pts,Legendre_polynomials
        ); CHKERRQ(ierr);
      }
    }
    // prism
    {
      int order = data.dataOnEntities[MBPRISM][0].getDataOrder();
      double *vertex_n  = &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0);
      double *diff_vertex_n = &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0);
      data.dataOnEntities[MBPRISM][0].getN(base).resize(nb_gauss_pts,NBVOLUMEPRISM_H1(order),false);
      data.dataOnEntities[MBPRISM][0].getDiffN(base).resize(nb_gauss_pts,3*NBVOLUMEPRISM_H1(order),false);
      if(NBVOLUMEPRISM_H1(order)>0) {
        ierr = H1_VolumeShapeFunctions_MBPRISM(
          order,
          vertex_n,
          diff_vertex_n,
          &data.dataOnEntities[MBPRISM][0].getN(base)(0,0),
          &data.dataOnEntities[MBPRISM][0].getDiffN(base)(0,0),
          nb_gauss_pts,
          Legendre_polynomials
        ); CHKERRQ(ierr);
      }
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FatPrismPolynomialBase::getValueL2(ublas::matrix<double> &pts) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
  PetscFunctionReturn(0);
}

PetscErrorCode FatPrismPolynomialBase::getValueHdiv(ublas::matrix<double> &pts) {
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
}

PetscErrorCode FatPrismPolynomialBase::getValueHCurl(ublas::matrix<double> &pts) {
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
}
