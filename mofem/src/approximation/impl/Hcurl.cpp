/** \file Hcurl.cpp

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

#include <petscsys.h>
#include <FTensor.hpp>
#include <h1_hdiv_hcurl_l2.h>
#include <Hcurl.hpp>
#include <definitions.h>

using namespace MoFEM;

#ifndef GENERATE_VTK_WITH_CURL_BASE

PetscErrorCode MoFEM::Hcurl_EdgeBaseFunctions_MBTET(
  int *sense,int *p,double *N,double *diffN,double *edgeN[],double *diff_edgeN[],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const int edges_nodes[6][2] = {
    {0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}
  };
  int P[6];
  for(int ee = 0;ee<6; ee++)
    P[ee] = NBEDGE_HCURL_AINSWORTH_COLE(p[ee]);

  FTensor::Index<'i',3> i;
  FTensor::Tensor1<double*,3> t_node_diff_ksi[4] = {
    FTensor::Tensor1<double*,3>(&diffN[0],&diffN[ 1],&diffN[ 2]),
    FTensor::Tensor1<double*,3>(&diffN[3],&diffN[ 4],&diffN[ 5]),
    FTensor::Tensor1<double*,3>(&diffN[6],&diffN[ 7],&diffN[ 8]),
    FTensor::Tensor1<double*,3>(&diffN[9],&diffN[10],&diffN[11])
  };
  double edge_diff_ksi[6][3];
  FTensor::Tensor1<double*,3> t_edge_diff_ksi[6] = {
    FTensor::Tensor1<double*,3>(&edge_diff_ksi[0][0],&edge_diff_ksi[0][1],&edge_diff_ksi[0][2]),
    FTensor::Tensor1<double*,3>(&edge_diff_ksi[1][0],&edge_diff_ksi[1][1],&edge_diff_ksi[1][2]),
    FTensor::Tensor1<double*,3>(&edge_diff_ksi[2][0],&edge_diff_ksi[2][1],&edge_diff_ksi[2][2]),
    FTensor::Tensor1<double*,3>(&edge_diff_ksi[3][0],&edge_diff_ksi[3][1],&edge_diff_ksi[3][2]),
    FTensor::Tensor1<double*,3>(&edge_diff_ksi[4][0],&edge_diff_ksi[4][1],&edge_diff_ksi[4][2]),
    FTensor::Tensor1<double*,3>(&edge_diff_ksi[5][0],&edge_diff_ksi[5][1],&edge_diff_ksi[5][2])
  };
  for(int ee = 0;ee!=6;ee++) {
    t_edge_diff_ksi[ee](i) = (
      t_node_diff_ksi[edges_nodes[ee][1]](i)-t_node_diff_ksi[edges_nodes[ee][0]](i)
    )*sense[ee];
  }

  FTensor::Tensor1<double*,3> t_edge_n[6] = {
    FTensor::Tensor1<double*,3>(&edgeN[0][0],&edgeN[0][1],&edgeN[0][2],3),
    FTensor::Tensor1<double*,3>(&edgeN[1][0],&edgeN[1][1],&edgeN[1][2],3),
    FTensor::Tensor1<double*,3>(&edgeN[2][0],&edgeN[2][1],&edgeN[2][2],3),
    FTensor::Tensor1<double*,3>(&edgeN[3][0],&edgeN[3][1],&edgeN[3][2],3),
    FTensor::Tensor1<double*,3>(&edgeN[4][0],&edgeN[4][1],&edgeN[4][2],3),
    FTensor::Tensor1<double*,3>(&edgeN[5][0],&edgeN[5][1],&edgeN[5][2],3)
  };
  FTensor::Tensor1<double,3> t_psi_e_0,t_psi_e_1;

  for(int ii = 0;ii!=GDIM;ii++) {

    const int node_shift = ii*4;
    for(int ee = 0;ee!=6;ee++) {

      t_psi_e_0(i) =
      N[node_shift+edges_nodes[ee][1]]*t_node_diff_ksi[edges_nodes[ee][0]](i)-
      N[node_shift+edges_nodes[ee][0]]*t_node_diff_ksi[edges_nodes[ee][1]](i);
      t_psi_e_1(i) =
      N[node_shift+edges_nodes[ee][1]]*t_node_diff_ksi[edges_nodes[ee][0]](i)+
      N[node_shift+edges_nodes[ee][0]]*t_node_diff_ksi[edges_nodes[ee][1]](i);
      (t_edge_n[ee])(i) = t_psi_e_0(i);
      ++(t_edge_n[ee]);
      (t_edge_n[ee])(i) = t_psi_e_1(i);
      ++(t_edge_n[ee]);

      if(p[ee]>1) {
        const double ksi_0i = N[node_shift+edges_nodes[ee][1]]-N[node_shift+edges_nodes[ee][0]];
        double psi_l[p[ee]+1],diff_psi_l[3*p[ee]+3];
        ierr = base_polynomials(
          p[ee],ksi_0i,&edge_diff_ksi[ee][0],psi_l,diff_psi_l,3
        ); CHKERRQ(ierr);
        for(int ll = 2;ll!=P[ee];ll++) {
          const double a = (double)(2*ll+1)/(double)(ll+1);
          const double b = (double)(ll)/(double)(ll+1);
          (t_edge_n[ee])(i) = a*psi_l[ll-1]*t_psi_e_1(i)-b*psi_l[ll-2]*t_psi_e_0(i);
          ++(t_edge_n[ee]);
        }
      }

    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode MoFEM::Hcurl_EdgeBasedFaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f_e[4][3],double *diff_phi_f_e[4][3],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  const int edges[3][2] = { {0,1}, {1,2}, {2,0} };

  FTensor::Index<'i',3> i;
  FTensor::Tensor1<double*,3> t_node_diff_ksi[4] = {
    FTensor::Tensor1<double*,3>(&diffN[0],&diffN[ 1],&diffN[ 2]),
    FTensor::Tensor1<double*,3>(&diffN[3],&diffN[ 4],&diffN[ 5]),
    FTensor::Tensor1<double*,3>(&diffN[6],&diffN[ 7],&diffN[ 8]),
    FTensor::Tensor1<double*,3>(&diffN[9],&diffN[10],&diffN[11])
  };
  FTensor::Tensor1<double,3> t_edge_diff_ksi;

  for(int ff = 0;ff!=4;ff++) {

    const int o_nodes[3] = {
      faces_nodes[3*ff+2],faces_nodes[3*ff+0],faces_nodes[3*ff+1]
    };
    FTensor::Tensor1<double*,3> t_o_diff[3] = {
      FTensor::Tensor1<double*,3>(&diffN[3*o_nodes[0]+0],&diffN[3*o_nodes[0]+1],&diffN[3*o_nodes[0]+2]),
      FTensor::Tensor1<double*,3>(&diffN[3*o_nodes[1]+0],&diffN[3*o_nodes[1]+1],&diffN[3*o_nodes[1]+2]),
      FTensor::Tensor1<double*,3>(&diffN[3*o_nodes[2]+0],&diffN[3*o_nodes[2]+1],&diffN[3*o_nodes[2]+2])
    };
    double psi_l[p[ff]+1],diff_psi_l[3*p[ff]+3];

    const int nb_base_fun_on_face = NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(p[ff]);

    for(int ee = 0;ee!=3;ee++) {

      FTensor::Tensor1<double*,3> t_face_edge_base(
        &phi_f_e[ff][ee][0],&phi_f_e[ff][ee][1],&phi_f_e[ff][ee][2],3
      );

      t_edge_diff_ksi(i) =
      t_node_diff_ksi[o_nodes[1]](i)-t_node_diff_ksi[o_nodes[0]](i);

      for(int ii = 0;ii!=GDIM;ii++) {

        const int node_shift = ii*4;
        const double n[] = {
          N[node_shift+faces_nodes[3*ff+0]],
          N[node_shift+faces_nodes[3*ff+1]],
          N[node_shift+faces_nodes[3*ff+2]]
        };
        const double ksi_0i = n[edges[ee][1]]-n[edges[ee][0]];
        ierr = base_polynomials(
          p[ff],ksi_0i,&t_edge_diff_ksi(0),psi_l,diff_psi_l,3
        ); CHKERRQ(ierr);
        const double beta_e = n[edges[ee][0]]*n[edges[ee][1]];

        for(int ll = 0;ll!=nb_base_fun_on_face;ll++) {
          t_face_edge_base(i) = beta_e*psi_l[ll]*t_o_diff[ee](i);
          ++t_face_edge_base;
        }

      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode MoFEM::Hcurl_BubbleFaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f[4],double *diff_phi_f[4],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  const double coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };

  FTensor::Index<'i',3> i;
  FTensor::Tensor1<double*,3> t_node_diff_ksi[4] = {
    FTensor::Tensor1<double*,3>(&diffN[0],&diffN[ 1],&diffN[ 2]),
    FTensor::Tensor1<double*,3>(&diffN[3],&diffN[ 4],&diffN[ 5]),
    FTensor::Tensor1<double*,3>(&diffN[6],&diffN[ 7],&diffN[ 8]),
    FTensor::Tensor1<double*,3>(&diffN[9],&diffN[10],&diffN[11])
  };
  FTensor::Tensor1<double,3> t_diff_ksi0i,t_diff_ksi0j;

  for(int ff = 0;ff!=4;ff++) {

    FTensor::Tensor1<double,3> t_coords0(
      3*coords[faces_nodes[3*ff+0]+0],
      3*coords[faces_nodes[3*ff+0]+1],
      3*coords[faces_nodes[3*ff+0]+2]
    );
    FTensor::Tensor1<double,3> tou_0i(
      coords[3*faces_nodes[3*ff+1]+0],
      coords[3*faces_nodes[3*ff+1]+1],
      coords[3*faces_nodes[3*ff+1]+2]
    );
    FTensor::Tensor1<double,3> tou_0j(
      coords[3*faces_nodes[3*ff+2]+0],
      coords[3*faces_nodes[3*ff+2]+1],
      coords[3*faces_nodes[3*ff+2]+2]
    );
    tou_0i(i) -= t_coords0(i);
    tou_0j(i) -= t_coords0(i);

    t_diff_ksi0i(i) =
    t_node_diff_ksi[faces_nodes[3*ff+1]](i)-t_node_diff_ksi[faces_nodes[3*ff+0]](i);
    t_diff_ksi0j(i) =
    t_node_diff_ksi[faces_nodes[3*ff+2]](i)-t_node_diff_ksi[faces_nodes[3*ff+0]](i);

    double psi_l_0i[p[ff]+1],diff_psi_l_0i[3*p[ff]+3];
    double psi_l_0j[p[ff]+1],diff_psi_l_0j[3*p[ff]+3];

    FTensor::Tensor1<double*,3> t_phi_f(
      &phi_f[ff][0],&phi_f[ff][1],&phi_f[ff][2],3
    );

    for(int ii = 0;ii!=GDIM;ii++) {

      const int node_shift = ii*4;
      const double beta_0ij =
      N[node_shift+faces_nodes[3*ff+0]]*
      N[node_shift+faces_nodes[3*ff+1]]*
      N[node_shift+faces_nodes[3*ff+2]];

      const double ksi_0i =
      N[node_shift+faces_nodes[3*ff+1]]-N[node_shift+faces_nodes[3*ff+0]];
      ierr = base_polynomials(
        p[ff],ksi_0i,&t_diff_ksi0i(0),psi_l_0i,diff_psi_l_0i,3
      ); CHKERRQ(ierr);

      const double ksi_0j =
      N[node_shift+faces_nodes[3*ff+2]]-N[node_shift+faces_nodes[3*ff+0]];
      ierr = base_polynomials(
        p[ff],ksi_0j,&t_diff_ksi0j(0),psi_l_0j,diff_psi_l_0j,3
      ); CHKERRQ(ierr);

      int cc = 0;
      for(int oo = 0;oo<=(p[ff]-3);oo++) {
        for(int pp0 = 0;pp0<=oo;pp0++) {
          const int pp1 = oo-pp0;
          if(pp1>=0) {
            const double a = beta_0ij*psi_l_0i[pp0]*psi_l_0j[pp1];
            t_phi_f(i) = a*tou_0i(i);
            ++t_phi_f;
            ++cc;
            t_phi_f(i) = a*tou_0j(i);
            ++t_phi_f;
            ++cc;
          }
        }
      }

      const int nb_base_fun_on_face = NBFACETRI_FACE_HCURL_AINSWORTH_COLE(p[ff]);
      if(cc!=nb_base_fun_on_face) {
        SETERRQ2(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "Wrong number of base functions %d != %d",
          cc,nb_base_fun_on_face
        );
      }

    }

  }
  PetscFunctionReturn(0);
}

#endif // Not GENERATE_VTK_WITH_CURL_BASE

#ifdef GENERATE_VTK_WITH_CURL_BASE

#include <MoFEM.hpp>
using namespace MoFEM;
using namespace boost::numeric;

PetscErrorCode VTK_Hcurl_MBTET(const string file_name) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  double base_coords[] = {
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1
  };

  moab::Core core_ref;
  moab::Interface& moab_ref = core_ref;

  EntityHandle nodes[4];
  for(int nn = 0;nn<4;nn++) {
    rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERRQ_MOAB(rval);
  }
  EntityHandle tet;
  rval = moab_ref.create_element(MBTET,nodes,4,tet); CHKERRQ_MOAB(rval);

  MoFEM::Core m_core_ref(moab_ref,PETSC_COMM_SELF,-2);
  MoFEM::Interface& m_field_ref = m_core_ref;

  ierr = m_field_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

  const int max_level = 4;
  for(int ll = 0;ll!=max_level;ll++) {
    Range edges;
    ierr = m_field_ref.get_entities_by_type_and_ref_level
    (BitRefLevel().set(ll),BitRefLevel().set(),MBEDGE,edges); CHKERRQ(ierr);
    Range tets;
    ierr = m_field_ref.get_entities_by_type_and_ref_level
    (BitRefLevel().set(ll),BitRefLevel(ll).set(),MBTET,tets); CHKERRQ(ierr);
    //refine mesh
    MeshRefinment& m_ref = m_core_ref;
    ierr = m_ref.add_verices_in_the_middel_of_edges(edges,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
    ierr = m_ref.refine_TET(tets,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
  }

  Range tets;
  ierr = m_field_ref.get_entities_by_type_and_ref_level(
    BitRefLevel().set(max_level),BitRefLevel().set(max_level),MBTET,tets
  ); CHKERRQ(ierr);

  // Use 10 node tets to print base
  if(1) {

    // Range edges;
    // rval = moab_ref.get_adjacencies(tets,1,true,edges); CHKERRQ_MOAB(rval);
    EntityHandle meshset;
    rval = moab_ref.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
    rval = moab_ref.add_entities(meshset,tets); CHKERRQ_MOAB(rval);
    rval = moab_ref.convert_entities(meshset,true,false,false); CHKERRQ_MOAB(rval);
    rval = moab_ref.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);

  }

  Range elem_nodes;
  rval = moab_ref.get_connectivity(tets,elem_nodes,false); CHKERRQ_MOAB(rval);

  const int nb_gauss_pts = elem_nodes.size();
  MatrixDouble gauss_pts(nb_gauss_pts,4);
  gauss_pts.clear();
  Range::iterator nit = elem_nodes.begin();
  for(int gg = 0;nit!=elem_nodes.end();nit++,gg++) {
    rval = moab_ref.get_coords(&*nit,1,&gauss_pts(gg,0)); CHKERRQ_MOAB(rval);
  }
  gauss_pts = trans(gauss_pts);

  MatrixDouble shape_fun;
  shape_fun.resize(nb_gauss_pts,4);
  ierr = ShapeMBTET(
    &*shape_fun.data().begin(),&gauss_pts(0,0),&gauss_pts(1,0),&gauss_pts(2,0),nb_gauss_pts
  ); CHKERRQ(ierr);

  double diff_shape_fun[12];
  ierr = ShapeDiffMBTET(diff_shape_fun); CHKERRQ(ierr);

  int edge_sense[6] = { 1,1,1, 1,1,1 };
  const int order = 4;
  int edge_order[6] = { order,order,order, order,order,order };

  MatrixDouble base_edge_functions(
    6,3*nb_gauss_pts*NBEDGE_HCURL_AINSWORTH_COLE(order)
  );
  double* edge_n[] = {
    &base_edge_functions(0,0),
    &base_edge_functions(1,0),
    &base_edge_functions(2,0),
    &base_edge_functions(3,0),
    &base_edge_functions(4,0),
    &base_edge_functions(5,0)
  };

  MatrixDouble diff_base_edge_functions(
    6,9*nb_gauss_pts*NBEDGE_HCURL_AINSWORTH_COLE(order)
  );
  double* diff_edge_n[] = {
    &diff_base_edge_functions(0,0),
    &diff_base_edge_functions(1,0),
    &diff_base_edge_functions(2,0),
    &diff_base_edge_functions(3,0),
    &diff_base_edge_functions(4,0),
    &diff_base_edge_functions(5,0)
  };

  ierr = Hcurl_EdgeBaseFunctions_MBTET(
    edge_sense,
    edge_order,
    &*shape_fun.data().begin(),
    diff_shape_fun,
    edge_n,
    diff_edge_n,
    nb_gauss_pts,
    Legendre_polynomials
  );

  double def_val[] = { 0,0,0 };

  for(int  ee = 0;ee!=6;ee++) {
    for(int ll = 0;ll!=NBEDGE_HCURL_AINSWORTH_COLE(order);ll++) {
      std::ostringstream ss;
      ss << "curl_edge_" << ee << "_" << ll;
      Tag th;
      rval = moab_ref.tag_get_handle(
        ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
      ); CHKERRQ_MOAB(rval);

      int gg = 0;
      for(Range::iterator nit = elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
        rval = moab_ref.tag_set_data(
          th,&*nit,1,&(edge_n[ee][gg*3*NBEDGE_HCURL_AINSWORTH_COLE(order)+ll*3])
        ); CHKERRQ_MOAB(rval);
      }
    }
  }

  int faces_order[] = { order,order,order,order };
  int faces_nodes[] = { 0,1,3, 1,2,3, 0,2,3, 0,1,2 };
  MatrixDouble base_face_edge_functions(
    4,3*3*NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(order)*nb_gauss_pts
  );
  MatrixDouble diff_base_face_edge_functions(
    4,3*9*NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(order)*nb_gauss_pts
  );
  double *phi_f_e[4][3];
  double *diff_phi_f_e[4][3];
  for(int ff = 0;ff!=4;ff++) {
    for(int ee = 0;ee!=3;ee++) {
      phi_f_e[ff][ee] = &base_face_edge_functions(ff,ee*3*NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(order)*nb_gauss_pts);
      diff_phi_f_e[ff][ee] = &diff_base_face_edge_functions(ff,ee*9*NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(order)*nb_gauss_pts);
    }
  }

  ierr = Hcurl_EdgeBasedFaceFunctions_MBTET(
    faces_nodes,
    faces_order,
    &*shape_fun.data().begin(),
    diff_shape_fun,
    phi_f_e,
    diff_phi_f_e,
    nb_gauss_pts,
    Legendre_polynomials
  ); CHKERRQ(ierr);

  for(int ff = 0;ff!=4;ff++) {
    for(int  ee = 0;ee!=3;ee++) {
      for(int ll = 0;ll!=NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(order);ll++) {
        std::ostringstream ss;
        ss << "curl_face_edge_" << ff << "_" << ee << "_" << ll;
        Tag th;
        rval = moab_ref.tag_get_handle(
          ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
        ); CHKERRQ_MOAB(rval);

        int gg = 0;
        for(Range::iterator nit = elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {

          int idx =
          3*NBFACETRI_EDGE_HCURL_AINSWORTH_COLE(order)*gg+ll*3;
          if(idx >= base_face_edge_functions.size2()) {
            cerr << ff << " " << ee << " " << ll << " " << gg << endl;
          }

          rval = moab_ref.tag_set_data(th,&*nit,1,&(phi_f_e[ff][ee][idx])); CHKERRQ_MOAB(rval);
        }
      }
    }
  }

  MatrixDouble base_face_bubble_functions(
    4,3*NBFACETRI_FACE_HCURL_AINSWORTH_COLE(order)*nb_gauss_pts
  );
  MatrixDouble diff_base_face_bubble_functions(
    4,9*NBFACETRI_FACE_HCURL_AINSWORTH_COLE(order)*nb_gauss_pts
  );
  double *phi_f[4];
  double *diff_phi_f[4];
  for(int ff=0;ff!=4;ff++) {
    phi_f[ff] = &base_face_bubble_functions(ff,0);
    diff_phi_f[ff] = &diff_base_face_bubble_functions(ff,0);
  }

  ierr = Hcurl_BubbleFaceFunctions_MBTET(
    faces_nodes,
    faces_order,
    &*shape_fun.data().begin(),
    diff_shape_fun,
    phi_f,
    diff_phi_f,
    nb_gauss_pts,
    Legendre_polynomials
  ); CHKERRQ(ierr);

  for(int ff = 0;ff!=4;ff++) {
    for(int ll = 0;ll!=NBFACETRI_FACE_HCURL_AINSWORTH_COLE(order);ll++) {
      std::ostringstream ss;
      ss << "curl_face_bubble_" << ff << "_" << ll;
      Tag th;
      rval = moab_ref.tag_get_handle(
        ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
      ); CHKERRQ_MOAB(rval);

      int gg = 0;
      for(Range::iterator nit = elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
        int idx = 3*NBFACETRI_FACE_HCURL_AINSWORTH_COLE(order)*gg+ll*3;
        rval = moab_ref.tag_set_data(th,&*nit,1,&(phi_f[ff][idx])); CHKERRQ_MOAB(rval);
      }
    }
  }

  EntityHandle meshset;
  rval = moab_ref.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
  rval = moab_ref.add_entities(meshset,tets); CHKERRQ_MOAB(rval);
  rval = moab_ref.write_file(file_name.c_str(),"VTK","",&meshset,1); CHKERRQ_MOAB(rval);

  PetscFunctionReturn(0);
}


static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscErrorCode ierr;
  ierr = VTK_Hcurl_MBTET("out_curl_vtk_base_on_tet.vtk"); CHKERRQ(ierr);

  PetscFinalize();

  return 0;
}

#endif // GENERATE_VTK_WITH_CURL_BASE
