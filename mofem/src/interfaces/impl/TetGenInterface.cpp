/** \file TetGenInterface.cpp
 * \brief TetGen interface for meshing/re-meshing and on the fly mesh creation
 *
 */

/**
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

#ifdef WITH_TETGEN
  #include <tetgen.h>
#endif

#ifdef WITH_TETGEN

// #define DEBUG_TETGEN
#ifdef DEBUG_TETGEN
  #include <FTensor.hpp>

  template<class T1,class T2>
  static inline MoFEMErrorCode determinantTensor3by3(
    T1 &t,T2 &det
  ) {
    MoFEMFunctionBeginHot;
    det =
    +t(0,0)*t(1,1)*t(2,2) + t(1,0)*t(2,1)*t(0,2)
    +t(2,0)*t(0,1)*t(1,2) - t(0,0)*t(2,1)*t(1,2)
    -t(2,0)*t(1,1)*t(0,2) - t(1,0)*t(0,1)*t(2,2);
    MoFEMFunctionReturnHot(0);
  }
#endif

  namespace MoFEM {

  MoFEMErrorCode
  TetGenInterface::query_interface(const MOFEMuuid &uuid,
                                   UnknownInterface **iface) const {
    MoFEMFunctionBeginHot;
    *iface = NULL;
    if (uuid == IDD_MOFEMTetGegInterface) {
      *iface = const_cast<TetGenInterface *>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");

    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode TetGenInterface::inData(
      Range &ents, tetgenio &in,
      std::map<EntityHandle, unsigned long> &moab_tetgen_map,
      std::map<unsigned long, EntityHandle> &tetgen_moab_map, Tag th) {
    MoFEMFunctionBeginHot;

    Interface &m_field = cOre;
    Range::iterator it;

    //
    ErrorCode rval;

    Tag th_marker;
    int def_marker = 0;
    rval = m_field.get_moab().tag_get_handle(
        "TETGEN_MARKER", 1, MB_TYPE_INTEGER, th_marker,
        MB_TAG_CREAT | MB_TAG_SPARSE, &def_marker);
    CHKERRQ_MOAB(rval);

    // All indices start from 0
    in.firstnumber = 0;

    Range points = ents.subset_by_dimension(0);
    in.numberofpoints = points.size();
    if (points.size() > 0) {
      in.pointlist = new double[in.numberofpoints * 3];
      in.pointmarkerlist = new int[in.numberofpoints];
      if (th) {
        rval = m_field.get_moab().tag_get_data(th, points, in.pointlist);
        CHKERRQ_MOAB(rval);
      } else {
        rval = m_field.get_moab().get_coords(points, in.pointlist);
        CHKERRQ_MOAB(rval);
      }
      rval = m_field.get_moab().tag_get_data(th_marker, points,
                                             in.pointmarkerlist);
      CHKERRQ_MOAB(rval);
      it = points.begin();
      for (int ii = 0; it != points.end(); it++, ii++) {
        unsigned long iii = MBVERTEX | (ii << MB_TYPE_WIDTH);
        tetgen_moab_map[iii] = *it;
        moab_tetgen_map[*it] = iii;
      }
    }

    Range tets = ents.subset_by_type(MBTET);
    in.numberoftetrahedra = tets.size();
    if (in.numberoftetrahedra > 0) {
      in.tetrahedronlist = new int[4 * ents.subset_by_type(MBTET).size()];
      it = tets.begin();
      for (int ii = 0; it != tets.end(); it++, ii++) {
        int num_nodes;
        const EntityHandle *conn;
        rval = m_field.get_moab().get_connectivity(*it, conn, num_nodes, true);
        CHKERRQ_MOAB(rval);
#ifdef DEBUG_TETGEN
        {
          double coords[12];
          if (th) {
            rval = m_field.tag_get_data(th, conn, num_nodes, coords);
            CHKERRQ_MOAB(rval);
          } else {
            rval = m_field.get_moab().get_coords(conn, num_nodes, coords);
            CHKERRQ_MOAB(rval);
          }
          double diff_n[12];
          ShapeDiffMBTET(diff_n);
          FTensor::Tensor1<double *, 3> t_diff_n(&diff_n[0], &diff_n[1],
                                                 &diff_n[2], 3);
          FTensor::Tensor1<double *, 3> t_coords(&coords[0], &coords[1],
                                                 &coords[2], 3);
          FTensor::Tensor2<double, 3, 3> jac;
          FTensor::Index<'i', 3> i;
          FTensor::Index<'j', 3> j;
          jac(i, j) = 0;
          for (int nn = 0; nn != 4; nn++) {
            jac(i, j) += t_coords(i) * t_diff_n(j);
            ++t_coords;
            ++t_diff_n;
          }
          double det;
          determinantTensor3by3(jac, det);
          if (det <= 0) {
            SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Negative volume det = %6.4e", det);
          }
        }
#endif
        tetgen_moab_map[MBTET | (ii << MB_TYPE_WIDTH)] = *it;
        moab_tetgen_map[*it] = MBTET | (ii << MB_TYPE_WIDTH);
        for (int nn = 0; nn != 4; nn++) {
          if (moab_tetgen_map.find(conn[nn]) == moab_tetgen_map.end()) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency between TetGen and MoAB");
          }
          in.tetrahedronlist[4 * ii + nn] =
              moab_tetgen_map[conn[nn]] >> MB_TYPE_WIDTH;
        }
      }
    }

    Range tris = ents.subset_by_type(MBTRI);
    in.numberoftrifaces = tris.size();
    if (in.numberoftrifaces) {
      in.trifacelist = new int[3 * in.numberoftrifaces];
      in.trifacemarkerlist = new int[in.numberoftrifaces];
      // std::fill(&in.trifacemarkerlist[0],&in.trifacemarkerlist[in.numberoftrifaces],1);
      rval = m_field.get_moab().tag_get_data(th_marker, tris,
                                             in.trifacemarkerlist);
      CHKERRQ_MOAB(rval);
      it = tris.begin();
      for (int ii = 0; it != tris.end(); it++, ii++) {
        int num_nodes;
        const EntityHandle *conn;
        rval = m_field.get_moab().get_connectivity(*it, conn, num_nodes, true);
        CHKERRQ_MOAB(rval);
        int order[] = {0, 1, 2};
        Range adj_tets;
        rval = m_field.get_moab().get_adjacencies(&*it, 1, 3, true, adj_tets);
        CHKERRQ_MOAB(rval);
        adj_tets = intersect(adj_tets, tets);
        if (adj_tets.size() == 1) {
          int side_number;
          int sense;
          int offset;
          rval = m_field.get_moab().side_number(adj_tets[0], *it, side_number,
                                                sense, offset);
          CHKERRQ_MOAB(rval);
          if (sense == -1) {
            order[0] = 1;
            order[1] = 0;
          }
        }
        tetgen_moab_map[MBTRI | (ii << MB_TYPE_WIDTH)] = *it;
        moab_tetgen_map[*it] = MBTRI | (ii << MB_TYPE_WIDTH);
        for (int nn = 0; nn < 3; nn++) {
          if (moab_tetgen_map.find(conn[order[nn]]) == moab_tetgen_map.end()) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency between TetGen and MoAB");
          }
          in.trifacelist[3 * ii + nn] =
              moab_tetgen_map[conn[order[nn]]] >> MB_TYPE_WIDTH;
        }
      }
    }

    Range edges = ents.subset_by_type(MBEDGE);
    in.numberofedges = edges.size();
    if (in.numberofedges > 0) {
      in.edgelist = new int[2 * in.numberofedges];
      in.edgemarkerlist = new int[in.numberofedges];
      // std::fill(&in.edgemarkerlist[0],&in.edgemarkerlist[in.numberofedges],1);
      rval =
          m_field.get_moab().tag_get_data(th_marker, edges, in.edgemarkerlist);
      CHKERRQ_MOAB(rval);
      it = edges.begin();
      for (int ii = 0; it != edges.end(); it++, ii++) {
        int num_nodes;
        const EntityHandle *conn;
        rval = m_field.get_moab().get_connectivity(*it, conn, num_nodes, true);
        CHKERRQ_MOAB(rval);
        tetgen_moab_map[MBEDGE | (ii << MB_TYPE_WIDTH)] = *it;
        moab_tetgen_map[*it] = MBEDGE | (ii << MB_TYPE_WIDTH);
        for (int nn = 0; nn < 2; nn++) {
          if (moab_tetgen_map.find(conn[nn]) == moab_tetgen_map.end()) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency between TetGen and MoAB");
          }
          in.edgelist[2 * ii + nn] = moab_tetgen_map[conn[nn]] >> MB_TYPE_WIDTH;
        }
      }
    }

    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode TetGenInterface::setGeomData(tetgenio &in,
                                              moabTetGen_Map &moab_tetgen_map,
                                              tetGenMoab_Map &tetgen_moab_map,
                                              std::map<int, Range> &type_ents) {
    MoFEMFunctionBeginHot;

    Interface &m_field = cOre;
    //
    // ErrorCode rval;
    in.pointparamlist = new tetgenio::pointparam[in.numberofpoints];
    // std::vector<bool> points_is_set(in.numberofpoints,false);
    std::map<int, Range>::iterator mit = type_ents.begin();
    for (; mit != type_ents.end(); mit++) {
      if (mit->first < 0 && mit->first > 3) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong TetGen point type");
      }
      Range::iterator it = mit->second.begin();
      for (; it != mit->second.end(); it++) {
        moabTetGen_Map::iterator miit = moab_tetgen_map.find(*it);
        if (miit == moab_tetgen_map.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Data inconsistency between TetGen and MoAB");
          continue;
        }
        int id = miit->second >> MB_TYPE_WIDTH;
        in.pointparamlist[id].uv[0] = 0;
        in.pointparamlist[id].uv[1] = 0;
        in.pointparamlist[id].type = mit->first;
        in.pointparamlist[id].tag = m_field.get_moab().id_from_handle(*it) + 1;
        // points_is_set[id] = true;
      }
    }

    // // Check only if tag and type is set to all points
    // for(
    //   std::vector<bool>::iterator bit = points_is_set.begin();
    //   bit!=points_is_set.end();bit++
    // ) {
    //   if(!*bit) {
    //     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Point type for
    //     TetGen is not set");
    //   }
    // }

    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode TetGenInterface::outData(
      tetgenio &in, tetgenio &out,
      std::map<EntityHandle, unsigned long> &moab_tetgen_map,
      std::map<unsigned long, EntityHandle> &tetgen_moab_map, Range *ents,
      bool id_in_tags, bool error_if_created) {
    MoFEMFunctionBeginHot;

    Interface &m_field = cOre;

    //
    ErrorCode rval;

    Tag th_marker;
    int def_marker = 0;
    rval = m_field.get_moab().tag_get_handle(
        "TETGEN_MARKER", 1, MB_TYPE_INTEGER, th_marker,
        MB_TAG_CREAT | MB_TAG_SPARSE, &def_marker);
    CHKERRQ_MOAB(rval);

    int num_nodes = 0;
    // std::vector<int> new_node_markes;
    // new_node_markes.reserve(out.numberofpoints);
    // std::vector<int> new_node_ii;
    // new_node_ii.reserve(out.numberofpoints);

    int ii = 0;
    for (; ii < out.numberofpoints; ii++) {
      if (ii < in.numberofpoints) {
        if (memcmp(&in.pointlist[3 * ii], &out.pointlist[3 * ii],
                   3 * sizeof(double)) == 0) {
          unsigned long iii = MBVERTEX | (ii << MB_TYPE_WIDTH);
          if (tetgen_moab_map.find(iii) != tetgen_moab_map.end()) {
            rval = m_field.get_moab().tag_set_data(
                th_marker, &tetgen_moab_map[iii], 1, &out.pointmarkerlist[ii]);
            CHKERRQ_MOAB(rval);
            continue;
          } else {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency between TetGen and MoAB");
          }
        }
      }
      if (id_in_tags) {
        if (out.pointparamlist[ii].tag > 0) {
          EntityHandle node;
          rval = m_field.get_moab().handle_from_id(
              MBVERTEX, in.pointparamlist[ii].tag - 1, node);
          CHKERRQ_MOAB(rval);
          if (moab_tetgen_map.find(node) != moab_tetgen_map.end()) {
            continue;
          }
        }
      }
      if (error_if_created) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "node should not be created");
      }
      num_nodes++;
      // new_node_ii.push_back(ii);
      // new_node_markes.push_back(out.pointmarkerlist[ii]);
      // EntityHandle node;
      // rval = m_field.get_moab().create_vertex(&out.pointlist[3*ii],node);
      // CHKERRQ_MOAB(rval);
      // rval =
      // m_field.get_moab().tag_set_data(th_marker,&node,1,&out.pointmarkerlist[ii]);
      // CHKERRQ_MOAB(rval);
      // moab_tetgen_map[node] = MBVERTEX|(ii<<MB_TYPE_WIDTH);
      // tetgen_moab_map[MBVERTEX|(ii<<MB_TYPE_WIDTH)] = node;
      // if(ents!=NULL) ents->insert(node);
    }

    ReadUtilIface *iface;
    rval = m_field.get_moab().query_interface(iface);
    CHKERRQ_MOAB(rval);

    if (num_nodes) {
      vector<double *> arrays;
      EntityHandle startv;
      rval = iface->get_node_coords(3, num_nodes, 0, startv, arrays);
      CHKERRQ_MOAB(rval);
      Range verts(startv, startv + num_nodes - 1);
      int ii = in.numberofpoints;
      for (Range::iterator vit = verts.begin(); vit != verts.end();
           vit++, ii++) {
        for (int nn = 0; nn != 3; nn++) {
          arrays[nn][ii - in.numberofpoints] = out.pointlist[3 * ii + nn];
        }
        moab_tetgen_map[*vit] = MBVERTEX | (ii << MB_TYPE_WIDTH);
        tetgen_moab_map[MBVERTEX | (ii << MB_TYPE_WIDTH)] = *vit;
        if (ents != NULL)
          ents->insert(*vit);
      }
      rval = m_field.get_moab().tag_set_data(
          th_marker, verts, &out.pointmarkerlist[in.numberofpoints]);
      CHKERRQ_MOAB(rval);
    }

    std::vector<int> tetgen_ii;

    // Build tets
    std::vector<EntityHandle> conn_seq_tet;
    conn_seq_tet.reserve(4 * out.numberoftetrahedra);
    tetgen_ii.reserve(out.numberoftetrahedra);
    conn_seq_tet.clear();
    tetgen_ii.clear();
    ii = 0;
    for (; ii < out.numberoftetrahedra; ii++) {
      unsigned long iii = MBTET | (ii << MB_TYPE_WIDTH);
      if (ii < in.numberoftetrahedra) {
        if (memcmp(&in.tetrahedronlist[4 * ii], &out.tetrahedronlist[4 * ii],
                   4 * sizeof(int)) == 0) {
          if (tetgen_moab_map.find(iii) != tetgen_moab_map.end()) {
            if (ents != NULL)
              ents->insert(tetgen_moab_map[iii]);
            continue;
          } else {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "data inconsistency between TetGen and MoAB");
          }
        }
      }
      EntityHandle conn[4];
      for (int nn = 0; nn < 4; nn++) {
        int nnn = out.tetrahedronlist[4 * ii + nn];
        if (tetgen_moab_map.find(MBVERTEX | (nnn << MB_TYPE_WIDTH)) ==
            tetgen_moab_map.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency between TetGen and MoAB");
        }
        conn[nn] = tetgen_moab_map.at(MBVERTEX | (nnn << MB_TYPE_WIDTH));
      }
      Range tets;
      rval = m_field.get_moab().get_adjacencies(conn, 4, 3, false, tets);
      CHKERRQ_MOAB(rval);
      EntityHandle tet;
      if (tets.empty()) {
        for (int nn = 0; nn != 4; nn++) {
          conn_seq_tet.push_back(conn[nn]);
        }
        tetgen_ii.push_back(ii);
        // rval = m_field.get_moab().create_element(MBTET,conn,4,tet);
        // CHKERRQ_MOAB(rval);
        // Range tet_nodes;
        // rval = m_field.get_moab().get_connectivity(&tet,1,tet_nodes,true);
        // CHKERRQ_MOAB(rval);
        // if(tet_nodes.size()!=4) {
        //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data
        //   inconsistency; tet should have 4 nodes");
        // }
        // Range edges;
        // rval = m_field.get_moab().get_adjacencies(&tet,1,1,true,edges);
        // Range tris;
        // rval = m_field.get_moab().get_adjacencies(&tet,1,2,true,tris);
        // ParallelComm* pcomm =
        // ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
        // BARRIER_RANK_START(pcomm)
        // cerr << pcomm->rank() << " : ";
        // cerr << tet << " : ";
        // // cerr << tet_nodes << endl;
        // cerr << conn[0] << " " << conn[1] << " " << conn[2] << " " << conn[3]
        // << " : ";
        // cerr << edges[0] << " " << edges[1] << " " << edges[2] << " " <<
        // edges[3] << " "
        // << edges[4] << " " << edges[5] << " : ";
        // cerr << tris[0] << " " << tris[1] << " " << tris[2] << " " << tris[3]
        // << endl;
        // BARRIER_RANK_END(pcomm)
      } else {
        // if(tets.size()!=1) {
        //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data
        //   inconsistency; expecting one element");
        // }
        tet = *tets.begin();
        moab_tetgen_map[tet] = iii;
        tetgen_moab_map[iii] = tet;
        if (ents != NULL)
          ents->insert(tet);
      }
    }

    Range new_tets;
    if (!conn_seq_tet.empty()) {
      EntityHandle starte; // Connectivity
      EntityHandle *conn;
      int num_el = tetgen_ii.size();
      rval = iface->get_element_connect(num_el, 4, MBTET, 0, starte, conn);
      CHKERRQ_MOAB(rval);
      std::copy(conn_seq_tet.begin(), conn_seq_tet.end(), conn);
      rval = iface->update_adjacencies(starte, num_el, 4, conn);
      CHKERRQ_MOAB(rval);
      new_tets = Range(starte, starte + num_el - 1);
      std::vector<int>::iterator ii_it = tetgen_ii.begin();
      int ii = 0;
      for (Range::iterator tit = new_tets.begin(); tit != new_tets.end();
           tit++, ii_it++, ii++) {
        // std::sort(&conn_seq_tet[4*ii],&conn_seq_tet[4*ii+4]);
        unsigned long iii = MBTET | ((*ii_it) << MB_TYPE_WIDTH);
        moab_tetgen_map[*tit] = iii;
        tetgen_moab_map[iii] = *tit;
      }
      if (ents != NULL)
        ents->merge(new_tets);
    }

    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::outData(
      tetgenio &in, tetgenio &out,
      std::map<EntityHandle, unsigned long> &moab_tetgen_map,
      std::map<unsigned long, EntityHandle> &tetgen_moab_map, BitRefLevel bit,
      bool id_in_tags, bool error_if_created) {
    MoFEMFunctionBeginHot;

    //
    Interface &m_field = cOre;

    // ParallelComm* pcomm =
    // ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
    //
    // BARRIER_RANK_START(pcomm)
    // {
    //   ostringstream ss;
    //   ss << "t0_" << pcomm->rank() << ".vtk";
    //   rval = m_field.get_moab().write_file(ss.str().c_str(),"VTK");
    //   CHKERRQ_MOAB(rval);
    // }
    // BARRIER_RANK_END(pcomm)

    Range ents;
    ierr = outData(in, out, moab_tetgen_map, tetgen_moab_map, &ents, id_in_tags,
                   error_if_created);
    CHKERRQ(ierr);

    // BARRIER_RANK_START(pcomm)
    // {
    //   ostringstream ss;
    //   ss << "t1_" << pcomm->rank() << ".vtk";
    //   rval = m_field.get_moab().write_file(ss.str().c_str(),"VTK");
    //   CHKERRQ_MOAB(rval);
    // }
    // BARRIER_RANK_END(pcomm)

    // std::cerr << ents.size() << std::endl;
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevel(
        ents.subset_by_type(MBTET), bit);
    CHKERRQ(ierr);

    // BARRIER_RANK_START(pcomm)
    // {
    //   ostringstream ss;
    //   ss << "t2_" << pcomm->rank() << ".vtk";
    //   rval = m_field.get_moab().write_file(ss.str().c_str(),"VTK");
    //   CHKERRQ_MOAB(rval);
    // }
    // BARRIER_RANK_END(pcomm)

    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::setFaceData(
      std::vector<std::pair<Range, int> > &markers, tetgenio &in,
      std::map<EntityHandle, unsigned long> &moab_tetgen_map,
      std::map<unsigned long, EntityHandle> &tetgen_moab_map) {
    MoFEMFunctionBeginHot;
    ErrorCode rval;
    Interface &m_field = cOre;
    in.numberoffacets = markers.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    std::vector<std::pair<Range, int> >::iterator mit = markers.begin();
    for (int ii = 0; mit != markers.end(); mit++, ii++) {
      in.facetmarkerlist[ii] = mit->second;
      Range &faces = mit->first;
      tetgenio::facet *f = &(in.facetlist[ii]);
      f->numberofpolygons = faces.size();
      f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
      int jj = 0;
      for (int dd = 3; dd >= 0; dd--) {
        Range dd_faces = faces.subset_by_dimension(dd);
        if (dd_faces.empty())
          continue;
        Range::iterator it = dd_faces.begin();
        for (; it != dd_faces.end(); it++, jj++) {
          int num_nodes;
          const EntityHandle *conn;
          tetgenio::polygon *p = &(f->polygonlist[jj]);
          switch (m_field.get_moab().type_from_handle(*it)) {
          case MBVERTEX: {
            p->numberofvertices = 1;
            conn = &*it;
          } break;
          default: {
            rval =
                m_field.get_moab().get_connectivity(*it, conn, num_nodes, true);
            CHKERRQ_MOAB(rval);
            p->numberofvertices = num_nodes;
          }
          }
          p->vertexlist = new int[p->numberofvertices];
          for (int nn = 0; nn < p->numberofvertices; nn++) {
            if (moab_tetgen_map.find(conn[nn]) == moab_tetgen_map.end()) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "data inconsistency between TetGen and MoAB");
            }
            p->vertexlist[nn] = moab_tetgen_map[conn[nn]] >> MB_TYPE_WIDTH;
          }
        }
      }
      // holes
      f->numberofholes = 0;
      f->holelist = NULL;
    }
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::getTriangleMarkers(
      std::map<EntityHandle, unsigned long> &tetgen_moab_map, tetgenio &out,
      Range *ents, idxRange_Map *ents_map, bool only_non_zero) {
    MoFEMFunctionBeginHot;
    ErrorCode rval;
    Interface &m_field = cOre;
    Tag th_marker;
    int def_marker = 0;
    rval = m_field.get_moab().tag_get_handle(
        "TETGEN_MARKER", 1, MB_TYPE_INTEGER, th_marker,
        MB_TAG_CREAT | MB_TAG_SPARSE, &def_marker);
    CHKERRQ_MOAB(rval);
    int ii = 0;
    for (; ii < out.numberoftrifaces; ii++) {
      if (only_non_zero) {
        if (out.trifacemarkerlist[ii] == 0) {
          continue;
        }
      }
      EntityHandle conn[3];
      for (int nn = 0; nn < 3; nn++) {
        int iii = MBVERTEX | (out.trifacelist[3 * ii + nn] << MB_TYPE_WIDTH);
        if (tetgen_moab_map.find(iii) == tetgen_moab_map.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency between TetGen and MoAB");
        } else {
          conn[nn] = tetgen_moab_map[iii];
        }
      }
      Range face;
      rval = m_field.get_moab().get_adjacencies(conn, 3, 2, false, face);
      CHKERRQ_MOAB(rval);
      face = face.subset_by_type(MBTRI);
      if (face.size() != 1) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "data inconsistency between TetGen and MoAB, %u", face.size());
      }
      if (ents != NULL)
        ents->merge(face);
      if (ents_map != NULL)
        (*ents_map)[out.trifacemarkerlist[ii]].merge(face);
      rval = m_field.get_moab().tag_set_data(th_marker, &*face.begin(), 1,
                                             &out.trifacemarkerlist[ii]);
      CHKERRQ_MOAB(rval);
    }
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::setRegionData(
      std::vector<std::pair<EntityHandle, int> > &regions, tetgenio &in,
      Tag th) {
    MoFEMFunctionBeginHot;
    ErrorCode rval;
    Interface &m_field = cOre;
    in.numberofregions = regions.size();
    in.regionlist = new double[5 * in.numberofregions];
    int kk = 0;
    std::vector<std::pair<EntityHandle, int> >::iterator it = regions.begin();
    for (int ii = 0; it != regions.end(); it++, ii++) {
      double coords[3];
      switch (m_field.get_moab().type_from_handle(it->first)) {
      case MBVERTEX: {
        if(th) {
          rval = m_field.get_moab().tag_get_data(th, &it->first, 1, coords);
          CHKERRQ_MOAB(rval);
        } else {
          rval = m_field.get_moab().get_coords(&it->first, 1, coords);
          CHKERRQ_MOAB(rval);
        }
      } break;
      case MBTET: {
        int num_nodes;
        const EntityHandle *conn;
        rval = m_field.get_moab().get_connectivity(it->first, conn, num_nodes,
                                                   true);
        CHKERRQ_MOAB(rval);
        double _coords[12];
        if (th) {
          rval = m_field.get_moab().tag_get_data(th, conn, num_nodes, _coords);
          CHKERRQ_MOAB(rval);
        } else {
          rval = m_field.get_moab().get_coords(conn, num_nodes, _coords);
          CHKERRQ_MOAB(rval);
        }
        coords[0] = (_coords[0] + _coords[3] + _coords[6] + _coords[9]) / 4.;
        coords[1] = (_coords[1] + _coords[4] + _coords[7] + _coords[10]) / 4.;
        coords[2] = (_coords[2] + _coords[5] + _coords[8] + _coords[11]) / 4.;
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }
      for (int nn = 0; nn < 3; nn++) {
        in.regionlist[kk++] = coords[nn];
      }
      in.regionlist[kk++] = it->second;
      in.regionlist[kk++] = it->second;
    }
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::getRegionData(
      std::map<EntityHandle, unsigned long> &tetgen_moab_map, tetgenio &out,
      Range *ents, idxRange_Map *ents_map) {
    MoFEMFunctionBeginHot;
    ErrorCode rval;
    Interface &m_field = cOre;
    int nbattributes = out.numberoftetrahedronattributes;
    if (nbattributes == 0) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "tetgen has no regions attributes");
    }
    Tag th_region;
    rval = m_field.get_moab().tag_get_handle("TETGEN_REGION", th_region);
    if (rval == MB_SUCCESS) {
      rval = m_field.get_moab().tag_delete(th_region);
      CHKERRQ_MOAB(rval);
    }
    double def_marker = 0;
    rval = m_field.get_moab().tag_get_handle(
        "TETGEN_REGION", nbattributes, MB_TYPE_DOUBLE, th_region,
        MB_TAG_CREAT | MB_TAG_SPARSE, &def_marker);
    CHKERRQ_MOAB(rval);
    int ii = 0;
    for (; ii < out.numberoftetrahedra; ii++) {
      int jj = 0;
      for (; jj < nbattributes; jj++) {
        double id = out.tetrahedronattributelist[ii * nbattributes + jj];
        int iii = MBTET | (ii << MB_TYPE_WIDTH);
        if (tetgen_moab_map.find(iii) == tetgen_moab_map.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency between TetGen and MoAB");
        }
        EntityHandle ent = tetgen_moab_map[iii];
        rval = m_field.get_moab().tag_set_data(th_region, &ent, 1, &id);
        CHKERRQ_MOAB(rval);
        if (ents != NULL)
          ents->insert(ent);
        if (ents_map != NULL)
          (*ents_map)[id].insert(ent);
      }
    }
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::tetRahedralize(char switches[], tetgenio &in,
                                                 tetgenio &out) {
    MoFEMFunctionBeginHot;
    tetgenbehavior a;
    a.parse_commandline(switches);
    tetrahedralize(&a, &in, &out);
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::loadPoly(char file_name[], tetgenio &in) {
    MoFEMFunctionBeginHot;
    if (!in.load_poly(file_name)) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "can not read TetGen poly file");
    }
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::checkPlanar_Trinagle(double coords[],
                                                       bool *result,
                                                       const double eps) {
    MoFEMFunctionBeginHot;
    double *pa = &coords[0];
    double *pb = &coords[3];
    double *pc = &coords[6];
    double *pd = &coords[9];
    double adx = pa[0] - pd[0];
    double bdx = pb[0] - pd[0];
    double cdx = pc[0] - pd[0];
    double ady = pa[1] - pd[1];
    double bdy = pb[1] - pd[1];
    double cdy = pc[1] - pd[1];
    double adz = pa[2] - pd[2];
    double bdz = pb[2] - pd[2];
    double cdz = pc[2] - pd[2];
    double v = adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) +
               cdx * (ady * bdz - adz * bdy);
    double l = sqrt(pow(pa[0] - pb[0], 2) + pow(pa[1] - pb[1], 2) +
                    pow(pa[2] - pb[2], 2)) +
               sqrt(pow(pa[0] - pc[0], 2) + pow(pa[1] - pc[1], 2) +
                    pow(pa[2] - pc[2], 2)) +
               sqrt(pow(pa[0] - pd[0], 2) + pow(pa[1] - pd[1], 2) +
                    pow(pa[2] - pd[2], 2)) +
               sqrt(pow(pb[0] - pc[0], 2) + pow(pb[1] - pc[1], 2) +
                    pow(pb[2] - pc[2], 2)) +
               sqrt(pow(pb[0] - pd[0], 2) + pow(pb[1] - pd[1], 2) +
                    pow(pb[2] - pd[2], 2)) +
               sqrt(pow(pc[0] - pd[0], 2) + pow(pc[1] - pd[1], 2) +
                    pow(pc[2] - pd[2], 2));
    // std::cerr << fabs(v/pow(l,3)) << " ";
    *result = fabs(v / pow(l, 3)) < eps ? true : false;
    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode
  TetGenInterface::groupPlanar_Triangle(Range &tris, std::vector<Range> &sorted,
                                        const double eps,Tag th) {
    MoFEMFunctionBeginHot;

    Interface &m_field = cOre;

    ErrorCode rval;
    Skinner skin(&m_field.get_moab());

    for (;;) {

      Range noplanar_to_anyother;
      std::vector<Range>::iterator vit = sorted.begin();

      do {

        bool repeat = false;

        // get edges on vit skin
        Range skin_edges;
        rval = skin.find_skin(0, *vit, false, skin_edges);
        CHKERRQ_MOAB(rval);

        // skin edge nodes
        Range skin_edges_nodes;
        rval = m_field.get_moab().get_connectivity(skin_edges, skin_edges_nodes,
                                                   true);
        CHKERRQ_MOAB(rval);

        // get tris adjacent to vit skin edges
        Range skin_edges_tris;
        rval = m_field.get_moab().get_adjacencies(
            skin_edges, 2, false, skin_edges_tris, moab::Interface::UNION);
        CHKERRQ_MOAB(rval);
        // get tris which are part of facet
        Range inner_tris = intersect(skin_edges_tris, *vit);
        Range outer_tris = intersect(skin_edges_tris, tris);

        // tris coplanar with vit tris
        Range coplanar_tris;

        Range::iterator tit = outer_tris.begin();
        for (; tit != outer_tris.end(); tit++) {
          Range tit_conn;
          rval = m_field.get_moab().get_connectivity(&*tit, 1, tit_conn, true);
          CHKERRQ_MOAB(rval);
          tit_conn = subtract(tit_conn, skin_edges_nodes);
          if (tit_conn.empty()) {
            coplanar_tris.insert(*tit);
            noplanar_to_anyother.erase(*tit);
            repeat = true;
          } else {
            Range tit_edges;
            rval = m_field.get_moab().get_adjacencies(
                &*tit, 1, 1, false, tit_edges, moab::Interface::UNION);
            CHKERRQ_MOAB(rval);
            tit_edges = intersect(tit_edges, skin_edges);
            if (tit_edges.size() != 1) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "data inconsistency");
            }
            Range inner_tri;
            rval = m_field.get_moab().get_adjacencies(
                tit_edges, 2, false, inner_tri, moab::Interface::UNION);
            CHKERRQ_MOAB(rval);
            inner_tri = intersect(inner_tri, inner_tris);
            if (inner_tri.size() != 1) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "data inconsistency");
            }
            // get connectivity
            int num_nodes;
            const EntityHandle *inner_tri_conn;
            rval = m_field.get_moab().get_connectivity(
                *inner_tri.begin(), inner_tri_conn, num_nodes, true);
            CHKERRQ_MOAB(rval);
            double coords[12];
            if(th) {
              rval = m_field.get_moab().tag_get_data(th, inner_tri_conn, 3,
                                                     coords);
              CHKERRQ_MOAB(rval);
              rval = m_field.get_moab().tag_get_data(th,&*tit_conn.begin(), 1,
                                                   &coords[9]);
              CHKERRQ_MOAB(rval);
            } else {
              rval = m_field.get_moab().get_coords(inner_tri_conn, 3, coords);
              CHKERRQ_MOAB(rval);
              rval = m_field.get_moab().get_coords(&*tit_conn.begin(), 1,
                                                   &coords[9]);
              CHKERRQ_MOAB(rval);
            }
            bool coplanar;
            ierr = checkPlanar_Trinagle(coords, &coplanar, eps);
            CHKERRQ(ierr);
            if (coplanar) {
              coplanar_tris.insert(*tit);
              noplanar_to_anyother.erase(*tit);
              repeat = true;
            } else {
              noplanar_to_anyother.insert(*tit);
            }
          }
        }

        vit->merge(coplanar_tris);
        tris = subtract(tris, *vit);

        if (repeat) {
          vit = sorted.begin();
        } else {
          vit++;
        }

      } while (vit != sorted.end());

      if (noplanar_to_anyother.empty()) {
        MoFEMFunctionReturnHot(0);
      } else {
        Range seed;
        seed.insert(noplanar_to_anyother[0]);
        tris.erase(noplanar_to_anyother[0]);
        sorted.push_back(seed);
      }
    }

    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode TetGenInterface::groupRegion_Triangle(
      Range &tris, std::vector<std::vector<Range> > &sorted, const double eps) {
    MoFEMFunctionBeginHot;

    // PetscAttachDebugger();

    Range seed;
    seed.insert(tris[0]);
    tris.erase(tris[0]);
    std::vector<Range> vec_seed;
    vec_seed.push_back(seed);
    sorted.push_back(vec_seed);

    for (;;) {
      std::vector<Range> &vec = sorted.back();
      ierr = groupPlanar_Triangle(tris, vec, eps);
      CHKERRQ(ierr);
      if (tris.empty()) {
        MoFEMFunctionReturnHot(0);
      } else {
        Range seed;
        seed.insert(tris[0]);
        tris.erase(tris[0]);
        std::vector<Range> vec_seed;
        vec_seed.push_back(seed);
        sorted.push_back(vec_seed);
      }
    }

    MoFEMFunctionReturnHot(0);
  }
  MoFEMErrorCode TetGenInterface::makePolygonFacet(Range &ents, Range &polygons,
                                                   bool reduce_edges,
                                                   Range *not_reducable_nodes,
                                                   const double eps,Tag th) {
    MoFEMFunctionBeginHot;
    // FIXME: assumes that are no holes

    if (ents.empty()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "no ents to build polygon");
    }

    Interface &m_field = cOre;

    //
    ErrorCode rval;
    Skinner skin(&m_field.get_moab());

    Range skin_edges;
    rval = skin.find_skin(0, ents, false, skin_edges);
    CHKERRQ_MOAB(rval);

    std::vector<EntityHandle> polygon_nodes;
    EntityHandle seed = skin_edges[0];
    Range seen_edges;
    seen_edges.insert(seed);
    skin_edges.erase(seed);
    int num_nodes;
    const EntityHandle *conn;
    rval = m_field.get_moab().get_connectivity(seed, conn, num_nodes, true);
    CHKERRQ_MOAB(rval);
    polygon_nodes.push_back(conn[0]);
    // std::cerr << std::endl;
    // std::cerr << conn[0] << " " << conn[1] << std::endl;
    do {
      EntityHandle last_node = polygon_nodes.back();
      Range adj_edges;
      rval = m_field.get_moab().get_adjacencies(&last_node, 1, 1, false,
                                                adj_edges);
      CHKERRQ_MOAB(rval);
      adj_edges = intersect(adj_edges, skin_edges);
      if (adj_edges.size() == 0) {
        break;
      }
      if (adj_edges.size() != 1) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "should be only one edge");
      }
      seen_edges.insert(adj_edges[0]);
      skin_edges.erase(adj_edges[0]);
      rval = m_field.get_moab().get_connectivity(adj_edges[0], conn, num_nodes,
                                                 true);
      CHKERRQ_MOAB(rval);
      EntityHandle add_node = (last_node == conn[0]) ? conn[1] : conn[0];
      polygon_nodes.push_back(add_node);
      // std::cerr << "\t" << add_node << std::endl;
    } while (1);

    if (reduce_edges) {
      // std::cerr << "polygon " << polygon_nodes.size();
      std::vector<EntityHandle>::iterator pit = polygon_nodes.begin();
      // std::cerr << std::endl;
      for (; pit != polygon_nodes.end();) {
        if (not_reducable_nodes != NULL) {
          if (not_reducable_nodes->find(*pit) != not_reducable_nodes->end()) {
            pit++;
            continue;
          }
        }
        EntityHandle mm;
        if (pit == polygon_nodes.begin()) {
          mm = polygon_nodes.back();
        } else {
          mm = *(pit - 1);
        }
        EntityHandle mc = *pit;
        EntityHandle mp;
        if (polygon_nodes.back() == mc) {
          mp = polygon_nodes[0];
        } else {
          mp = *(pit + 1);
        }
        double coords[9];
        if (th) {
          rval = m_field.get_moab().tag_get_data(th, &mm, 1, &coords[3 * 0]);
          CHKERRQ_MOAB(rval);
          rval = m_field.get_moab().tag_get_data(th, &mc, 1, &coords[3 * 1]);
          CHKERRQ_MOAB(rval);
          rval = m_field.get_moab().tag_get_data(th, &mp, 1, &coords[3 * 2]);
          CHKERRQ_MOAB(rval);
        } else {
          rval = m_field.get_moab().get_coords(&mm, 1, &coords[3 * 0]);
          CHKERRQ_MOAB(rval);
          rval = m_field.get_moab().get_coords(&mc, 1, &coords[3 * 1]);
          CHKERRQ_MOAB(rval);
          rval = m_field.get_moab().get_coords(&mp, 1, &coords[3 * 2]);
          CHKERRQ_MOAB(rval);
        }
        cblas_daxpy(3, -1, &coords[3 * 1], 1, &coords[3 * 0], 1); // mm = mm -
                                                                  // mc
        cblas_daxpy(3, -1, &coords[3 * 1], 1, &coords[3 * 2], 1); // mp = mp -
                                                                  // mc
        double spin[9];
        ierr = Spin(spin, &coords[3 * 0]);
        CHKERRQ(ierr);
        double l0 = cblas_dnrm2(3, &coords[3 * 0], 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1. / l0, spin, 3,
                    &coords[3 * 2], 1, 0., &coords[3 * 1], 1);
        double dot = cblas_dnrm2(3, &coords[3 * 1], 1);
        // std::cerr << mm << " " << mc << " " << mp << " " << dot << std::endl;
        if (dot < eps) {
          polygon_nodes.erase(pit);
          pit = polygon_nodes.begin();
          // std::cerr << std::endl;
        } else {
          pit++;
        }
      }
    }
    // std::cerr << " " << polygon_nodes.size() << std::endl;
    /*pit = polygon_nodes.begin();
    for(;pit!=polygon_nodes.end();pit++) {
      double coords[3];
      rval = m_field.get_moab().get_coords(&*pit,1,coords); CHKERRQ_MOAB(rval);
      std::cerr << *pit << " " << coords[0] << " " << coords[1] << " " <<
    coords[2] << std::endl;
    }*/

    Range existing_polygon;
    rval = m_field.get_moab().get_adjacencies(
        &polygon_nodes[0], polygon_nodes.size(), 2, true, existing_polygon);
    CHKERRQ_MOAB(rval);
    if (existing_polygon.empty()) {
      EntityHandle polygon;
      rval = m_field.get_moab().create_element(MBPOLYGON, &polygon_nodes[0],
                                               polygon_nodes.size(), polygon);
      CHKERRQ_MOAB(rval);
      polygons.insert(polygon);
    } else {
      polygons.merge(existing_polygon);
    }

    MoFEMFunctionReturnHot(0);
  }
  }

#endif //TETGEN
