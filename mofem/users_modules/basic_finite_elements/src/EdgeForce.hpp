
/** \file EdgeForce.hpp
  \ingroup mofem_static_boundary_conditions
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

#ifndef __EDGE_FORCE_HPP__
#define __EDGE_FORCE_HPP__

/** \brief Force on edges and lines
*/
struct EdgeForce {

  MoFEM::Interface &mField;
  EdgeForce(MoFEM::Interface &m_field): mField(m_field),fe(m_field,1){}

  struct MyFE: public MoFEM::EdgeElementForcesAndSurcesCore {
    int addToRule;
    MyFE(MoFEM::Interface &m_field,int add_to_rule):
    EdgeElementForcesAndSurcesCore(m_field),
    addToRule(add_to_rule)
    {}
      int getRule(int order) { return order+addToRule; };
    };

    MyFE fe;
    MyFE& getLoopFe() { return fe; }

    struct bCForce {
      ForceCubitBcData data;
      Range eDges;
    };
    std::map<int,bCForce> mapForce;

    boost::ptr_vector<MethodForForceScaling> methodsOp;

    struct OpEdgeForce: public MoFEM::EdgeElementForcesAndSurcesCore::UserDataOperator {

      Vec F;
      bCForce &dAta;
      boost::ptr_vector<MethodForForceScaling> &methodsOp;
      bool useSnesF;

      OpEdgeForce(
        const std::string field_name,Vec f,bCForce &data,
        boost::ptr_vector<MethodForForceScaling> &methods_op,
        bool use_snes_f = false
      );


      VectorDouble wEights;
      VectorDouble Nf;

      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);


    };

    PetscErrorCode addForce(const std::string field_name,Vec F,int ms_id,bool use_snes_f = false);

  };

  struct MetaEdgeForces {

    /// Add element taking information from NODESET
    static PetscErrorCode addElement (MoFEM::Interface &m_field,const std::string field_name) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = m_field.add_finite_element("FORCE_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("FORCE_FE",field_name); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("FORCE_FE",field_name); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("FORCE_FE",field_name); CHKERRQ(ierr);
      if(m_field.check_field("MESH_NODE_POSITIONS")) {
        ierr = m_field.modify_finite_element_add_field_data("FORCE_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      }
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
        Range tris;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
        Range edges;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERRQ_MOAB(rval);
        Range tris_edges;
        rval = m_field.get_moab().get_adjacencies(tris,1,false,tris_edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
        edges = subtract(edges,tris_edges);
        ierr = m_field.add_ents_to_finite_element_by_EDGEs(edges,"FORCE_FE"); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

    /// Set integration point operators
    static PetscErrorCode setOperators(
      MoFEM::Interface &m_field,
      boost::ptr_map<std::string,EdgeForce> &edge_forces,
      Vec F,const std::string field_name
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      string fe_name = "FORCE_FE";
      edge_forces.insert(fe_name,new EdgeForce(m_field));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
        ierr = edge_forces.at(fe_name).addForce(field_name,F,it->getMeshsetId());  CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };


#endif //__EDGE_FORCE_HPP__
