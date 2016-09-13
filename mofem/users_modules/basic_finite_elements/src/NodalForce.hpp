/** \file NodalForce.hpp
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

#ifndef __NODAL_FORCES_HPP__
#define __NODAL_FORCES_HPP__

/** \brief Force applied to nodes
  * \ingroup mofem_static_boundary_conditions
  */
struct NodalForce {

  MoFEM::Interface &mField;
  NodalForce(MoFEM::Interface &m_field): mField(m_field),fe(m_field) {}

  struct MyFE: public MoFEM::VertexElementForcesAndSourcesCore {
    MyFE(MoFEM::Interface &m_field);
  };

  MyFE fe;
  MyFE& getLoopFe() { return fe; }

  struct bCForce {
    ForceCubitBcData data;
    Range nOdes;
  };
  std::map<int,bCForce> mapForce;

  boost::ptr_vector<MethodForForceScaling> methodsOp;

  /// \brief Operator to assemble nodal force into right hand side vector
  struct OpNodalForce: public MoFEM::VertexElementForcesAndSourcesCore::UserDataOperator {

    Vec F;
    bool useSnesF;
    bCForce &dAta;
    boost::ptr_vector<MethodForForceScaling> &methodsOp;

    OpNodalForce(const std::string field_name,Vec _F,bCForce &data,
      boost::ptr_vector<MethodForForceScaling> &methods_op,
      bool use_snes_f = false);

    ublas::vector<FieldData> Nf;

    /** Exeuted for each entity on element, i.e. in this case Vertex element
      has only one entity, that is vertex
    */
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

  };

  PetscErrorCode addForce(const std::string field_name,Vec F,int ms_id,bool use_snes_f = false);

};

struct MetaNodalForces {

  /** \brief Scale force based on tag value "_LoadFactor_Scale_"

    This is obsolete, is kept to have back compatibility with fracture code.

  */
  struct TagForceScale: public MethodForForceScaling {
    MoFEM::Interface &mField;
    double *sCale;
    Tag thScale;

    TagForceScale(MoFEM::Interface &m_field);
    PetscErrorCode scaleNf(const FEMethod *fe,ublas::vector<FieldData> &Nf);

  };

  /// Add element taking information from NODESET
  static PetscErrorCode addElement (MoFEM::Interface &m_field,const std::string field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    ierr = m_field.add_finite_element("FORCE_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("FORCE_FE",field_name); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      Range tris;
      rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
      Range edges;
      rval = m_field.get_moab().get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERRQ_MOAB(rval);
      Range tris_nodes;
      rval = m_field.get_moab().get_connectivity(tris,tris_nodes); CHKERRQ_MOAB(rval);
      Range edges_nodes;
      rval = m_field.get_moab().get_connectivity(edges,edges_nodes); CHKERRQ_MOAB(rval);
      Range nodes;
      rval = m_field.get_moab().get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERRQ_MOAB(rval);
      nodes = subtract(nodes,tris_nodes);
      nodes = subtract(nodes,edges_nodes);
      ierr = m_field.add_ents_to_finite_element_by_VERTICEs(nodes,"FORCE_FE"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  /// Set integration point operators
  static PetscErrorCode setOperators(
    MoFEM::Interface &m_field, boost::ptr_map<std::string,NodalForce> &nodal_forces, Vec F,const std::string field_name
  ) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    string fe_name = "FORCE_FE";
    nodal_forces.insert(fe_name,new NodalForce(m_field));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = nodal_forces.at(fe_name).addForce(field_name,F,it->getMeshsetId());  CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

};

#endif //__NODAL_FORCES_HPP__
