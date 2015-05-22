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

namespace MoFEM {

  /** \brief Force on edges and lines
  */
  struct EdgeForce {

    FieldInterface &mField;
    EdgeForce(FieldInterface &m_field): mField(m_field),fe(m_field) {}

    struct MyFE: public EdgeElementForcesAndSurcesCore {
      MyFE(FieldInterface &m_field): EdgeElementForcesAndSurcesCore(m_field) {}
    };

    MyFE fe;
    MyFE& getLoopFe() { return fe; }

    struct bCForce {
      ForceCubitBcData data;
      Range eDges;
    };
    map<int,bCForce> mapForce;

    boost::ptr_vector<MethodsForOp> methodsOp;

    struct OpEdgeForce: public EdgeElementForcesAndSurcesCore::UserDataOperator {

      Vec &F;
      bCForce &dAta;
      boost::ptr_vector<MethodsForOp> &methodsOp;
      bool useSnesF;

      OpEdgeForce(
        const string field_name,Vec &f,bCForce &data,
        boost::ptr_vector<MethodsForOp> &methods_op,
        bool use_snes_f = false
      );


      VectorDouble wEights;
      VectorDouble Nf;

      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);


    };

    PetscErrorCode addForce(const string field_name,Vec &F,int ms_id,bool use_snes_f = false);

  };

  struct MetaEdgeForces {

    /// Add element taking information from NODESET
    static PetscErrorCode addEdgeForceElement (FieldInterface &mField,const string field_name) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = mField.add_finite_element("FORCE_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("FORCE_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("FORCE_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("FORCE_FE",field_name); CHKERRQ(ierr);
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
        Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
        Range edges;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
        Range tris_edges;
        rval = mField.get_moab().get_adjacencies(tris,1,false,tris_edges,Interface::UNION); CHKERR_PETSC(rval);
        edges = subtract(edges,tris_edges);
        ierr = mField.add_ents_to_finite_element_by_EDGEs(edges,"FORCE_FE"); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

    /// Set integration point operators
    static PetscErrorCode setEdgeForceElementOperators(
      FieldInterface &mField,
      boost::ptr_map<string,EdgeForce> &edge_forces,
      Vec &F,const string field_name
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      string fe_name = "FORCE_FE";
      edge_forces.insert(fe_name,new EdgeForce(mField));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
        ierr = edge_forces.at(fe_name).addForce(field_name,F,it->get_msId());  CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

}

#endif //__EDGE_FORCE_HPP__
