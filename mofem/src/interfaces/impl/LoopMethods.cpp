/** \file LoopMethods.cpp
\brief methods for managing loops
*/

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#include <MoFEM.hpp>

namespace MoFEM {

constexpr PetscData::Switches PetscData::CtxSetNone;
constexpr PetscData::Switches PetscData::CtxSetF;
constexpr PetscData::Switches PetscData::CtxSetA;
constexpr PetscData::Switches PetscData::CtxSetB;
constexpr PetscData::Switches PetscData::CtxSetX;
constexpr PetscData::Switches PetscData::CtxSetX_T;
constexpr PetscData::Switches PetscData::CtxSetX_TT;
constexpr PetscData::Switches PetscData::CtxSetTime;

// PetscData
MoFEMErrorCode
PetscData::query_interface(boost::typeindex::type_index type_index,
                           UnknownInterface **iface) const {
  *iface = const_cast<PetscData *>(this);
  return 0;
}

PetscData::PetscData()
    : f(PETSC_NULL), A(PETSC_NULL), B(PETSC_NULL), x(PETSC_NULL),
      x_t(PETSC_NULL), x_tt(PETSC_NULL) {}

// KSP
MoFEMErrorCode
KspMethod::query_interface(boost::typeindex::type_index type_index,
                           UnknownInterface **iface) const {
  *iface = const_cast<KspMethod *>(this);
  return 0;
};

KspMethod::KspMethod()
    : ksp_ctx(CTX_KSPNONE), ksp(PETSC_NULL), ksp_f(PetscData::f),
      ksp_A(PetscData::A), ksp_B(PetscData::B) {}

MoFEMErrorCode KspMethod::copyKsp(const KspMethod &ksp) {
  MoFEMFunctionBeginHot;
  this->ksp_ctx = ksp.ksp_ctx;
  this->ksp = ksp.ksp;
  this->ksp_f = ksp.ksp_f;
  this->ksp_A = ksp.ksp_A;
  this->ksp_B = ksp.ksp_B;
  MoFEMFunctionReturnHot(0);
}

// SNES
MoFEMErrorCode SnesMethod::query_interface(boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<SnesMethod *>(this);
  return 0;
}

SnesMethod::SnesMethod()
    : snes_ctx(CTX_SNESNONE), snes_x(PetscData::x), snes_f(PetscData::f),
      snes_A(PetscData::A), snes_B(PetscData::B) {}

MoFEMErrorCode SnesMethod::copySnes(const SnesMethod &snes) {
  MoFEMFunctionBeginHot;
  this->snes_ctx = snes.snes_ctx;
  this->snes = snes.snes;
  this->snes_x = snes.snes_x;
  this->snes_f = snes.snes_f;
  this->snes_A = snes.snes_A;
  this->snes_B = snes.snes_B;
  MoFEMFunctionReturnHot(0);
}

// TS
MoFEMErrorCode
TSMethod::query_interface(boost::typeindex::type_index type_index,
                          UnknownInterface **iface) const {
  *iface = const_cast<TSMethod *>(this);
  return 0;
}

TSMethod::TSMethod()
    : ts_ctx(CTX_TSNONE), ts_step(-1), ts_a(0), ts_t(0), ts_u(PetscData::x),
      ts_u_t(PetscData::x_t), ts_u_tt(PetscData::x_tt), ts_F(PetscData::f),
      ts_A(PetscData::A), ts_B(PetscData::B) {}

MoFEMErrorCode TSMethod::copyTs(const TSMethod &ts) {
  MoFEMFunctionBeginHot;
  this->ts_ctx = ts.ts_ctx;
  this->ts = ts.ts;
  this->ts_u = ts.ts_u;
  this->ts_u_t = ts.ts_u_t;
  this->ts_u_tt = ts.ts_u_tt;
  this->ts_F = ts.ts_F;
  this->ts_A = ts.ts_A;
  this->ts_B = ts.ts_B;
  this->ts_step = ts.ts_step;
  this->ts_a = ts.ts_a;
  this->ts_t = ts.ts_t;
  MoFEMFunctionReturnHot(0);
}

// BasicMethod
BasicMethod::BasicMethod()
    : PetscData(), KspMethod(), SnesMethod(), TSMethod(), nInTheLoop(0),
      loopSize(0), rAnk(-1), sIze(-1), refinedEntitiesPtr(nullptr),
      refinedFiniteElementsPtr(nullptr), problemPtr(nullptr),
      fieldsPtr(nullptr), entitiesPtr(nullptr), dofsPtr(nullptr),
      finiteElementsPtr(nullptr), finiteElementsEntitiesPtr(nullptr),
      adjacenciesPtr(nullptr) {}

MoFEMErrorCode BasicMethod::copyBasicMethod(const BasicMethod &basic) {
  MoFEMFunctionBeginHot;

  this->nInTheLoop = basic.nInTheLoop;
  this->loopSize = basic.loopSize;
  this->rAnk = basic.rAnk;
  this->sIze = basic.sIze;
  this->refinedEntitiesPtr = basic.refinedEntitiesPtr;
  this->refinedFiniteElementsPtr = basic.refinedFiniteElementsPtr;
  this->problemPtr = basic.problemPtr;
  this->fieldsPtr = basic.fieldsPtr;
  this->entitiesPtr = basic.entitiesPtr;
  this->dofsPtr = basic.dofsPtr;
  this->finiteElementsPtr = basic.finiteElementsPtr;
  this->finiteElementsEntitiesPtr = basic.finiteElementsEntitiesPtr;
  this->adjacenciesPtr = basic.adjacenciesPtr;

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BasicMethod::preProcess() {
  MoFEMFunctionBeginHot;
  if (preProcessHook) {
    ierr = preProcessHook();
    CHKERRG(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "should be implemented by user in derived class (preProcess)");
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode BasicMethod::postProcess() {
  MoFEMFunctionBeginHot;
  if (postProcessHook) {
    ierr = postProcessHook();
    CHKERRG(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "should be implemented by user in derived class (postProcess)");
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode BasicMethod::operator()() {
  MoFEMFunctionBeginHot;
  if (operatorHook) {
    ierr = operatorHook();
    CHKERRG(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "should be implemented by user in derived class (operator)");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FEMethod::getNodeData(const std::string field_name,
                                     VectorDouble &data,
                                     const bool reset_dofs) {
  MoFEMFunctionBegin;

  // TODO: [CORE-60] Fix implementation not to use DOFs

  auto get_nodes_field_data =
      [&](boost::shared_ptr<FEDofEntity_multiIndex> &&dofs,
          VectorDouble &nodes_data) {
        MoFEMFunctionBegin;

        auto bit_number = getFieldBitNumber(field_name);
        auto &dofs_by_uid = dofs->get<Unique_mi_tag>();
        auto dit =
            dofs_by_uid.lower_bound(FieldEntity::getLocalUniqueIdCalculate(
                bit_number, get_id_for_min_type<MBVERTEX>()));
        auto hi_dit =
            dofs_by_uid.upper_bound(FieldEntity::getLocalUniqueIdCalculate(
                bit_number, get_id_for_max_type<MBVERTEX>()));

        if (dit != hi_dit) {
          auto &first_dof = **dit;
          const int num_nodes = getNumberOfNodes();
          const int nb_dof_idx = first_dof.getNbOfCoeffs();
          const int max_nb_dofs = nb_dof_idx * num_nodes;
          nodes_data.resize(max_nb_dofs, false);
          nodes_data.clear();

          for (; dit != hi_dit;) {
            const auto &dof_ptr = *dit;
            const auto &dof = *dof_ptr;
            const auto &sn = *dof.getSideNumberPtr();
            const int side_number = sn.side_number;

            int pos = side_number * nb_dof_idx;
            auto ent_filed_data_vec = dof.getEntFieldData();
            for (int ii = 0; ii != nb_dof_idx; ++ii) {
              nodes_data[pos] = ent_filed_data_vec[ii];
              ++pos;
              ++dit;
            }
          }

        } else if (reset_dofs) {
          nodes_data.resize(0, false);
        }

        MoFEMFunctionReturn(0);
      };

  return get_nodes_field_data(getDataDofsPtr(), data);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
