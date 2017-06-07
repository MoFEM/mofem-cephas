/** \file ProblemCore.cpp
 * \brief Managing complexities for problem
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

  PetscErrorCode Simple::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMSimple) {
      *iface = dynamic_cast<Simple*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  Simple::Simple(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
  meshSet(0) {
    PetscLogEventRegister("Simple",0,&USER_EVENT_Simple);
  }
  Simple::~Simple() {
  }

  PetscErrorCode Simple::addDomainField(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_of_cooficients,
    const TagType tag_type,
    const enum MoFEMTypes bh,
    int verb
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    ierr = m_field.add_field(
      name, space, base, nb_of_cooficients, tag_type, bh, verb
    ); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

}
