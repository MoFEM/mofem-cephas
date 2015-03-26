/** \file PCMGSetUpViaApproxOrders.cpp
 * \brief useful compiler directives and definitions
 * 
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

#include <MoFEM.hpp>
using namespace MoFEM;
#include <PCMGSetUpViaApproxOrders.hpp>

struct PCMGSetUpViaApproxOrdersCtx {

  FieldInterface *mField_ptr;		///< MoFEM interface
  string problemName;			///< Problem name

  PCMGSetUpViaApproxOrdersCtx(FieldInterface *mfield_ptr,string &problem_name): 
    mField_ptr(m_field),problemName(problem_name) {
  }

  PetscErrorCode getOptions() {
    PetscFunctionBegin;
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","MOFEM Multi-Grid (Orders) pre-conditioner","none"); CHKERRQ(ierr);

    ierr = PetscOptionsInt("-my_order",
      "default approximation order","",
      2,&order,PETSC_NULL); CHKERRQ(ierr);


    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

PetcErrorCode PCMGSetUpViaApproxOrders(PC pc,FieldInterface *m_fild_ptr,const char problem_name[]) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscFunctionBegin;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)pc,&comm); CHKERRQ(ierr);
  int result = 0;
  MPI_Comm_compare(comm,m_field_ptr->get_comm(),&result);
  if(result > MPI_CONGRUENT) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"MoFEM and PC have to use the same communicator");
  }
  
  rval = PCMGSetUpViaApproxOrdersCtx ctx(m_field_ptr,problem_name); CHKERR(rval);


  
  PetscFunctionReturn(0);
}
