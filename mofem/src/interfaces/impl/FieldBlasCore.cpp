/** \file FieldBlasCore.cpp
 * \brief Implementation of field algebra
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

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode Core::field_axpy(
  const double alpha,const std::string& field_name_x,const std::string& field_name_y,
  bool error_if_missing,bool creat_if_missing) {
  return FieldBlas(*this).fieldAxpy(alpha,field_name_x,field_name_y,error_if_missing,creat_if_missing);
}
PetscErrorCode Core::set_field(const double val,const EntityType type,const std::string& field_name) {
  return FieldBlas(*this).setField(val,type,field_name);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::set_field(const double val,const EntityType type,const Range &ents,const std::string& field_name) {
  return FieldBlas(*this).setField(val,type,ents,field_name);
}
PetscErrorCode Core::field_scale(const double alpha,const std::string& field_name) {
  return FieldBlas(*this).fieldScale(alpha,field_name);
  MoFEMFunctionReturnHot(0);
}

}
