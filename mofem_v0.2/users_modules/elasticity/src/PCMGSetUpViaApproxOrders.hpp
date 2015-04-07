/** \file PCMGSetUpViaApproxOrders.hpp
 * \brief header of multi-grid solver for p- adaptivity 
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

#ifndef __PCMGSETUP_VIA_APPROX_ORDERS_HPP__
#define __PCMGSETUP_VIA_APPROX_ORDERS_HPP__

PetscErrorCode PCMGSetUpViaApproxOrders(PC pc,FieldInterface *mfield_ptr,const char problem_name[],int verb = 0);

#endif //__PCMGSETUP_VIA_APPROX_ORDERS_HPP__

