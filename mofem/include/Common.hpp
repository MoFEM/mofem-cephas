/** \file Common.hpp
 * \brief Basic structures and data
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

#ifndef __COMMON_HPP__
#define __COMMON_HPP__

namespace MoFEM {

/**
 * @brief No entity handle is indicated by zero handle, i.e. root meshset
 *
 */
const EntityHandle no_handle = 0;

} // namespace MoFEM

#include <Exceptions.hpp>
#include <Types.hpp>
#include <Templates.hpp>

#endif //__COMMON_HPP__

/**
 * \defgroup mofem MoFEM
 */
