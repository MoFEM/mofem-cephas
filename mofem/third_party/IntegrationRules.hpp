/**
 * \file IntegrationRules.hpp
 *
 * Integration rules
 *
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

namespace MoFEM {

namespace IntRules {

namespace NCC {

#include <triangle_ncc_rule.h>
#include <tetrahedron_ncc_rule.h>

} // namespace NCC

namespace NCO {

#include <triangle_nco_rule.h>
#include <tetrahedron_nco_rule.h>

} // namespace NCO

} // namespace IntRules

}