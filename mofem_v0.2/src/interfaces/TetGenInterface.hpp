/** \file TetGebInterface.hpp
 * \brief TetGen interface 
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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

#ifndef __TETGENINTERFACE_HPP__
#define __TETGENINTERFACE_HPP__

#include "FieldUnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMFieldInterface = MOFEMuuid( BitIntefaceId(TETGEN_INTERFACE) );

struct TetGenInterface: public FieldUnknownInterface {

  FieldInterface& mField;
  TetGenInterface(FieldInterface& m_field): mField(m_field) {};

  PetscErrorCode tEtraedralize(
    char *switches,
    EntityHandle in_meshset,const BitRefLevel &bit_in,
    EntityHandle out_meshset,const BitRefLevel &bit_out);


};

}
