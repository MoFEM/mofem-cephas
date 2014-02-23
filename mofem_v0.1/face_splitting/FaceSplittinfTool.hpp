/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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


#ifndef __FACESPLITTINGTOOL_HPP__
#define __FACESPLITTINGTOOL_HPP__

#include "FieldInterface.hpp"
#include "CoreDataStructures.hpp"

struct FaceSplittingTools {

  FieldInterface& mField;

  //This is moab mesh to do work
  Interface& moab_work;
  Core mb_instance_work;

  FaceSplittingTools (FieldInterface& _mField): 
    mField(_mField),moab_work(mb_instance_work) {
  }

  


}

#endif // __FACESPLITTINGTOOL_HPP__

