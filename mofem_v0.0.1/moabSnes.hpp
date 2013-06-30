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

#ifndef __MOABSNES_HPP__
#define __MOABSNES_HPP__

#include "moabField.hpp"
#include <petsc.h>
#include <petscmat.h>
#include <petscsnes.h>


namespace MoFEM {

struct moabSnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  moabField &mField;
  Interface &moab;
  ParallelComm* pcomm;

  moabSnesCtx(moabField &_mField): 
    mField(_mField),moab(_mField.get_moab()) {
    pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  }

  const moabField& get_mField() const { return mField; }
  const Interface& get_moab() const { return moab; }

};

}

#endif // __MOABSNES_HPP__
