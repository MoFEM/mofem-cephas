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


#ifndef __DMMMOFEM_H
#define __DMMMOFEM_H

#define DMMOFEM "mofem"

PetscErrorCode DMMoFEMCreate(MPI_Comm comm, DM *mofem);

PetscErrorCode DMMoFEMSetParallelComm(DM dm,moab::ParallelComm *pcomm);
PetscErrorCode DMMoFEMGetParallelComm(DM dm,moab::ParallelComm **pcomm);
PetscErrorCode DMMoabSetInterface(DM dm,moab::Interface *iface);
PetscErrorCode DMMoabGetInterface(DM dm,moab::Interface **iface);



#endif //__DMMMOFEM_H

