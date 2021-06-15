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

template <> void intrusive_ptr_release<Vec>(Vec obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = VecDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<Mat>(Mat obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = MatDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<DM>(DM obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = DMDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<IS>(IS obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = ISDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<AO>(AO obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = AODestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<KSP>(KSP obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = KSPDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<SNES>(SNES obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = SNESDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}

template <> void intrusive_ptr_release<TS>(TS obj) {
  int cnt = 0;
  PetscErrorCode ierr =
      PetscObjectGetReference(MoFEM::getPetscObject(obj), &cnt);
  if (!ierr) {
    if (cnt) {
      ierr = TSDestroy(&obj);
      auto comm = PetscObjectComm(MoFEM::getPetscObject(obj));
      CHKERRABORT(comm, ierr);
    }
  }
}