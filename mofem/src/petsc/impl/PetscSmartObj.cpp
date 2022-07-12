/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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