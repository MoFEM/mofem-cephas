

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