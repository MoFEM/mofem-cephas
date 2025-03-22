/** \file mesh_data_transfer.cpp

  \brief Transfer data from source mesh to target mesh using MOAB

*/

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    char mesh_source_file[255] = "source.h5m";
    char mesh_target_file[255] = "target.h5m";
    char mesh_out_file[255] = "out.h5m";
    char iterp_tag_name[255] = "INTERNAL_STRESS";

    int interp_order = 1;
    int interp_tag_len = 9;

    PetscBool set_source_data = PETSC_FALSE;
    PetscBool use_target_verts = PETSC_TRUE;
    PetscBool atom_test = PETSC_FALSE;

    double toler = 5.e-10;

    CHKERR PetscOptionsBegin(m_field.get_comm(), "", "Read MED tool", "none");
    CHKERR PetscOptionsString("-source_file", "source mesh file name", "",
                              "source.h5m", mesh_source_file, 255, PETSC_NULL);
    CHKERR PetscOptionsString("-target_file", "target mesh file name", "",
                              "target.h5m", mesh_target_file, 255, PETSC_NULL);
    CHKERR PetscOptionsString("-output_file", "output mesh file name", "",
                              "out.h5m", mesh_out_file, 255, PETSC_NULL);
    CHKERR PetscOptionsString("-interp_tag", "Interpolation tag name", "",
                              "INTERNAL_STRESS", iterp_tag_name, 255,
                              PETSC_NULL);
    CHKERR PetscOptionsBool("-set_source_data", "Set source data", "",
                            set_source_data, &set_source_data, PETSC_NULL);
    CHKERR PetscOptionsBool("-use_target_verts",
                            "use target vertices for interpolation", "",
                            use_target_verts, &use_target_verts, PETSC_NULL);
    CHKERR PetscOptionsBool("-atom_test", "is atom test", "", atom_test,
                            &atom_test, PETSC_NULL);
    CHKERR PetscOptionsInt("-interp_order", "interpolation order", "", 0,
                           &interp_order, PETSC_NULL);
    CHKERR PetscOptionsInt("-interp_tag_length", "interpolation tag length", "",
                           9, &interp_tag_len, PETSC_NULL);

    ierr = PetscOptionsEnd();
    CHKERRQ(ierr);

    Coupler::Method method = Coupler::CONSTANT;

    if (interp_order == 1) {
      method = Coupler::LINEAR_FE;
    }

    std::vector<std::string> meshFiles(2);
    meshFiles[0] = string(mesh_source_file);
    meshFiles[1] = string(mesh_target_file);

    std::string readOpts, writeOpts;

    int nprocs, rank;
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    CHKERRQ(ierr);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRQ(ierr);

    readOpts =
        (nprocs > 1
             ? "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARTITION_"
               "DISTRIBUTE;"
               "PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1"
             : "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARTITION_"
               "DISTRIBUTE;"
               "PARALLEL_RESOLVE_SHARED_ENTS");

    writeOpts = (nprocs > 1 ? "PARALLEL=WRITE_PART" : "");

    std::vector<ParallelComm *> pcs(meshFiles.size());

    EntityHandle *roots =
        (EntityHandle *)malloc(sizeof(EntityHandle) * meshFiles.size());

    moab::Interface *mbImpl = &moab;

    for (unsigned int i = 0; i < meshFiles.size(); i++) {
      pcs[i] = new ParallelComm(mbImpl, MPI_COMM_WORLD);
      int index = pcs[i]->get_id();
      std::string newReadopts;
      std::ostringstream extraOpt;
      extraOpt << ";PARALLEL_COMM=" << index;
      newReadopts = readOpts + extraOpt.str();

      CHKERR mbImpl->create_meshset(MESHSET_SET, roots[i]);

      CHKERR mbImpl->load_file(meshFiles[i].c_str(), &roots[i],
                               newReadopts.c_str());
    }

    Range src_elems, targ_elems, targ_verts;
    CHKERR pcs[0]->get_part_entities(src_elems, 3);

    Tag interp_tag;
    double def_val[9];
    bzero(def_val, 9 * sizeof(double));
    CHKERR mbImpl->tag_get_handle(iterp_tag_name, interp_tag_len,
                                  MB_TYPE_DOUBLE, interp_tag,
                                  MB_TAG_CREAT | MB_TAG_DENSE, &def_val);

    // Set source tag to create sslv116 test mesh
    if (set_source_data) {
      std::vector<double> spos;
      spos.resize(3 * src_elems.size());
      CHKERR mbImpl->get_coords(src_elems, &spos[0]);

      std::vector<double> stress(9 * src_elems.size(), 0.0);
      auto t_stress = getFTensor2FromPtr<3, 3>(&stress[0]);

      auto temperature = [](double z) -> double {
        return (z < -1)   ? (120.0 / -9.0) * (z + 10) + 250
               : (z <= 2) ? (50.0 / -3.0) * (z + 1) + 130
                          : (170.0 / 8.0) * (z - 2) + 80;
      };

      double young_modulus = 2e11;
      double poisson_ratio = 0.3;
      double bulk_modulus_K = young_modulus / (3 * (1 - 2 * poisson_ratio));
      double shear_modulus_G = young_modulus / (2 * (1 + poisson_ratio));
      double alpha = 1e-5;
      double temp_0 = 250;

      FTensor::Index<'i', 3> i;
      FTensor::Index<'j', 3> j;
      FTensor::Index<'k', 3> k;
      FTensor::Index<'l', 3> l;

      MatrixDouble mat_D_ptr;
      mat_D_ptr.resize(6, 6, false);
      auto t_D = getFTensor4DdgFromPtr<3, 3, 0>(&*(mat_D_ptr.data().begin()));
      constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
      t_D(i, j, k, l) = 2 * shear_modulus_G * ((t_kd(i, k) ^ t_kd(j, l)) / 4.) +
                        (bulk_modulus_K - (2. / 3.) * shear_modulus_G) *
                            t_kd(i, j) * t_kd(k, l);

      for (auto ielem = 0; ielem < src_elems.size(); ++ielem) {
        double z = spos[3 * ielem + 2];
        double temp = temperature(z);
        t_stress(i, j) =
            -t_D(i, j, k, l) * t_kd(k, l) * (temp - temp_0) * alpha;
        ++t_stress;
      }

      CHKERR mbImpl->tag_set_data(interp_tag, src_elems, &stress[0]);
    }

    Coupler mbc(mbImpl, pcs[0], src_elems, 0, true);

    std::vector<double> vpos; // the positions we are interested in
    int numPointsOfInterest = 0;

    Range tmp_verts;

    // First get all vertices adj to partition entities in target mesh
    CHKERR pcs[1]->get_part_entities(targ_elems, 3);
    if (!use_target_verts) {
      targ_verts = targ_elems;
    } else {
      CHKERR mbImpl->get_adjacencies(targ_elems, 0, false, targ_verts,
                                     moab::Interface::UNION);
    }

    // Then get non-owned verts and subtract
    CHKERR pcs[1]->get_pstatus_entities(0, PSTATUS_NOT_OWNED, tmp_verts);
    targ_verts = subtract(targ_verts, tmp_verts);

    // get position of these entities; these are the target points
    numPointsOfInterest = (int)targ_verts.size();
    vpos.resize(3 * targ_verts.size());
    CHKERR mbImpl->get_coords(targ_verts, &vpos[0]);

    // Locate those points in the source mesh
    std::cout << "rank " << pcs[0]->proc_config().proc_rank()
              << " points of interest: " << numPointsOfInterest << "\n";
    CHKERR mbc.locate_points(&vpos[0], numPointsOfInterest, 0, toler);

    std::vector<double> source_data(interp_tag_len * src_elems.size(), 0.0);
    std::vector<double> target_data(interp_tag_len * numPointsOfInterest, 0.0);

    CHKERR mbImpl->tag_get_data(interp_tag, src_elems, &source_data[0]);

    Tag scalar_tag, adj_count_tag;
    double def_scl = 0;
    string scalar_tag_name = string(iterp_tag_name) + "_COMP";
    CHKERR mbImpl->tag_get_handle(scalar_tag_name.c_str(), 1, MB_TYPE_DOUBLE,
                                  scalar_tag, MB_TAG_CREAT | MB_TAG_DENSE,
                                  &def_scl);

    string adj_count_tag_name = "ADJ_COUNT";
    double def_adj = 0;
    CHKERR mbImpl->tag_get_handle(adj_count_tag_name.c_str(), 1, MB_TYPE_DOUBLE,
                                  adj_count_tag, MB_TAG_CREAT | MB_TAG_DENSE,
                                  &def_adj);

    std::vector<double> source_data_scalar(src_elems.size());

    for (int itag = 0; itag < interp_tag_len; itag++) {

      for (int ielem = 0; ielem < src_elems.size(); ielem++) {
        source_data_scalar[ielem] = source_data[itag + ielem * interp_tag_len];
      }

      CHKERR mbImpl->tag_set_data(scalar_tag, src_elems,
                                  &source_data_scalar[0]);

      if (interp_order) {

        Range src_verts;
        CHKERR mbImpl->get_connectivity(src_elems, src_verts, true);

        CHKERR mbImpl->tag_clear_data(scalar_tag, src_verts, &def_val);
        CHKERR mbImpl->tag_clear_data(adj_count_tag, src_verts, &def_adj);

        for (auto &tet : src_elems) {
          double tet_data = 0;
          CHKERR mbImpl->tag_get_data(scalar_tag, &tet, 1, &tet_data);

          Range adj_verts;
          CHKERR mbImpl->get_connectivity(&tet, 1, adj_verts, true);

          std::vector<double> adj_vert_data(adj_verts.size(), 0.0);
          std::vector<double> adj_vert_count(adj_verts.size(), 0.0);

          CHKERR mbImpl->tag_get_data(scalar_tag, adj_verts, &adj_vert_data[0]);
          CHKERR mbImpl->tag_get_data(adj_count_tag, adj_verts,
                                      &adj_vert_count[0]);

          for (int ivert = 0; ivert < adj_verts.size(); ivert++) {
            adj_vert_data[ivert] += tet_data;
            adj_vert_count[ivert] += 1;
          }

          CHKERR mbImpl->tag_set_data(scalar_tag, adj_verts, &adj_vert_data[0]);
          CHKERR mbImpl->tag_set_data(adj_count_tag, adj_verts,
                                      &adj_vert_count[0]);
        }

        std::vector<Tag> tags;
        tags.push_back(scalar_tag);
        tags.push_back(adj_count_tag);
        pcs[0]->reduce_tags(tags, tags, MPI_SUM, src_verts);

        std::vector<double> src_vert_data(src_verts.size(), 0.0);
        std::vector<double> src_vert_adj_count(src_verts.size(), 0.0);

        CHKERR mbImpl->tag_get_data(scalar_tag, src_verts, &src_vert_data[0]);
        CHKERR mbImpl->tag_get_data(adj_count_tag, src_verts,
                                    &src_vert_adj_count[0]);

        for (int ivert = 0; ivert < src_verts.size(); ivert++) {
          src_vert_data[ivert] /= src_vert_adj_count[ivert];
        }
        CHKERR mbImpl->tag_set_data(scalar_tag, src_verts, &src_vert_data[0]);
      }

      std::vector<double> target_data_scalar(numPointsOfInterest, 0.0);
      CHKERR mbc.interpolate(method, scalar_tag_name, &target_data_scalar[0]);

      for (int ielem = 0; ielem < numPointsOfInterest; ielem++) {
        target_data[itag + ielem * interp_tag_len] = target_data_scalar[ielem];
      }
    }

    CHKERR mbImpl->tag_delete(scalar_tag);
    CHKERR mbImpl->tag_delete(adj_count_tag);

    // Use original tag
    CHKERR mbImpl->tag_set_data(interp_tag, targ_verts, &target_data[0]);

    if (use_target_verts) {

      std::vector<Tag> tags;
      tags.push_back(interp_tag);
      pcs[1]->reduce_tags(tags, tags, MPI_SUM, targ_verts);

      for (auto &tet : targ_elems) {
        Range adj_verts;
        CHKERR mbImpl->get_connectivity(&tet, 1, adj_verts, true);

        std::vector<double> adj_vert_data(adj_verts.size() * interp_tag_len,
                                          0.0);
        CHKERR mbImpl->tag_get_data(interp_tag, adj_verts, &adj_vert_data[0]);

        std::vector<double> tet_data(interp_tag_len, 0.0);
        for (int itag = 0; itag < interp_tag_len; itag++) {
          for (int i = 0; i < adj_verts.size(); i++) {
            tet_data[itag] += adj_vert_data[i * interp_tag_len + itag];
          }
          tet_data[itag] /= adj_verts.size();
        }
        CHKERR mbImpl->tag_set_data(interp_tag, &tet, 1, &tet_data[0]);
      }
    }

    if (atom_test) {

      auto compute_tet_volume = [](const std::vector<double> &vpos) -> double {
        // Extract vertices
        double x0 = vpos[0], y0 = vpos[1], z0 = vpos[2];
        double x1 = vpos[3], y1 = vpos[4], z1 = vpos[5];
        double x2 = vpos[6], y2 = vpos[7], z2 = vpos[8];
        double x3 = vpos[9], y3 = vpos[10], z3 = vpos[11];

        // Compute vectors
        double ax = x1 - x0, ay = y1 - y0, az = z1 - z0;
        double bx = x2 - x0, by = y2 - y0, bz = z2 - z0;
        double cx = x3 - x0, cy = y3 - y0, cz = z3 - z0;

        // Compute cross product b x c
        double cross_x = by * cz - bz * cy;
        double cross_y = bz * cx - bx * cz;
        double cross_z = bx * cy - by * cx;

        // Compute dot product a . (b x c)
        double dot = ax * cross_x + ay * cross_y + az * cross_z;

        // Compute volume
        return std::abs(dot) / 6.0;
      };

      std::vector<double> data_integ(interp_tag_len, 0.0);
      for (auto &tet : targ_elems) {

        std::vector<double> tet_data(interp_tag_len, 0.0);
        CHKERR mbImpl->tag_get_data(interp_tag, &tet, 1, &tet_data[0]);

        Range tmp_verts;
        CHKERR mbImpl->get_connectivity(&tet, 1, tmp_verts);
        std::vector<double> vpos(3 * tmp_verts.size());
        CHKERR mbImpl->get_coords(tmp_verts, &vpos[0]);
        double vol = compute_tet_volume(vpos);

        for (int itag = 0; itag < interp_tag_len; itag++) {
          data_integ[itag] += tet_data[itag] * vol;
        }
      }

      std::vector<double> global_data_integ(interp_tag_len, 0.0);
      MPI_Allreduce(&data_integ[0], &global_data_integ[0], interp_tag_len,
                    MPI_DOUBLE, MPI_SUM, m_field.get_comm());

      if (0 == rank) {
        std::cout << "Integrated stress for sslv116 test: ";
        for (auto &val : global_data_integ) {
          std::cout << val << " ";
        }
        std::cout << "\n";
      }

      double non_zero_val = 1.655e12;
      double non_zero_tol = 4e-3;
      if (interp_order || use_target_verts) {
        non_zero_tol = 2e-3;
      }
      int non_zero_inds[] = {0, 4, 8};
      bool non_zero_check =
          all_of(begin(non_zero_inds), end(non_zero_inds), [&](int i) {
            return abs(global_data_integ[i] - non_zero_val) / non_zero_val <
                   non_zero_tol;
          });
      int zero_inds[] = {1, 2, 3, 5, 6, 7};
      bool zero_check = all_of(begin(zero_inds), end(zero_inds), [&](int i) {
        return abs(global_data_integ[i]) < 1e-12;
      });
      if (!non_zero_check || !zero_check) {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                "Wrong value of the integrated stress");
      }
    }

    Range partSets;
    // Only save the target mesh
    partSets.insert((EntityHandle)roots[1]);
    std::string newwriteOpts;
    std::ostringstream extraOpt;
    if (nprocs > 1) {
      int index = pcs[1]->get_id();
      extraOpt << ";PARALLEL_COMM=" << index;
    }
    newwriteOpts = writeOpts + extraOpt.str();
    CHKERR mbImpl->write_file(mesh_out_file, NULL, newwriteOpts.c_str(),
                              partSets);
    if (0 == rank) {
      std::cout << "Wrote file " << mesh_out_file << std::endl;
    }

    free(roots);

    for (unsigned int i = 0; i < meshFiles.size(); i++)
      delete pcs[i];
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
