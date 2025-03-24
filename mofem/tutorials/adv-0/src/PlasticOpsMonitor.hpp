

/** \file PlasticOpsMonitor.hpp
 * \example PlasticOpsMonitor.hpp
 */

namespace PlasticOps {

template <int DIM> struct Monitor : public FEMethod {

  Monitor(
      SmartPetscObj<DM> dm,
      std::pair<boost::shared_ptr<PostProcEle>,
                boost::shared_ptr<SkinPostProcEle>>
          pair_post_proc_fe,
      boost::shared_ptr<DomainEle> reaction_fe,
      std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> ux_scatter,
      std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uy_scatter,
      std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uz_scatter,
      std::array<double, DIM> pass_field_eval_coords,
      boost::shared_ptr<SetPtsData> pass_field_eval_data,
      boost::shared_ptr<std::map<std::string, boost::shared_ptr<VectorDouble>>>
          scalar_field_ptrs,
      boost::shared_ptr<std::map<std::string, boost::shared_ptr<MatrixDouble>>>
          vec_field_ptrs,
      boost::shared_ptr<std::map<std::string, boost::shared_ptr<MatrixDouble>>>
          sym_tensor_field_ptrs,
      boost::shared_ptr<std::map<std::string, boost::shared_ptr<MatrixDouble>>>
          tensor_field_ptrs)
      : dM(dm), reactionFe(reaction_fe), uXScatter(ux_scatter),
        uYScatter(uy_scatter), uZScatter(uz_scatter),
        fieldEvalCoords(pass_field_eval_coords),
        fieldEvalData(pass_field_eval_data), scalarFieldPtrs(scalar_field_ptrs),
        vecFieldPtrs(vec_field_ptrs), symTensorFieldPtrs(sym_tensor_field_ptrs),
        tensorFieldPtrs(tensor_field_ptrs) {
    postProcFe = pair_post_proc_fe.first;
    skinPostProcFe = pair_post_proc_fe.second;
  };

  MoFEMErrorCode preProcess() { return 0; }
  MoFEMErrorCode operator()() { return 0; }

private:
  SmartPetscObj<DM> dM;
  boost::shared_ptr<PostProcEle> postProcFe;
  boost::shared_ptr<SkinPostProcEle> skinPostProcFe;
  boost::shared_ptr<DomainEle> reactionFe;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uXScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uYScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uZScatter;
  std::array<double, DIM> fieldEvalCoords;
  boost::shared_ptr<SetPtsData> fieldEvalData;
  boost::shared_ptr<std::map<std::string, boost::shared_ptr<VectorDouble>>>
      scalarFieldPtrs;
  boost::shared_ptr<std::map<std::string, boost::shared_ptr<MatrixDouble>>>
      vecFieldPtrs;
  boost::shared_ptr<std::map<std::string, boost::shared_ptr<MatrixDouble>>>
      symTensorFieldPtrs;
  boost::shared_ptr<std::map<std::string, boost::shared_ptr<MatrixDouble>>>
      tensorFieldPtrs;

protected:
  MoFEMErrorCode postProcess();
};

template <int DIM> MoFEMErrorCode Monitor<DIM>::postProcess() {
  MoFEMFunctionBegin;

  MoFEM::Interface *m_field_ptr;
  CHKERR DMoFEMGetInterfacePtr(dM, &m_field_ptr);
  auto *simple = m_field_ptr->getInterface<Simple>();

  if (do_eval_field) {

    CHKERR m_field_ptr->getInterface<FieldEvaluatorInterface>()
        ->evalFEAtThePoint<DIM>(
            fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
            simple->getDomainFEName(), fieldEvalData,
            m_field_ptr->get_comm_rank(), m_field_ptr->get_comm_rank(),
            getCacheWeakPtr().lock(), MF_EXIST, QUIET);

    auto process_scalar_field =
        [](const std::string label,
           const boost::shared_ptr<VectorDouble> scalarFieldPtr) {
          if (scalarFieldPtr->size()) {
            auto t_scalar_holder = getFTensor0FromVec(*scalarFieldPtr);

            MOFEM_LOG("SYNC", Sev::inform)
                << "For " << label << " field, " << t_scalar_holder;
          }
        };

    auto process_vector_field =
        [](const std::string label,
           const boost::shared_ptr<MatrixDouble> vecFieldPtr) {
          if (vecFieldPtr->size1()) {
            auto t_vec_holder = getFTensor1FromMat<DIM>(*vecFieldPtr);

            std::string vec_holder_string = "";
            for (int i = 0; i < DIM; i++) {
              vec_holder_string +=
                  " " + boost::lexical_cast<std::string>(t_vec_holder(i));
            }

            MOFEM_LOG("SYNC", Sev::inform)
                << "For " << label << " field," << vec_holder_string;
          }
        };

    auto process_sym_tensor_field =
        [](const std::string label,
           const boost::shared_ptr<MatrixDouble> symTensorFieldPtr) {
          if (symTensorFieldPtr->size1()) {
            auto t_sym_tensor_holder =
                getFTensor2SymmetricFromMat<DIM>(*symTensorFieldPtr);

            std::string sym_tensor_holder_string = "";
            for (int i = 0; i < DIM; i++) {
              for (int j = i; j < DIM; j++) {
                sym_tensor_holder_string +=
                    ", entry " + boost::lexical_cast<std::string>(i) +
                    boost::lexical_cast<std::string>(j) + " = ";
                sym_tensor_holder_string +=
                    boost::lexical_cast<std::string>(t_sym_tensor_holder(i, j));
              }
            }

            MOFEM_LOG("SYNC", Sev::inform)
                << "For " << label << " field" << sym_tensor_holder_string;
          }
        };

    auto process_tensor_field =
        [](const std::string label,
           const boost::shared_ptr<MatrixDouble> tensorFieldPtr) {
          if (tensorFieldPtr->size1()) {
            auto t_tensor_holder =
                getFTensor2FromMat<DIM, DIM>(*tensorFieldPtr);

            std::string tensor_holder_string = "";
            for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                tensor_holder_string +=
                    ", entry " + boost::lexical_cast<std::string>(i) +
                    boost::lexical_cast<std::string>(j) + " = ";
                tensor_holder_string +=
                    boost::lexical_cast<std::string>(t_tensor_holder(i, j));
              }
            }

            MOFEM_LOG("SYNC", Sev::inform)
                << "For " << label << " field" << tensor_holder_string;
          }
        };

    auto process_fields = [](auto &fieldPtr, auto process_field) {
      for (const auto &pair : fieldPtr) {
        const std::string &label = pair.first;
        const boost::shared_ptr<MatrixDouble> &ptr = pair.second;

        if (ptr) {
          process_field(label, ptr);
        }
      }
    };

    if (scalarFieldPtrs) {
      for (const auto &pair : *scalarFieldPtrs) {
        const std::string &label = pair.first;
        const boost::shared_ptr<VectorDouble> &ptr = pair.second;

        if (ptr) {
          process_scalar_field(label, ptr);
        }
      }
    }

    if (vecFieldPtrs) {
      for (const auto &pair : *vecFieldPtrs) {
        const std::string &label = pair.first;
        const boost::shared_ptr<MatrixDouble> &ptr = pair.second;

        if (ptr) {
          process_vector_field(label, ptr);
        }
      }
    }

    if (symTensorFieldPtrs) {
      for (const auto &pair : *symTensorFieldPtrs) {
        const std::string &label = pair.first;
        const boost::shared_ptr<MatrixDouble> &ptr = pair.second;

        if (ptr) {
          process_sym_tensor_field(label, ptr);
        }
      }
    }

    if (tensorFieldPtrs) {
      for (const auto &pair : *tensorFieldPtrs) {
        const std::string &label = pair.first;
        const boost::shared_ptr<MatrixDouble> &ptr = pair.second;

        if (ptr) {
          process_tensor_field(label, ptr);
        }
      }
    }
  }

  MOFEM_LOG_SEVERITY_SYNC(m_field_ptr->get_comm(), Sev::inform);

  auto make_vtk = [&]() {
    MoFEMFunctionBegin;
    if (postProcFe) {
      CHKERR DMoFEMLoopFiniteElements(dM, "dFE", postProcFe, getCacheWeakPtr());
      CHKERR postProcFe->writeFile(
          "out_plastic_" + boost::lexical_cast<std::string>(ts_step) + ".h5m");
    }
    if (skinPostProcFe) {
      CHKERR DMoFEMLoopFiniteElements(dM, "bFE", skinPostProcFe,
                                      getCacheWeakPtr());
      CHKERR skinPostProcFe->writeFile(
          "out_skin_plastic_" + boost::lexical_cast<std::string>(ts_step) +
          ".h5m");
    }
    MoFEMFunctionReturn(0);
  };

  auto calculate_reaction = [&]() {
    MoFEMFunctionBegin;
    auto r = createDMVector(dM);
    reactionFe->f = r;
    CHKERR VecZeroEntries(r);
    CHKERR DMoFEMLoopFiniteElements(dM, "dFE", reactionFe, getCacheWeakPtr());

#ifndef NDEBUG
    auto post_proc_residual = [&](auto dm, auto f_res, auto out_name) {
      MoFEMFunctionBegin;
      auto post_proc_fe =
          boost::make_shared<PostProcBrokenMeshInMoab<DomainEle>>(*m_field_ptr);
      using OpPPMap = OpPostProcMapInMoab<DIM, DIM>;
      auto u_vec = boost::make_shared<MatrixDouble>();
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<DIM>("U", u_vec, f_res));
      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

              {},

              {{"RES", u_vec}},

              {}, {})

      );

      CHKERR DMoFEMLoopFiniteElements(dM, "dFE", post_proc_fe);
      post_proc_fe->writeFile("res.h5m");
      MoFEMFunctionReturn(0);
    };

    CHKERR post_proc_residual(dM, r, "reaction");
#endif // NDEBUG

    MoFEMFunctionReturn(0);
  };

  auto print_max_min = [&](auto &tuple, const std::string msg) {
    MoFEMFunctionBegin;
    CHKERR VecScatterBegin(std::get<1>(tuple), ts_u, std::get<0>(tuple),
                           INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecScatterEnd(std::get<1>(tuple), ts_u, std::get<0>(tuple),
                         INSERT_VALUES, SCATTER_FORWARD);
    double max, min;
    CHKERR VecMax(std::get<0>(tuple), PETSC_NULL, &max);
    CHKERR VecMin(std::get<0>(tuple), PETSC_NULL, &min);
    MOFEM_LOG_C("PLASTICITY", Sev::inform, "%s time %3.4e min %3.4e max %3.4e",
                msg.c_str(), ts_t, min, max);
    MoFEMFunctionReturn(0);
  };

  int se = 1;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-save_every", &se, PETSC_NULL);
  if (!(ts_step % se)) {
    CHKERR make_vtk();
  }
  if (reactionFe)
    CHKERR calculate_reaction();
  CHKERR print_max_min(uXScatter, "Ux");
  CHKERR print_max_min(uYScatter, "Uy");
  if (DIM == 3)
    CHKERR print_max_min(uZScatter, "Uz");

  MoFEMFunctionReturn(0);
}

}; // namespace PlasticOps