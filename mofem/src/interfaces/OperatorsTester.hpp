/**
 * @file OperatorsTester.hpp
 * @brief
 * @version 0.13.2
 * @date 2022-12-04
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef __OPERATORS_TESTER_HPP__
#define __OPERATORS_TESTER_HPP__

namespace MoFEM {

struct OperatorsTester : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  OperatorsTester(const MoFEM::Core &core);
  virtual ~OperatorsTester() = default;

  using RandomFieldData = std::pair<std::string, std::array<double, 2>>;

  SmartPetscObj<Vec> setRandomFields(SmartPetscObj<DM> dm,
                                     std::vector<RandomFieldData> random_fields);

  SmartPetscObj<Vec> assembleVec(SmartPetscObj<DM> dm, std::string fe_name,
                                 boost::shared_ptr<FEMethod> pipeline,
                                 SmartPetscObj<Vec> x,
                                 SmartPetscObj<Vec> delta_x,
                                 SmartPetscObj<Vec> delta2_x, double time,
                                 double delta_t, CacheTupleWeakPtr cache_ptr);

  SmartPetscObj<Mat> assembleMat(SmartPetscObj<DM> dm, std::string fe_name,
                                 boost::shared_ptr<FEMethod> pipeline,
                                 SmartPetscObj<Vec> x,
                                 SmartPetscObj<Vec> delta_x,
                                 SmartPetscObj<Vec> delta2_x, double time,
                                 double delta_t, CacheTupleWeakPtr cache_ptr);

  SmartPetscObj<Vec> directionalCentralFiniteDiffence(
      SmartPetscObj<DM> dm, std::string fe_name,
      boost::shared_ptr<FEMethod> pipeline, SmartPetscObj<Vec> x,
      SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
      SmartPetscObj<Vec> diff_x, double time, double delta_t, double eps,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

  SmartPetscObj<Vec> checkCentralFiniteDifference(
      SmartPetscObj<DM> dm, std::string fe_name,
      boost::shared_ptr<FEMethod> pipeline_rhs,
      boost::shared_ptr<FEMethod> pipeline_lhs, SmartPetscObj<Vec> x,
      SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
      SmartPetscObj<Vec> diff_x, double time, double delta_t, double eps,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

private:
  MoFEM::Core &cOre;
  std::pair<SmartPetscObj<Vec>, SmartPetscObj<Vec>>
  setPipelineX(boost::shared_ptr<FEMethod> pipeline, SmartPetscObj<Vec> x,
               SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
               double delta_t);
};

} // namespace MoFEM

#endif //__OPERATORS_TESTER_HPP__