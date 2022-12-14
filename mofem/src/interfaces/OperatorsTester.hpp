/**
 * @file OperatorsTester.hpp
 * @brief Used to calculate directives and testing tangent matrix, i.e. hessian.
 * @version 0.13.2
 * @date 2022-12-04
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef __OPERATORS_TESTER_HPP__
#define __OPERATORS_TESTER_HPP__

namespace MoFEM {

/**
 * @brief Calculate directional derivative of the right hand side and compare it
 * with tangent matrix derivative.
 *
 */
struct OperatorsTester : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  OperatorsTester(const MoFEM::Core &core);
  virtual ~OperatorsTester() = default;

  using RandomFieldData = std::pair<std::string, std::array<double, 2>>;

  /**
   * @brief Generate random fileds
   *
   * Example: generate random vector for DM (problem) from simple interface,
   where FIELD random values of DOFs are in range from -1 to 1, and FIELD2
   random values are in range from 0 to 1.
   * @code {.cpp}
   *  auto x = opt->setRandomFields(simple->getDM(),
                                  {{"FIELD1", {-1, 1}}, {"FIELD2", {0,1}}});
   * @endcode
   *
   *
   * TODO: Set random field to specific entities, and potentially order for
   testing proposes to dissect error in tangent matrix.
   *
   * @param dm
   * @param random_fields look at definition @ref RandomFieldData
   * @return SmartPetscObj<Vec> smart vector
   */
  SmartPetscObj<Vec> setRandomFields(SmartPetscObj<DM> dm,
                                     std::vector<RandomFieldData> random_fields);

  /**
   * @brief Assemble the right hand side vector
   *
   * @param dm
   * @param fe_name // fe name
   * @param pipeline // pipeline, i.e. fe instance
   * @param x // problem (dm) vector
   * @param delta_x // vector for x rate, can be null, i.e, SmartPetscObj<Vec>()
   * @param delta2_x // vector for x second rate, i.e. acceleration
   * @param time // time
   * @param delta_t // time increment
   * @param cache_ptr // finite element data cache, can be null
   * @return SmartPetscObj<Vec>
   */
  SmartPetscObj<Vec> assembleVec(SmartPetscObj<DM> dm, std::string fe_name,
                                 boost::shared_ptr<FEMethod> pipeline,
                                 SmartPetscObj<Vec> x,
                                 SmartPetscObj<Vec> delta_x,
                                 SmartPetscObj<Vec> delta2_x, double time,
                                 double delta_t, CacheTupleWeakPtr cache_ptr);

  /**
   * @brief Assemble the left hand side vector
   *
   * @param dm
   * @param fe_name // fe name
   * @param pipeline // pipeline, i.e. fe instance
   * @param x // problem (dm) vector
   * @param delta_x // vector for x rate, can be null, i.e, SmartPetscObj<Vec>()
   * @param delta2_x // vector for x second rate, i.e. acceleration
   * @param time // time
   * @param delta_t // time increment
   * @param cache_ptr // finite element data cache, can be null
   * @return SmartPetscObj<Vec>
   */
  SmartPetscObj<Mat> assembleMat(SmartPetscObj<DM> dm, std::string fe_name,
                                 boost::shared_ptr<FEMethod> pipeline,
                                 SmartPetscObj<Vec> x,
                                 SmartPetscObj<Vec> delta_x,
                                 SmartPetscObj<Vec> delta2_x, double time,
                                 double delta_t, CacheTupleWeakPtr cache_ptr);

  /**
   * @brief Calculate directional directive using finite difference
   * 
   * @param dm 
   * @param fe_name 
   * @param pipeline 
   * @param x 
   * @param delta_x 
   * @param delta2_x 
   * @param diff_x // direction of derivative
   * @param time 
   * @param delta_t 
   * @param eps 
   * @param cache_ptr 
   * @return SmartPetscObj<Vec> 
   */
  SmartPetscObj<Vec> directionalCentralFiniteDifference(
      SmartPetscObj<DM> dm, std::string fe_name,
      boost::shared_ptr<FEMethod> pipeline, SmartPetscObj<Vec> x,
      SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
      SmartPetscObj<Vec> diff_x, double time, double delta_t, double eps,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

  /**
   * @brief Check consistency between directional derivative with matrix
   * 
   * @param dm 
   * @param fe_name 
   * @param pipeline_rhs 
   * @param pipeline_lhs 
   * @param x 
   * @param delta_x 
   * @param delta2_x 
   * @param diff_x 
   * @param time 
   * @param delta_t 
   * @param eps 
   * @param cache_ptr 
   * @return SmartPetscObj<Vec> 
   */
  SmartPetscObj<Vec> checkCentralFiniteDifference(
      SmartPetscObj<DM> dm, std::string fe_name,
      boost::shared_ptr<FEMethod> pipeline_rhs,
      boost::shared_ptr<FEMethod> pipeline_lhs, SmartPetscObj<Vec> x,
      SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
      SmartPetscObj<Vec> diff_x, double time, double delta_t, double eps,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr());

private:
  MoFEM::Core &cOre;

  /**
   * @brief Set vectors, x, x_t, and x_tt to finite element instance
   *
   * Finite element instance is a pipeline. x_t and x_tt are evaluated for given
   * delta_x, delta2_x and delta_t.
   *
   * @param pipeline
   * @param x
   * @param delta_x
   * @param delta2_x
   * @param delta_t
   * @return std::pair<SmartPetscObj<Vec>, SmartPetscObj<Vec>>
   */
  std::pair<SmartPetscObj<Vec>, SmartPetscObj<Vec>>
  setPipelineX(boost::shared_ptr<FEMethod> pipeline, SmartPetscObj<Vec> x,
               SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
               double delta_t);
};

} // namespace MoFEM

#endif //__OPERATORS_TESTER_HPP__