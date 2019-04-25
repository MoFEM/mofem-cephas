/** \file FieldEvaluator.hpp
 * \brief Field Evaluator
 *
 * Evaluate field at given coordinate 
 *
 *
 * \ingroup field_evaluator
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __FIELD_EVALUATOR_HPP__
#define __FIELD_EVALUATOR_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMFieldEvaluator =
    MOFEMuuid(BitIntefaceId(FIELDEVALUATOR_INTERFACE));

/** \brief Field evaluator interface

  * \ingroup field_evaluator
  */
struct FieldEvaluatorInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  FieldEvaluatorInterface(const MoFEM::Core &core)
      : cOre(const_cast<MoFEM::Core &>(core)) {}

  MoFEMErrorCode buildTree3D(const std::string finite_element);

  MoFEMErrorCode
  evalFEAtThePoint3D(const double *const point, const double distance,
                     const double *const eval_points, const int nb_eval_points,
                     const double eps, const std::string problem,
                     const std::string finite_element,
                     MoFEM::ForcesAndSourcesCore &fe_method, int lower_rank,
                     int upper_rank, MoFEMTypes bh = MF_EXIST,
                     VERBOSITY_LEVELS verb = QUIET);

private:

  EntityHandle rooTreeSet;
  boost::shared_ptr<BVHTree> treePtr;

  boost::function<int(int order_row, int order_col, int order_data)>
      setGaussPtsHook;

  struct SetGaussPts {

    SetGaussPts(MoFEM::ForcesAndSourcesCore &fe_method,
                const double *eval_points, const int nb_eval_points,
                const double eps,
                VERBOSITY_LEVELS verb)
        : feMethod(fe_method), evalPoints(eval_points),
          nbEvalPoints(nb_eval_points), eps(eps), verb(verb) {
      localCoords.resize(nbEvalPoints, 3);
      shapeFunctions.resize(nbEvalPoints, 4);
    }

    MoFEMErrorCode operator()(int order_row, int order_col, int order_data);

  private:

    MoFEM::ForcesAndSourcesCore &feMethod;
    const double *evalPoints;
    const int nbEvalPoints;
    const double eps;
    VERBOSITY_LEVELS verb;

    MatrixDouble localCoords;
    MatrixDouble shapeFunctions;

  };


};

}

#endif // __FIELD_EVALUATOR_HPP__


/**
 * \defgroup field_evaluator Field Evaluator
 * \brief Evaluate field at the point
 *
 * \ingroup mofem
 */