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

#ifndef NDEBUG
#define BOOST_DISABLE_ASSERTS ///< If not set makes BOOST freezing, check what
                              ///< happen fi boost recompiled with cxx17
                              ///< standard
#define BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE
#endif

#ifndef __MOFEM_HPP__
#define __MOFEM_HPP__

// Include system and libraries files
#include <Includes.hpp>

// SRC APPROXIMATION
#include <config.h>
#include <definitions.h>

// FTensor
#include <FTensor.hpp>

namespace MoFEM {
// FIXME: All operators in FTensor move to FTensor namespace.
using FTensor::operator<<;
using FTensor::operator>>;
} // namespace MoFEM

#include <Common.hpp>
#include <UnknownInterface.hpp>
#include <DeprecatedPetsc.hpp>

// SRC/APPROXIMATION
#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <Hdiv.hpp>
#include <Hcurl.hpp>
#include <BernsteinBezier.hpp>
#include <fem_tools.h>
#include <BaseFunction.hpp>
#include <LegendrePolynomial.hpp>
#include <LobattoPolynomial.hpp>
#include <JacobiPolynomial.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <QuadPolynomialBase.hpp> // Base functions on quad
#include <EdgePolynomialBase.hpp>
#include <TriPolynomialBase.hpp>
#include <TetPolynomialBase.hpp>
#include <HexPolynomialBase.hpp>
#include <FatPrismPolynomialBase.hpp>
#include <FlatPrismPolynomialBase.hpp>
#include <EdgeQuadHexPolynomials.hpp>

// SRC/MULTI-INDICES
#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <RefEntsMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <FieldEntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <RefElementMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <SeriesMultiIndices.hpp>

// SRC/INTERFACES
#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <Core.hpp>

#include <AuxPETSc.hpp>

#include <LogManager.hpp>
#include <BitRefManager.hpp>
#include <Tools.hpp>
#include <CommInterface.hpp>
#include <ISManager.hpp>
#include <VecManager.hpp>
#include <FieldBlas.hpp>
#include <ProblemsManager.hpp>
#include <MatrixManager.hpp>
#include <Simple.hpp>
#include <MeshRefinement.hpp>
#include <SeriesRecorder.hpp>
#include <PrismInterface.hpp>
#include <MeshsetsManager.hpp>
#ifdef WITH_TETGEN
#include <TetGenInterface.hpp>
#endif // WITH_TETGEN
#ifdef WITH_MED
#include <MedInterface.hpp>
#endif // WITH_MED
#include <CutMeshInterface.hpp>
#include <NodeMerger.hpp>
#include <PrismsFromSurfaceInterface.hpp>

// SRC/PETSC
#include <KspCtx.hpp>
#include <SnesCtx.hpp>
#include <TsCtx.hpp>
#include <DMMoFEM.hpp>

// SRC/FINITE_ELEMENTS
#include <EntitiesFieldData.hpp>
#include <DataOperators.hpp>
#include <ForcesAndSourcesCore.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>
#include <EdgeElementForcesAndSourcesCore.hpp>
#include <VertexElementForcesAndSourcesCore.hpp>
#include <FlatPrismElementForcesAndSourcesCore.hpp>
#include <ContactPrismElementForcesAndSourcesCore.hpp>
#include <FatPrismElementForcesAndSourcesCore.hpp>
#include <VolumeElementForcesAndSourcesCoreOnSide.hpp>
#include <VolumeElementForcesAndSourcesCoreOnContactPrismSide.hpp>
#include <FaceElementForcesAndSourcesCoreOnSide.hpp>
#include <FaceElementForcesAndSourcesCoreOnParent.hpp>
#include <EdgeElementForcesAndSourcesCoreOnParent.hpp>
#include <Projection10NodeCoordsOnField.hpp>
#include <UserDataOperators.hpp>
#include <HODataOperators.hpp> // Manage HO order geometry
#include <MeshProjectionDataOperators.hpp> // Operators for projections between bit ref levels
#include <BaseDerivativesDataOperators.hpp> // Operators to calculate HO direcarives
#include <FormsIntegrators.hpp>
#include <LinearFormsIntegrators.hpp>
#include <BiLinearFormsIntegrators.hpp>

// More interfaces

#include <PipelineManager.hpp>
#include <FieldEvaluator.hpp>
#include <BcManager.hpp>

#endif // MOFEM_HPP__
