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

#ifndef __MOFEM_HPP__
#define __MOFEM_HPP__

//Include system and libraries files
#include <Includes.hpp>

//SRC APPROXIMATION
#include <version.h>
#include <config.h>
#include <definitions.h>

//FTensor
#include <FTensor.hpp>
#include <Common.hpp>
#include <UnknownInterface.hpp>

//SRC/APPROXIMATION
#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>
#include <BaseFunction.hpp>
#include <LegendrePolynomial.hpp>
#include <LobattoPolynomial.hpp>

//SRC/MULTI-INDICES
#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <SeriesMultiIndices.hpp>

//SRC/INTERFACES
#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

//SRC/PETSC
#include <KspCtx.hpp>
#include <SnesCtx.hpp>
#include <TsCtx.hpp>
#include <DMMoFEM.hpp>

//SRC/FINITE_ELEMENTS
#include <DataStructures.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <TetPolynomialBase.hpp> // Base functions on tet
#include <TriPolynomialBase.hpp> // Base functions on tri
#include <EdgePolynomialBase.hpp> // Base functions on tri
#include <FlatPrismPolynomialBase.hpp> // Base functions on tri
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <UserDataOperators.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>
#include <EdgeElementForcesAndSurcesCore.hpp>
#include <VertexElementForcesAndSourcesCore.hpp>
#include <FlatPrismElementForcesAndSurcesCore.hpp>
#include <FatPrismElementForcesAndSurcesCore.hpp>
#include <Projection10NodeCoordsOnField.hpp>

// interfaces
#ifdef WITH_TETGEN
  #include <TetGenInterface.hpp>
#endif //WITH_TETGEN
#include <BitLevelCoupler.hpp>
#include <NodeMerger.hpp>
#include <PrismsFromSurfaceInterface.hpp>

#endif //MOFEM_HPP__
