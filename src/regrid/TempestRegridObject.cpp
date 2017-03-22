///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestRegridObject.cpp
///	\author  Paul Ullrich
///	\version March 10, 2017
///
///	<remarks>
///		Copyright 2016- Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "TempestRegridObject.h"
#include "GridObject.h"
#include "FiniteElementTools.h"
#include "DataArray3D.h"

#include "OverlapMesh.h"
#include "LinearRemapSE0.h"

namespace HRegrid {

///////////////////////////////////////////////////////////////////////////////
// TempestRegridObjectConstructor
///////////////////////////////////////////////////////////////////////////////

std::string TempestRegridObjectConstructor::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	TempestRegridObject * pobjTempestRegrid = new TempestRegridObject("");
	if (pobjTempestRegrid == NULL) {
		_EXCEPTIONT("Unable to initialize TempestRegridObject");
	}

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjTempestRegrid;
	}

	// Initialize the grid with given parameters
	std::string strError =
		pobjTempestRegrid->Initialize(
			objreg,
			vecCommandLine,
			vecCommandLineType);

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjTempestRegrid;
	} else {
		delete pobjTempestRegrid;
	}

	return strError;
}

///////////////////////////////////////////////////////////////////////////////
// TempestRegridObject
///////////////////////////////////////////////////////////////////////////////

std::string TempestRegridObject::Initialize(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecFuncArguments,
	const std::vector<ObjectType> & vecFuncArgumentsType
) {

	// Set the return value
	if (vecFuncArguments.size() != 3) {
		return std::string("ERROR: Invalid tempestregrid() arguments:\n")
			+ std::string("Expected tempestremap([source grid], "
				"[target grid], [parameters])");
	}

	// Source grid
	GridObject * pobjSourceGrid =
		dynamic_cast<GridObject *>(
			objreg.GetObject(vecFuncArguments[0]));
	if (pobjSourceGrid == NULL) {
		return std::string("ERROR: First argument of tempestregrid() must be"
			" of type grid");
	}

	// Target grid
	GridObject * pobjTargetGrid =
		dynamic_cast<GridObject *>(
			objreg.GetObject(vecFuncArguments[1]));
	if (pobjTargetGrid == NULL) {
		return std::string("ERROR: Second argument of tempestregrid() must be"
			" of type grid");
	}

	// Parameters
	Object * pobjParameters =
		dynamic_cast<Object *>(
			objreg.GetObject(vecFuncArguments[2]));
	if (pobjTargetGrid == NULL) {
		return std::string("ERROR: Third argument to tempestregrid() must be"
			" of type parameter_list");
	}

	// Monotone type
	int nMonotoneType = 0;

	// No conservation
	bool fNoConservation = false;

	// Always bubble
	bool fBubble = true;

	// Concave faces in mesh
	bool fSourceConcave = false;

	// Concave faces in mesh
	bool fTargetConcave = false;

	// Generate EdgeMap on source mesh if not available
	Mesh & meshSource = pobjSourceGrid->GetMesh();
	meshSource.ConstructEdgeMap();
	if (meshSource.vecFaceArea.GetRows() != meshSource.dLat.GetRows()) {
		meshSource.CalculateFaceAreas(fSourceConcave);
	}

	// Generate EdgeMap on target mesh if not available
	Mesh & meshTarget = pobjTargetGrid->GetMesh();
	meshTarget.ConstructEdgeMap();
	if (meshTarget.vecFaceArea.GetRows() != meshTarget.dLat.GetRows()) {
		meshTarget.CalculateFaceAreas(fTargetConcave);
	}

	// Generate the overlap mesh
	Mesh * pMeshOverlap = new Mesh();

	AnnounceStartBlock("Generate overlap mesh");
	GenerateOverlapMesh_v2(
		pobjSourceGrid->GetMesh(),
		pobjTargetGrid->GetMesh(),
		*pMeshOverlap,
		HRegrid::OverlapMeshMethod_Fuzzy);
	AnnounceEndBlock("Done");

	// Generate offline map
	AnnounceStartBlock("Generate offline map");

	if ((meshSource.eDataLayout == Mesh::DataLayout_SpectralElementGLL) ||
		(meshSource.eDataLayout == Mesh::DataLayout_DiscontinuousGLL)
	) {
		bool fContinuousIn =
			(meshSource.eDataLayout == Mesh::DataLayout_SpectralElementGLL);

		if (meshTarget.eDataLayout == Mesh::DataLayout_Volumetric) {
			DataArray3D<int> dataGLLNodes;
			DataArray3D<double> dataGLLJacobian;

			double dNumericalArea =
				GenerateMetaData(
					meshSource,
					meshSource.nP,
					fBubble,
					dataGLLNodes,
					dataGLLJacobian);

			LinearRemapSE4(
				meshSource,
				meshTarget,
				*pMeshOverlap,
				dataGLLNodes,
				dataGLLJacobian,
				nMonotoneType,
				fContinuousIn,
				fNoConservation,
				m_mapRemap);

		} else {
			_EXCEPTIONT("Not implemented");
		}

	} else {
		_EXCEPTIONT("Not implemented");
	}

	AnnounceEndBlock("Done");

	delete pMeshOverlap;

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

};

