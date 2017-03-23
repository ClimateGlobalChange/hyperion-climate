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

std::string TempestRegridObject::Call(
	const ObjectRegistry & objreg,
	const std::string & strFunctionName,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	if (strFunctionName == "regrid") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"regrid\"");
		}
		return std::string("Not implemented");
	}

	// Output the map as a netcdf file
	if (strFunctionName == "outputnetcdf") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"outputnetcdf\"");
		}

		m_mapRemap.Write(vecCommandLine[0]);

		return std::string("");

	} else if (strFunctionName == "regrid") {
	}

	return
		Object::Call(
			objreg,
			strFunctionName,
			vecCommandLine,
			vecCommandLineType,
			ppReturn);
}

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
	meshSource.ConstructReverseNodeArray();
	meshSource.ConstructEdgeMap();
	double dSourceArea = meshSource.CalculateFaceAreas(fSourceConcave);

	if ((meshSource.vecDimSizes.size() == 0) ||
	    (meshSource.vecDimNames.size() == 0)
	) {
		_EXCEPTIONT("Source mesh missing vecDimSizes and vecDimNames");
	}

	m_mapRemap.InitializeSourceDimensions(
		meshSource.vecDimSizes,
		meshSource.vecDimNames);

	DataArray3D<int> dataSourceGLLNodes;
	DataArray3D<double> dataSourceGLLJacobian;

	if (meshSource.eDataLayout == Mesh::DataLayout_Volumetric) {
		m_mapRemap.SetSourceAreas(meshSource.vecFaceArea);
		m_mapRemap.InitializeSourceCoordinatesFromMeshFV(meshSource);

	} else {
		double dSourceNumericalArea =
			GenerateMetaData(
				meshSource,
				meshSource.nP,
				fBubble,
				dataSourceGLLNodes,
				dataSourceGLLJacobian);

		m_mapRemap.InitializeSourceCoordinatesFromMeshFE(
			meshSource, meshSource.nP, dataSourceGLLNodes);

		if (meshSource.eDataLayout == Mesh::DataLayout_SpectralElementGLL) {
			GenerateUniqueJacobian(
				dataSourceGLLNodes,
				dataSourceGLLJacobian,
				m_mapRemap.GetSourceAreas());

		} else {
			GenerateDiscontinuousJacobian(
				dataSourceGLLJacobian,
				m_mapRemap.GetSourceAreas());
		}

	}

	// Generate EdgeMap on target mesh if not available
	Mesh & meshTarget = pobjTargetGrid->GetMesh();
	meshTarget.ConstructReverseNodeArray();
	meshTarget.ConstructEdgeMap();
	double dTargetArea = meshTarget.CalculateFaceAreas(fTargetConcave);

	if ((meshTarget.vecDimSizes.size() == 0) ||
	    (meshTarget.vecDimNames.size() == 0)
	) {
		_EXCEPTIONT("Target mesh missing vecDimSizes and vecDimNames");
	}

	m_mapRemap.InitializeTargetDimensions(
		meshTarget.vecDimSizes,
		meshTarget.vecDimNames);

	DataArray3D<int> dataTargetGLLNodes;
	DataArray3D<double> dataTargetGLLJacobian;

 	if (meshTarget.eDataLayout == Mesh::DataLayout_Volumetric) {
		m_mapRemap.SetTargetAreas(meshTarget.vecFaceArea);
		m_mapRemap.InitializeTargetCoordinatesFromMeshFV(meshTarget);
	} else {
		double dTargetNumericalArea =
			GenerateMetaData(
				meshTarget,
				meshTarget.nP,
				fBubble,
				dataTargetGLLNodes,
				dataTargetGLLJacobian);

		m_mapRemap.InitializeTargetCoordinatesFromMeshFE(
			meshTarget, meshTarget.nP, dataTargetGLLNodes);

		if (meshTarget.eDataLayout == Mesh::DataLayout_SpectralElementGLL) {
			GenerateUniqueJacobian(
				dataTargetGLLNodes,
				dataTargetGLLJacobian,
				m_mapRemap.GetTargetAreas());

		} else {
			GenerateDiscontinuousJacobian(
				dataTargetGLLJacobian,
				m_mapRemap.GetTargetAreas());
		}
	}

	// Output mesh areas
	Announce("Source mesh area: %1.15e", dSourceArea);
	Announce("Target mesh area: %1.15e", dTargetArea);

	// Generate the overlap mesh
	Mesh * pMeshOverlap = new Mesh();

	// This can be parallelized over the faces of meshSource
	AnnounceStartBlock("Generate overlap mesh");
	if (dSourceArea < dTargetArea) {
		GenerateOverlapMesh_v2(
			meshSource,
			meshTarget,
			*pMeshOverlap,
			HRegrid::OverlapMeshMethod_Fuzzy);
	} else {
		GenerateOverlapMesh_v2(
			meshTarget,
			meshSource,
			*pMeshOverlap,
			HRegrid::OverlapMeshMethod_Fuzzy);
	}
	AnnounceEndBlock("Done");

	// Switch overlap mesh correspondence

	//pMeshOverlap->Write("ov_test.g");

	// Generate offline map
	AnnounceStartBlock("Generate offline map");

	if ((meshSource.eDataLayout == Mesh::DataLayout_SpectralElementGLL) ||
		(meshSource.eDataLayout == Mesh::DataLayout_DiscontinuousGLL)
	) {
		bool fContinuousIn =
			(meshSource.eDataLayout == Mesh::DataLayout_SpectralElementGLL);

		if (meshTarget.eDataLayout == Mesh::DataLayout_Volumetric) {

			LinearRemapSE4(
				meshSource,
				meshTarget,
				*pMeshOverlap,
				dataSourceGLLNodes,
				dataSourceGLLJacobian,
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

