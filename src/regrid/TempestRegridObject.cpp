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
#include "RecapConfigObject.h"
#include "GridObject.h"
#include "FileListObject.h"
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

void TempestRegridObject::PrepareDelete() {
	if (m_pobjSourceConfig != NULL) {
		m_pobjSourceConfig->ReleaseLock();
	}
	if (m_pobjTargetConfig != NULL) {
		m_pobjTargetConfig->ReleaseLock();
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string TempestRegridObject::Call(
	const ObjectRegistry & objreg,
	const std::string & strFunctionName,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	// Output the map as a netcdf file
	if (strFunctionName == "outputnetcdf") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"outputnetcdf\"");
		}

		m_mapRemap.Write(vecCommandLine[0]);

		return std::string("");

	}

	// Regrid data
	if (strFunctionName == "regrid") {
		if (vecCommandLineType.size() != 1) {
			return std::string("ERROR: Invalid parameters to function \"regrid\"");
		}

		// Variables
		std::vector<std::string> vecVariables;

		ListObject * pobjVariableList =
			dynamic_cast<ListObject *>(
				objreg.GetObject(vecCommandLine[0]));

		if (pobjVariableList == NULL) {
			StringObject * pobjVariableString =
				dynamic_cast<StringObject *>(
					objreg.GetObject(vecCommandLine[0]));
			if (pobjVariableString == NULL) {
				vecVariables.push_back(vecCommandLine[0]);
			} else {
				vecVariables.push_back(pobjVariableString->Value());
			}

		} else {
			for (size_t i = 0; i < pobjVariableList->Count(); i++) {
				StringObject * pobjVariableString =
					dynamic_cast<StringObject *>(
						objreg.GetObject(pobjVariableList->ChildName(i)));

				if (pobjVariableString == NULL) {
					return std::string("ERROR: First argument of regrid() must be"
						" a list of string objects");
				}

				vecVariables.push_back(pobjVariableString->Value());
			}
		}

		// Get FileListObjects
		FileListObject * pobjSourceFileList =
			m_pobjSourceConfig->GetFileList();
		FileListObject * pobjTargetFileList =
			m_pobjTargetConfig->GetFileList();

		// Check if a reduce target is available
		bool fReduceTargetAvailable =
			pobjTargetFileList->HasReduceTarget();

		// Check if times are compatible across FileLists
		bool fFileListCompatible =
			pobjSourceFileList->IsCompatible(pobjTargetFileList);

		// Check the size of all variables
		for (int v = 0; v < vecVariables.size(); v++) {
			Variable * pvarSource;
			std::string strError =
				m_pobjSourceConfig->GetVariable(
					vecVariables[v],
					&pvarSource);

			if (strError != "") {
				return strError;
			}

			if ((pvarSource->IsReductionOp()) && (!fReduceTargetAvailable)) {
				return std::string("ERROR: No reduction target available "
					"for output of variable reduction \"")
					+ pvarSource->Name() + std::string("\"");
			}
			if ((!pvarSource->IsReductionOp()) && (!fFileListCompatible)) {
				return std::string("ERROR: Source and target configurations "
					"are not compatible for regrid operation");
			}

			Variable * pvarTarget;
			strError =
				m_pobjTargetConfig->AddVariableFromTemplate(
					m_pobjSourceConfig,
					pvarSource,
					&pvarTarget);

			if (strError != "") {
				return strError;
			}

			// Apply the regridding operation to a reduction
			if (pvarSource->IsReductionOp()) {

				pvarSource->LoadGridData(
					Variable::SingleTimeIndex);

				pvarTarget->AllocateGridData(
					Variable::SingleTimeIndex);

				const DataArray1D<float> & dataSource =
					pvarSource->GetData();

				DataArray1D<float> & dataTarget =
					pvarTarget->GetData();

				m_mapRemap.ApplyFloat(
					dataSource,
					dataTarget);

				pvarTarget->WriteGridData(
					Variable::SingleTimeIndex);

			// Apply the regridding operation to all times
			} else {
				for (size_t t = 0; t < pobjSourceFileList->GetTimeCount(); t++) {
					pvarSource->LoadGridData(t);

					pvarTarget->AllocateGridData(t);

					const DataArray1D<float> & dataSource =
						pvarSource->GetData();

					DataArray1D<float> & dataTarget =
						pvarTarget->GetData();

					m_mapRemap.ApplyFloat(
						dataSource,
						dataTarget);

					pvarTarget->WriteGridData(t);
				}
			}
		}

/*
		for (size_t f = 0; f < sSourceFilenameCount; f++) {
			Announce("Input File:  %s", pobjSourceFileList->GetFilename(f).c_str());
			Announce("Output File: %s", pobjTargetFileList->GetFilename(f).c_str());
			m_mapRemap.Apply(
				pobjSourceFileList->GetFilename(f),
				pobjTargetFileList->GetFilename(f),
				vecVariables,
				std::string("ncol"));
		}
*/
		return std::string("");
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

	// Source configuration
	RecapConfigObject * pobjSourceConfig =
		dynamic_cast<RecapConfigObject *>(
			objreg.GetObject(vecFuncArguments[0]));
	if (pobjSourceConfig == NULL) {
		return std::string("ERROR: First argument of tempestregrid() must be"
			" of type recap_configuration");
	}
	if (!pobjSourceConfig->IsValid()) {
		return std::string("ERROR: recap_configuration must have both "
			" \"data\" and \"grid\" attributes before being used");
	}

	// Target grid
	RecapConfigObject * pobjTargetConfig =
		dynamic_cast<RecapConfigObject *>(
			objreg.GetObject(vecFuncArguments[1]));
	if (pobjTargetConfig == NULL) {
		return std::string("ERROR: Second argument of tempestregrid() must be"
			" of type recap_configuration");
	}
	if (!pobjTargetConfig->IsValid()) {
		return std::string("ERROR: recap_configuration must have both "
			" \"data\" and \"grid\" attributes before being used");
	}

	// Parameters
	Object * pobjParameters =
		dynamic_cast<Object *>(
			objreg.GetObject(vecFuncArguments[2]));
	if (pobjParameters == NULL) {
		return std::string("ERROR: Third argument to tempestregrid() must be"
			" of type parameter_list");
	}

	// Store pointer
	m_pobjSourceConfig = pobjSourceConfig;
	m_pobjTargetConfig = pobjTargetConfig;

	// Add locks to RecapConfigs
	pobjSourceConfig->AddLock();
	pobjTargetConfig->AddLock();

	// Get the grid
	GridObject * pobjSourceGrid = pobjSourceConfig->GetGrid();
	GridObject * pobjTargetGrid = pobjTargetConfig->GetGrid();

	// Pointers to Meshes
	Mesh * pmeshSource = &(pobjSourceGrid->GetMesh());
	Mesh * pmeshTarget = &(pobjTargetGrid->GetMesh());

	// Monotone type
	int nMonotoneType = 0;

	// No conservation
	bool fNoConservation = false;

	// Always bubble
	bool fBubble = true;

	// Concave faces in mesh
	bool fSourceConcave = pmeshSource->fConcave;

	// Concave faces in mesh
	bool fTargetConcave = pmeshTarget->fConcave;

	// Convexify the source mesh if it contains concave elements
	if (fSourceConcave) {
		Announce("Convexifying source mesh");
		pmeshSource = new Mesh;
		ConvexifyMesh(pobjSourceGrid->GetMesh(), *pmeshSource);
	}
	Mesh & meshSource = (*pmeshSource);

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
		if (meshSource.vecMultiFaceMap.size() != 0) {
			DataArray1D<double> vecSourceAreas(
				pobjSourceGrid->GetMesh().faces.size());

			for (size_t f = 0; f < meshSource.vecMultiFaceMap.size(); f++) {
				vecSourceAreas[meshSource.vecMultiFaceMap[f]] +=
					meshSource.vecFaceArea[f];
			}
			m_mapRemap.SetSourceAreas(vecSourceAreas);

		} else {
			m_mapRemap.SetSourceAreas(meshSource.vecFaceArea);
		}

		m_mapRemap.InitializeSourceCoordinatesFromMeshFV(
			pobjSourceGrid->GetMesh());

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
	if (fTargetConcave) {
		Announce("Convexifying target mesh");
		pmeshTarget = new Mesh;
		ConvexifyMesh(pobjTargetGrid->GetMesh(), *pmeshTarget);
	}
	Mesh & meshTarget = (*pmeshTarget);
	Announce("Target mesh has %i faces", meshTarget.faces.size());
	Announce("Target mesh has %i nodes", meshTarget.nodes.size());

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
		if (meshTarget.vecMultiFaceMap.size() != 0) {
			DataArray1D<double> vecTargetAreas(
				pobjTargetGrid->GetMesh().faces.size());
			for (size_t f = 0; f < meshTarget.vecMultiFaceMap.size(); f++) {
				vecTargetAreas[meshTarget.vecMultiFaceMap[f]] +=
					meshTarget.vecFaceArea[f];
			}
			m_mapRemap.SetTargetAreas(vecTargetAreas);

		} else {
			m_mapRemap.SetTargetAreas(meshTarget.vecFaceArea);
		}
		m_mapRemap.InitializeTargetCoordinatesFromMeshFV(
			pobjTargetGrid->GetMesh());
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

		pMeshOverlap->ExchangeFirstAndSecondMesh();
	}
	AnnounceEndBlock("Done");

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

	// Delete the convexified Mesh
	if (fSourceConcave) {
		delete pmeshSource;
	}
	if (fTargetConcave) {
		delete pmeshTarget;
	}

	// Delete the overlap mesh
	delete pMeshOverlap;

	AnnounceEndBlock("Done");

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

};

