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

#include <cfloat>

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
	std::string strError;

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

		strError =
			ArgumentToStringVector(
				objreg, vecCommandLine[0], vecVariables);

		if (strError != "") {
			return strError;
		}

		// Perform regridding
		strError = Regrid(vecVariables);
		if (strError != "") {
			return strError;
		}

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
	std::string strError;

	// Clear level values
	m_strLevelsUnits = "";
	m_vecLevelsValues.clear();

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

	// Check if vertical remapping is requested
	Object * pobjLevels = pobjParameters->GetChild("levels");
	if (pobjLevels != NULL) {
		std::vector<std::string> vecLevels;

		// Get a list of output levels
		strError = ArgumentToStringVector(objreg, pobjLevels->Name(), vecLevels);
		if (strError != "") {
			return strError;
		}

		if (vecLevels.size() == 0) {
			return std::string("ERROR: levels must contain at least one entry");
		}

		// Check if levels are specified in height units
		double dValue;
		bool fIsHeightUnit =
			StringToValueUnit(
				vecLevels[0],
				std::string("m"),
				dValue,
				false);

		// Height units
		if (fIsHeightUnit) {

			// Check for geopotential height at surface
			Variable * pvarPHIS = NULL;
			pobjSourceConfig->GetVariable("PHIS", &pvarPHIS);
			if (pvarPHIS == NULL) {
				return std::string("ERROR: A surface height variable cannot be found"
					" but is needed for vertical interpolation");
			}

			// Check for geopotential height variable
			Variable * pvarZ3 = NULL;
			pobjSourceConfig->GetVariable("Z3", &pvarZ3);
			if (pvarZ3 != NULL) {
				m_eVertInterpMethod = VerticalInterp_PHISZ3toZ;
			}

			// Check for surface pressure and temperature
			Variable * pvarPS = NULL;
			Variable * pvarT = NULL;
			pobjSourceConfig->GetVariable("PS", &pvarPS);
			pobjSourceConfig->GetVariable("T", &pvarT);
			if ((pvarPS != NULL) && (pvarT != NULL)) {
				m_eVertInterpMethod = VerticalInterp_PHISPSTtoZ;
			}
		}

		// Extract remaining levels
		bool fExtractSuccess =
			ExtractValueUnit(
				vecLevels[0],
				dValue,
				m_strLevelsUnits);

		if (!fExtractSuccess) {
			_EXCEPTIONT("Logic error");
		}

		m_vecLevelsValues.push_back(dValue);

		for (int k = 1; k < vecLevels.size(); k++) {
			fExtractSuccess =
				StringToValueUnit(
					vecLevels[k],
					m_strLevelsUnits,
					dValue,
					false);

			if (!fExtractSuccess) {
				return std::string("ERROR: All levels must have the same units");
			}

			m_vecLevelsValues.push_back(dValue);
		}
/*
		std::cout << "Unit: " << m_strLevelsUnits << std::endl;
		for (int k = 0; k < m_vecLevelsValues.size(); k++) {
			std::cout << "Level " << k << ": " << m_vecLevelsValues[k] << std::endl;
		}

		_EXCEPTION();
*/
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

	// Check if the two grids are equal to one another

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

std::string TempestRegridObject::Regrid(
	const std::vector<std::string> & vecVariables
) {
	std::string strError;

	// Get FileListObjects
	FileListObject * pobjSourceFileList =
		m_pobjSourceConfig->GetFileList();
	FileListObject * pobjTargetFileList =
		m_pobjTargetConfig->GetFileList();

	// Get GridObject
	GridObject * pobjSourceGrid =
		m_pobjSourceConfig->GetGrid();

	// Check if a reduce target is available
	bool fReduceTargetAvailable =
		pobjTargetFileList->HasReduceTarget();

	// Check if times are compatible across FileLists
	bool fFileListCompatible =
		pobjSourceFileList->IsCompatible(pobjTargetFileList);

	// Check the size of all variables
	for (int v = 0; v < vecVariables.size(); v++) {

		bool fInterpolateVertical = (m_vecLevelsValues.size() != 0);

		// Source variable
		Variable * pvarSource;
		strError =
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

		if ((pvarSource->IsReductionOp()) && (m_vecLevelsValues.size() != 0)) {
			_EXCEPTIONT("Vertical interpolation on reduction operators not yet supported");
		}

		if (pvarSource->VerticalDimIx() == (-1)) {
			fInterpolateVertical = false;
		}

		// Perform vertical interpolation
		if (fInterpolateVertical) {

			// Initialize reduced iterators
			VariableAuxIndexIterator iterAux = pvarSource->GetAuxIndexBegin();
			VariableAuxIndexIterator iterEnd = pvarSource->GetAuxIndexEnd();
			//std::cout << iterAux.size() << " / " << pvarSource->VerticalDimIx() << std::endl;
			//std::cout << iterEnd.size() << " / " << pvarSource->VerticalDimIx() << std::endl;
			iterAux.RemoveDim(pvarSource->VerticalDimIx());
			iterEnd.RemoveDim(pvarSource->VerticalDimIx());

			// Target variable
			Variable * pvarTarget = NULL;

			// Get PHIS, if needed
			Variable * pvarPHIS = NULL;
			if ((m_eVertInterpMethod == VerticalInterp_PHISZ3toZ) ||
				(m_eVertInterpMethod == VerticalInterp_PHISPSTtoZ)
			) {
				strError = m_pobjSourceConfig->GetVariable("PHIS", &pvarPHIS);
				if (strError != "") return strError;

				strError = pobjTargetFileList->AddVerticalDimension(
					std::string("z"),
					m_vecLevelsValues,
					m_strLevelsUnits);
				if (strError != "") return strError;

				strError =
					m_pobjTargetConfig->AddVariableFromTemplateWithNewVerticalDim(
						m_pobjSourceConfig,
						pvarSource,
						std::string("z"),
						&pvarTarget);
				if (strError != "") return strError;

			} else {
				return std::string("ERROR: This type of vertical interpolation is not implemented");
			}

			// Interpolate vertically using Z3
			if (m_eVertInterpMethod == VerticalInterp_PHISZ3toZ) {

				bool fReloaded = false;

				// Get Z3 variable
				Variable * pvarZ3;
				strError = m_pobjSourceConfig->GetVariable("Z3", &pvarZ3);
				if (strError != "") {
					return strError;
				}

				if (pvarZ3->VerticalDimSize() < 2) {
					return std::string("ERROR: Z3 must have at least two vertical levels");
				}

				// Allow variables to load more than one data instance
				pvarZ3->RemoveSingleDataInstanceLimitation();

				// Loop through all auxiliary indices
				for (; iterAux != iterEnd; iterAux++) {

					// Create new auxiliary index with vertical level
					VariableAuxIndex ixAux = iterAux;
					ixAux.insert(ixAux.begin() + pvarSource->VerticalDimIx(), 0);

					// Set of vertical levels required for interpolation
					std::set<long> setVerticalLevelsUsed;

					// Load PHIS data
					DataArray1D<float> * pdataPHIS = NULL;
					strError = pvarPHIS->LoadGridData(iterAux, &pdataPHIS);
					if (strError != "") {
						return strError;
					}

					// Unload existing Z3 data
					pvarZ3->UnloadAllGridData();

					// Data instances storing data on bounding levels
					DataArray1D<float> * pdataZ3Last = NULL;
					DataArray1D<float> * pdataZ3Next = NULL;

					// Load Z3 from the lowest model level
					ixAux[pvarZ3->VerticalDimIx()] = 0;
					strError = pvarZ3->LoadGridData(ixAux, &pdataZ3Next, &fReloaded);
					if (strError != "") return strError;

					// Remove surface component from Z3
					if (fReloaded) {
						for (int i = 0; i < pdataZ3Next->GetRows(); i++) {
							(*pdataZ3Next)[i] -= (*pdataPHIS)[i] / 9.80616;
						}
					}

					// Get variable order
					double dOrder = pvarZ3->VerticalDimOrder();

					// Determine which vertical levels are needed from the source data
					std::set<int> setVerticalLevelsNeeded;

					for (long l = -1; l < pvarZ3->VerticalDimSize(); l++) {

						// Load in Z3 data
						pdataZ3Last = pdataZ3Next;
						if (l != pvarZ3->VerticalDimSize()-1) {
							ixAux[pvarSource->VerticalDimIx()] = l+1;
							strError = pvarZ3->LoadGridData(ixAux, &pdataZ3Next, &fReloaded);
							if (strError != "") return strError;

							// Remove surface component from Z3
							if (fReloaded) {
								for (int i = 0; i < pdataZ3Next->GetRows(); i++) {
									(*pdataZ3Next)[i] -= (*pdataPHIS)[i] / 9.80616;
								}
							}
						}

						// Loop through all output levels
						for (size_t k = 0; k < m_vecLevelsValues.size(); k++) {

							// Loop through all cells on this level
							for (int i = 0; i < pdataZ3Last->GetRows(); i++) {
								double dLastZ;
								double dNextZ;

								if (l == -1) {
									dLastZ = dOrder * (-DBL_MAX);
								} else {
									dLastZ = (*pdataZ3Last)[i];
								}

								if (l == pvarZ3->VerticalDimSize()-1) {
									dNextZ = dOrder * (+DBL_MAX);
								} else {
									dNextZ = (*pdataZ3Next)[i];
								}

								// Level value found in interval [dLastZ, dNextZ)
								if ((dOrder * m_vecLevelsValues[k] >= dOrder * dLastZ) &&
									(dOrder * m_vecLevelsValues[k] < dOrder * dNextZ)
								) {
									if (l == -1) {
										setVerticalLevelsNeeded.insert(0);
									} else if (l == pvarZ3->VerticalDimSize()-1) {
										setVerticalLevelsNeeded.insert(l);
									} else {
										setVerticalLevelsNeeded.insert(l);
										setVerticalLevelsNeeded.insert(l+1);
									}
									break;
								}
							}
						}
					}

					// Buffer data vector, storing remapped data on levels
					// where data is needed for vertical interpolation
					std::vector< DataArray1D<float> > vecdataInterp;
					vecdataInterp.resize(pvarZ3->VerticalDimSize());

					// Horizontally remap data on needed levels
					for (
						std::set<int>::iterator iterLevNeeded = setVerticalLevelsNeeded.begin();
						iterLevNeeded != setVerticalLevelsNeeded.end();
						iterLevNeeded++
					) {
						// Allocate buffer data
						vecdataInterp[*iterLevNeeded].Allocate(m_mapRemap.GetTargetSize());

						// Load source data
						DataArray1D<float> * pdataSource = NULL;
						ixAux[pvarSource->VerticalDimIx()] = *iterLevNeeded;
						strError = pvarSource->LoadGridData(ixAux, &pdataSource);
						if (strError != "") return strError;

						// Apply the horizontal map on this level
						m_mapRemap.ApplyFloat(
							(*pdataSource),
							vecdataInterp[*iterLevNeeded]);
					}

					// Loop through all output levels
					for (size_t k = 0; k < m_vecLevelsValues.size(); k++) {

						// Load in target data array
						DataArray1D<float> * pdataTarget = NULL;
						ixAux[pvarTarget->VerticalDimIx()] = k;
						strError = pvarTarget->AllocateGridData(ixAux, &pdataTarget);
						if (strError != "") return strError;

						// Perform the vertical interpolation
						DataArray1D<float> * pdataZ3Last = NULL;
						DataArray1D<float> * pdataZ3Next = NULL;

						// Load Z3 from the lowest model level
						ixAux[pvarZ3->VerticalDimIx()] = 0;
						strError = pvarZ3->LoadGridData(ixAux, &pdataZ3Next, &fReloaded);
						if (strError != "") return strError;

						// Remove surface component from Z3
						if (fReloaded) {
							for (int i = 0; i < pdataZ3Next->GetRows(); i++) {
								(*pdataZ3Next)[i] -= (*pdataPHIS)[i] / 9.80616;
							}
						}

						// Loop through all input levels
						for (long l = -1; l < pvarZ3->VerticalDimSize(); l++) {
							pdataZ3Last = pdataZ3Next;
							if (l != pvarZ3->VerticalDimSize()-1) {
								ixAux[pvarSource->VerticalDimIx()] = l+1;
								pvarZ3->LoadGridData(ixAux, &pdataZ3Next, &fReloaded);

								// Remove surface component from Z3
								if (fReloaded) {
									for (int i = 0; i < pdataZ3Next->GetRows(); i++) {
										(*pdataZ3Next)[i] -= (*pdataPHIS)[i] / 9.80616;
									}
								}
							}

							// Loop through all cells
							for (int i = 0; i < pdataTarget->GetRows(); i++) {
								double dLastZ;
								double dNextZ;

								if (l == -1) {
									dLastZ = dOrder * (-DBL_MAX);
								} else {
									dLastZ = (*pdataZ3Last)[i];
								}

								if (l == pvarZ3->VerticalDimSize()-1) {
									dNextZ = dOrder * (+DBL_MAX);
								} else {
									dNextZ = (*pdataZ3Next)[i];
								}

								// Level value found in interval [dLastZ, dNextZ)
								if ((dOrder * m_vecLevelsValues[k] >= dOrder * dLastZ) &&
									(dOrder * m_vecLevelsValues[k] < dOrder * dNextZ)
								) {
									if (l == -1) {
										(*pdataTarget)[i] = vecdataInterp[0][i];
									} else if (l == pvarZ3->VerticalDimSize()-1) {
										(*pdataTarget)[i] = vecdataInterp[l][i];
									} else {
										double dAlpha =
											(m_vecLevelsValues[k] - (*pdataZ3Last)[i])
											/ ((*pdataZ3Next)[i] - (*pdataZ3Last)[i]);
										(*pdataTarget)[i] =
											dAlpha * vecdataInterp[l+1][i]
											+ (1.0 - dAlpha) * vecdataInterp[l][i];
									}
								}
							}
						}

						// Write target data
						ixAux[pvarTarget->VerticalDimIx()] = k;
						strError = pvarTarget->WriteGridData(ixAux);
						if (strError != "") return strError;

						pvarTarget->UnloadAllGridData();
					}

					// Unload all data
					pvarZ3->UnloadAllGridData();
					pvarSource->UnloadAllGridData();
				}
			}

		// No vertical interpolation
		} else {

			// Add the output variable to the target configuration
			Variable * pvarTarget;
			strError =
				m_pobjTargetConfig->AddVariableFromTemplate(
					m_pobjSourceConfig,
					pvarSource,
					&pvarTarget);

			if (strError != "") {
				return strError;
			}

			// Loop over all auxiliary indices and remap data
			//std::cout << pvarSource->GetAuxIndexBegin().ToString() << std::endl;
			//std::cout << pvarSource->GetAuxIndexEnd().ToString() << std::endl;
			VariableAuxIndexIterator iterAux = pvarSource->GetAuxIndexBegin();
			for (; iterAux != pvarSource->GetAuxIndexEnd(); iterAux++) {

				//std::cout << iterAux.ToString() << std::endl;

				// Load source data
				DataArray1D<float> * pdataSource = NULL;
				strError = pvarTarget->LoadGridData(iterAux, &pdataSource);
				if (strError != "") return strError;

				// Allocate target data
				DataArray1D<float> * pdataTarget = NULL;
				strError = pvarTarget->AllocateGridData(iterAux, &pdataTarget);
				if (strError != "") return strError;

				// Apply the map
				m_mapRemap.ApplyFloat(
					(*pdataSource),
					(*pdataTarget));

				// Write target data
				pvarTarget->WriteGridData(iterAux);
			}
		}
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

};

