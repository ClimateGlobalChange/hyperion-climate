///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridObject.cpp
///	\author  Paul Ullrich
///	\version March 16, 2017
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

#include "GridObject.h"

///////////////////////////////////////////////////////////////////////////////
// GridObjectConstructor
///////////////////////////////////////////////////////////////////////////////

std::string GridObjectConstructor::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	GridObject * pobjGrid = new GridObject("");
	if (pobjGrid == NULL) {
		_EXCEPTIONT("Unable to initialize GridObject");
	}

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjGrid;
	}

	// Initialize the grid with given parameters
	std::string strError =
		pobjGrid->Initialize(
			objreg,
			vecCommandLine,
			vecCommandLineType);

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjGrid;
	} else {
		delete pobjGrid;
	}

	return strError;
}

///////////////////////////////////////////////////////////////////////////////
// GridObject
///////////////////////////////////////////////////////////////////////////////

std::string GridObject::Initialize(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecFuncArguments,
	const std::vector<ObjectType> & vecFuncArgumentsType
) {
	// Create a blank grid
	if (vecFuncArguments.size() == 0) {
		return std::string("");
	}

	// At most two arguments expected
	if (vecFuncArguments.size() > 2) {
		return std::string("ERROR: Invalid arguments to grid()");
	}

	// Initialize from file name
	if (vecFuncArgumentsType[0] == ObjectType_String) {
		if (vecFuncArguments.size() > 2) {
			return std::string("ERROR: Too many arguments to grid()");
		}
		if ((vecFuncArguments.size() == 2) && 
		    (vecFuncArgumentsType[1] != ObjectType_String)
		) {
			return std::string("ERROR: Second argument to grid([filename]) "
				"must specify grid type");
		}

		// First argument specifies the grid filename
		m_mesh.Read(vecFuncArguments[0]);

		// Second argument specifies finite element type
		if (vecFuncArguments.size() == 2) {
			if (vecFuncArguments[1].substr(0,6) == "sem_np") {
				int nP = atoi(vecFuncArguments[1].substr(6).c_str());
				if (nP < 2) {
					return std::string("ERROR: sem_np# must be at least 2");
				}
				m_mesh.InitializeAsFiniteElement(
					Mesh::DataLayout_SpectralElementGLL,
					nP);

			} else if (vecFuncArguments[1].substr(0,8) == "dggll_np") {
				int nP = atoi(vecFuncArguments[1].substr(8).c_str());
				if (nP < 1) {
					return std::string("ERROR: dggll_np# must be at least 1");
				}
				m_mesh.InitializeAsFiniteElement(
					Mesh::DataLayout_DiscontinuousGLL,
					nP);

			} else {
				return std::string("ERROR: Second argument to grid() must be "
					"\"sem_np#\" or \"dggll_np#\"");
			}
		}

	// Initialize from parameter list
	} else if (vecFuncArgumentsType[0] == ObjectType_Token) {

		// Check for GridType
		Object * pobjParam =
			dynamic_cast<Object *>(
				objreg.GetObject(vecFuncArguments[0]));
		if (pobjParam == NULL) {
			return std::string("ERROR: Invalid first argument to grid(): "
				"Expected Token or String");
		}

		StringObject * pobjGridType =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobjParam->ChildName("gridtype")));
		if (pobjParam == NULL) {
			return std::string("ERROR: parameter_list argument must specify gridtype");
		}

		if (pobjGridType->Value() == "rll") {
			return InitializeRLLGrid(objreg, pobjParam);
		}

	// Invalid arguments
	} else {
		return std::string("ERROR: Invalid arguments to grid()");
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string GridObject::InitializeRLLGrid(
	const ObjectRegistry & objreg,
	const Object * pobjParameters
) {
	double dLatBegin = -90.0;
	double dLatEnd = 90.0;
	double dLonBegin = 0.0;
	double dLonEnd = 360.0;

	int nLatitudes = 0;
	int nLongitudes = 0;

	// Minimum latitude
	StringObject * pobjMinimumLatitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParameters->ChildName("minlat")));
	if (pobjMinimumLatitude != NULL) {
		bool fSuccessMag =
			pobjMinimumLatitude->ToUnit(
				"deg", &(dLatBegin));
		if (!fSuccessMag) {
			return std::string("Cannot convert ")
				+ pobjParameters->ChildName("minlat")
				+ std::string(" to \"deg\"");
		}
	}

	// Maximum latitude
	StringObject * pobjMaximumLatitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParameters->ChildName("maxlat")));
	if (pobjMaximumLatitude != NULL) {
		bool fSuccess =
			pobjMaximumLatitude->ToUnit(
				"deg", &(dLatEnd));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParameters->ChildName("maxlat")
				+ std::string(" to \"deg\"");
		}
	}

	// Minimum longitude
	StringObject * pobjMinimumLongitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParameters->ChildName("minlon")));
	if (pobjMinimumLongitude != NULL) {
		bool fSuccess =
			pobjMinimumLongitude->ToUnit(
				"deg", &(dLonBegin));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParameters->ChildName("minlon")
				+ std::string(" to \"deg\"");
		}
	}

	// Maximum longitude
	StringObject * pobjMaximumLongitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParameters->ChildName("maxlon")));
	if (pobjMaximumLongitude != NULL) {
		bool fSuccess =
			pobjMaximumLongitude->ToUnit(
				"deg", &(dLonEnd));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParameters->ChildName("maxlon")
				+ std::string(" to \"deg\"");
		}
	}

	// Number of latitudes
	IntegerObject * pobjLatitudeCount =
		dynamic_cast<IntegerObject *>(
			objreg.GetObject(pobjParameters->ChildName("nlat")));
	if (pobjLatitudeCount == NULL) {
		return std::string("ERROR: Missing mandatory parameter \"nlat\""
			" or invalid type in ")
				+ pobjParameters->Name();
	}
	nLatitudes = pobjLatitudeCount->Value();

	// Number of longitudes
	IntegerObject * pobjLongitudeCount =
		dynamic_cast<IntegerObject *>(
			objreg.GetObject(pobjParameters->ChildName("nlon")));
	if (pobjLongitudeCount == NULL) {
		return std::string("ERROR: Missing mandatory parameter \"nlon\""
			" or invalid type in ")
				+ pobjParameters->Name();
	}
	nLongitudes = pobjLongitudeCount->Value();

	// Initialize the Mesh
	m_mesh.InitializeAsRLL(
		dLatBegin,
		dLatEnd,
		dLonBegin,
		dLonEnd,
		nLatitudes,
		nLongitudes);

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string GridObject::Call(
	const ObjectRegistry & objreg,
	const std::string & strFunctionName,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	// Output information about the FileList to a CSV file
	if (strFunctionName == "append") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"append\"");
		}

		m_mesh.Read(vecCommandLine[0], true);

		return std::string("");
	}

	// Write the mesh to file
	if (strFunctionName == "write") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"write\"");
		}

		m_mesh.Write(vecCommandLine[0]);

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

