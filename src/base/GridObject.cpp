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

std::string GridObject::Initialize(
	const std::vector<std::string> & vecFuncArguments,
	const std::vector<ObjectType> & vecFuncArgumentsType
) {
	// Check arguments
	if (vecFuncArguments.size() == 0) {
		return std::string("ERROR: grid filename argument missing");
	}
	if (vecFuncArguments.size() > 2) {
		return std::string("ERROR: Too many arguments to grid()");
	}
	if (vecFuncArgumentsType[0] != ObjectType_String) {
		return std::string("ERROR: First argument to grid() must be of type string");
	}
	if ((vecFuncArguments.size() == 2) && 
	    (vecFuncArgumentsType[1] != ObjectType_String)
	) {
		return std::string("ERROR: Second argument to grid() must be of type string");
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
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

