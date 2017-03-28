///////////////////////////////////////////////////////////////////////////////
///
///	\file    RecapConfigObject.cpp
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

#include "RecapConfigObject.h"

///////////////////////////////////////////////////////////////////////////////
// RecapConfigObjectConstructor
///////////////////////////////////////////////////////////////////////////////

std::string RecapConfigObjectConstructor::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	RecapConfigObject * pobjRecapConfig = new RecapConfigObject("");
	if (pobjRecapConfig == NULL) {
		_EXCEPTIONT("Unable to initialize RecapConfigObject");
	}

	// Check command line
	if (vecCommandLine.size() > 1) {
		return std::string("ERROR: Too many arguments to recap_configuration()");
	}
	if (vecCommandLine.size() == 1) {
		if (vecCommandLineType[0] != ObjectType_String) {
			return std::string("ERROR: First argument to recap_configuration() "
				"must be either \"read\", \"write\" or \"create\"");
		}

		if (vecCommandLine[0] == "read") {
			pobjRecapConfig->m_eAccessMode = AccessMode_ReadOnly;

		} else if (vecCommandLine[0] == "write") {
			pobjRecapConfig->m_eAccessMode = AccessMode_ReadWrite;

		} else if (vecCommandLine[0] == "create") {
			pobjRecapConfig->m_eAccessMode = AccessMode_Create;

		} else {
			return std::string("ERROR: First argument to recap_configuration() "
				"must be either \"read\", \"write\" or \"create\"");
		}
	}

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjRecapConfig;
	} else {
		delete pobjRecapConfig;
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////
// RecapConfigObject
///////////////////////////////////////////////////////////////////////////////

std::string RecapConfigObject::Call(
	const ObjectRegistry & objreg,
	const std::string & strFunctionName,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	if (!IsValid()) {
		return std::string("ERROR: recap_configuration must have both "
			" \"data\" and \"grid\" attributes before being used");
	}

	// Add a new file to the FileList for reductions
	if (strFunctionName == "set_reduce_target") {

		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid arguments to function set_reduce_target()");
		}

		if (m_eAccessMode == AccessMode_ReadOnly) {
			return std::string("ERROR: set_reduce_target() requires "
				"write-permissions in recap_configuration");
		}

		FileListObject * pobjFileList = GetFileList();

		std::string strError =
			pobjFileList->CreateFileNoTime(
				vecCommandLine[0],
				GetGrid());

		if (strError != "") {
			return strError;
		}

		return
			pobjFileList->SetReduceTarget(
				vecCommandLine[0]);

	}
/*
	// Add a new file to the FileList for single time slices
	if (strFunctionName == "add_data_file") {
		if (!IsValid()) {
			return std::string("ERROR: recap_configuration must have both "
				" \"data\" and \"grid\" attributes before being used");
		}

		if ((vecCommandLineType.size() != 2) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"add_file\"");
		}

		StringObject * pobjFilename =
			objreg.GetObject(vecCommandLine[0]);
		if (pobjFilename == NULL) {
			return std::string("ERROR: First argument to add_file() "
				"must be of type string");
		}

		// Variables
		std::vector<std::string> vecVariables;

		ListObject * pobjVariableList =
			dynamic_cast<ListObject *>(
				objreg.GetObject(vecCommandLine[1]));

		if (pobjVariableList == NULL) {
			StringObject * pobjVariableString =
				dynamic_cast<StringObject *>(
					objreg.GetObject(vecCommandLine[1]));
			if (pobjVariableString == NULL) {
				return std::string("ERROR: Second argument of add_file() must be"
					" a list of strings or singular string");
			}

			vecVariables.push_back(pobjVariableString->Value());

		} else {
			for (size_t i = 0; i < pobjVariableList->Count(); i++) {
				StringObject * pobjVariableString =
					dynamic_cast<StringObject *>(
						objreg.GetObject(pobjVariableList->ChildName(i)));

				if (pobjVariableString == NULL) {
					return std::string("ERROR: Second argument of add_file() must be"
						" a list of strings or singular string");
				}

				vecVariables.push_back(pobjVariableString->Value());
			}
		}

		return
			GetFileList()->AddDataFile(
				pobjFilename->Value(),
				GetGridObject(),
				vecVariables);
	}
*/
	return
		Object::Call(
			objreg,
			strFunctionName,
			vecCommandLine,
			vecCommandLineType,
			ppReturn);
}

///////////////////////////////////////////////////////////////////////////////

bool RecapConfigObject::AddChild(
	Object * pChild
) {
	if (pChild == NULL) {
		_EXCEPTIONT("Invalid child pointer");
	}

	// When the GridObject is added, update the FileListObject
	if (pChild->Name() == ChildName("grid")) {
		GridObject * pobjGrid = dynamic_cast<GridObject *>(pChild);
		if (pobjGrid == NULL) {
			Announce("ERROR: %s must be of type grid",
				pChild->Name().c_str());
			return false;
		}

		FileListObject * pobjFileList =
			 dynamic_cast<FileListObject *>(GetChild("data"));

		if (pobjFileList != NULL) {
			const Mesh & mesh = pobjGrid->GetMesh();

			for (int d = 0; d < mesh.vecDimNames.size(); d++) {
				pobjFileList->AddDimension(
					mesh.vecDimNames[d],
					mesh.vecDimSizes[d],
					true);
			}
		}
		return Object::AddChild(pChild);
	}

	// When the FileListObject is added, update using the GridObject
	if (pChild->Name() == ChildName("data")) {
		FileListObject * pobjFileList = dynamic_cast<FileListObject *>(pChild);
		if (pobjFileList == NULL) {
			Announce("ERROR: %s must be of type file_list",
				pChild->Name().c_str());
			return false;
		}

		pobjFileList->SetRecapConfigObject(this);

		GridObject * pobjGrid =
			dynamic_cast<GridObject *>(GetChild("grid"));

		if (pobjGrid != NULL) {
			const Mesh & mesh = pobjGrid->GetMesh();

			for (int d = 0; d < mesh.vecDimNames.size(); d++) {
				pobjFileList->AddDimension(
					mesh.vecDimNames[d],
					mesh.vecDimSizes[d],
					true);
			}
		}
		return Object::AddChild(pChild);
	}

	// Variable lookup table
	if (pChild->Name() == ChildName("lookup")) {
		VariableLookupObject * pobj =
			dynamic_cast<VariableLookupObject *>(pChild);
		if (pobj == NULL) {
			Announce("ERROR: %s must be of type variable_lookup",
				pChild->Name().c_str());
			return false;
		}
		return Object::AddChild(pChild);
	}

	// Invalid child object
	Announce("ERROR: Only \"data\", \"grid\" and \"lookup\" allowed as "
		"children of recap_configuration()");

	return false;
}

///////////////////////////////////////////////////////////////////////////////

std::string RecapConfigObject::GetVariable(
	const std::string & strVariableName,
	Variable ** ppvar
) {
	// Convert the name from native to standard
	std::string strStandardName = strVariableName;
	VariableLookupObject * pobjLookup = GetVariableLookup();
	if (pobjLookup != NULL) {
		strStandardName = pobjLookup->Convert(strStandardName);
	}

	// Find the Variable in the VariableRegistry
	if (ppvar == NULL) {
		_EXCEPTIONT("Invalid argument");
	}

	std::string strError =
		m_varreg.FindOrRegister(
			strStandardName,
			ppvar);

	// Variable could not be found or registered
	if (strError != "") {
		return strError;
	}
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string RecapConfigObject::AddVariableFromTemplate(
	const RecapConfigObject * pobjSourceRecapConfig,
	const Variable * pvar,
	Variable ** ppvarOut
) {
	if (!IsValid()) {
		_EXCEPTIONT("Invalid RecapConfigObject");
	}
	if (pobjSourceRecapConfig == NULL) {
		_EXCEPTIONT("Invalid RecapConfigObject");
	}
	if (!pobjSourceRecapConfig->IsValid()) {
		_EXCEPTIONT("Invalid RecapConfigObject");
	}
	if (pvar == NULL) {
		_EXCEPTIONT("Invalid Variable template");
	}

	// Add the variable to the FileList
	VariableInfo * pvarinfo = NULL;

	std::string strError =
		GetFileList()->AddVariableFromTemplate(
			pobjSourceRecapConfig->GetFileList(),
			pvar,
			&pvarinfo);

	if (strError != "") {
		return strError;
	}

	// Add the variable to the VariableRegistry
	return
		m_varreg.FindOrRegister(
			pvar->Name(),
			ppvarOut);
}

///////////////////////////////////////////////////////////////////////////////

