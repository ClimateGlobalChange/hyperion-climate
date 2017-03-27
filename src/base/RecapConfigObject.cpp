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
			return std::string("ERROR: Invalid arguments to function \"set_reduce_target\"");
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

