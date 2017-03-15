///////////////////////////////////////////////////////////////////////////////
///
///	\file    RecapConfigObject.h
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

#ifndef _RECAPCONFIGOBJECT_H_
#define _RECAPCONFIGOBJECT_H_

#include "Announce.h"
#include "Object.h"
#include "FileListObject.h"
#include "GridObject.h"
#include "VariableLookupObject.h"
#include "Variable.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing a list of files.
///	</summary>
class RecapConfigObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RecapConfigObject(
		const std::string & strName
	) :
		Object(strName)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new RecapConfigObject(strDuplicateName),
			objreg);
	}

	///	<summary>
	///		Check validity of configuration.
	///	</summary>
	bool IsValid() const {
		if ((GetChild("data") != NULL) &&
		    (GetChild("grid") != NULL)
		) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Get file_list information from configuration.
	///	</summary>
	FileListObject * GetFileList() const {
		FileListObject * pobj =
			 dynamic_cast<FileListObject *>(GetChild("data"));
		if (pobj == NULL) {
			_EXCEPTION();
		}
		return pobj;
	}

	///	<summary>
	///		Get grid information from configuration.
	///	</summary>
	GridObject * GetGrid() const {
		GridObject * pobj =
			 dynamic_cast<GridObject *>(GetChild("grid"));
		if (pobj == NULL) {
			_EXCEPTION();
		}
		return pobj;
	}

	///	<summary>
	///		Get variable lookup information from configuration.
	///	</summary>
	VariableLookupObject * GetVariableLookup() const {
		VariableLookupObject * pobj =
			 dynamic_cast<VariableLookupObject *>(GetChild("lookup"));
		return pobj;
	}

	///	<summary>
	///		Get the Variable with the specified name.
	///	</summary>
	std::string GetVariable(
		const std::string & strVariableName,
		Variable ** pVariable
	) {
		// Convert the name from native to standard
		std::string strStandardName = strVariableName;
		VariableLookupObject * pobjLookup = GetVariableLookup();
		if (pobjLookup != NULL) {
			strStandardName = pobjLookup->Convert(strStandardName);
		}

		// Find the Variable in the VariableRegistry
		if (pVariable == NULL) {
			_EXCEPTIONT("Invalid argument");
		}

		(*pVariable) =
			varreg.FindOrRegister(
				strStandardName,
				GetFileList());

		// Variable could not be found or registered
		if ((*pVariable) == NULL) {
			return std::string("Invalid variable \"")
				+ strVariableName + std::string("\"");
		}
		return std::string("");
	}

protected:
	///	<summary>
	///		RecapConfigObject only supports three child objects.
	///	</summary>
	virtual bool AddChild(Object * pChild) {
		if (pChild == NULL) {
			_EXCEPTIONT("Invalid child pointer");
		}
		if (pChild->Name() == ChildName("data")) {
			FileListObject * pobj = dynamic_cast<FileListObject *>(pChild);
			if (pobj == NULL) {
				Announce("ERROR: %s must be of type file_list",
					pChild->Name().c_str());
				return false;
			}
			return Object::AddChild(pChild);
		}
		if (pChild->Name() == ChildName("grid")) {
			GridObject * pobj = dynamic_cast<GridObject *>(pChild);
			if (pobj == NULL) {
				Announce("ERROR: %s must be of type grid",
					pChild->Name().c_str());
				return false;
			}
			return Object::AddChild(pChild);
		}
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

		Announce("ERROR: Only \"data\", \"grid\" and \"lookup\" allowed as "
			"children of recap_configuration()");

		return false;
	}

protected:
	///	<summary>
	///		Variable registry used in this configuration.
	///	</summary>
	VariableRegistry varreg;
};

///////////////////////////////////////////////////////////////////////////////

#endif

