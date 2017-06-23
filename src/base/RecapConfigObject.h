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
#include "AccessMode.h"
#include "Object.h"
#include "FileListObject.h"
#include "GridObject.h"
#include "VariableLookupObject.h"
#include "Variable.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A GlobalFunction that builds a new RecapConfigObject.
///	</summary>
class RecapConfigObjectConstructor : public GlobalFunction {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RecapConfigObjectConstructor(const std::string & strName) :
		GlobalFunction(strName)
	{ }

public:
	///	<summary>
	///		Call a member function of this GlobalFunction.
	///	</summary>
	virtual std::string Call(
		const ObjectRegistry & objreg,
		const std::vector<std::string> & vecCommandLine,
		const std::vector<ObjectType> & vecCommandLineType,
		Object ** ppReturn
	);
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing a list of files.
///	</summary>
class RecapConfigObject : public Object {

friend class RecapConfigObjectConstructor;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RecapConfigObject(
		const std::string & strName
	) :
		m_eAccessMode(AccessMode_ReadOnly),
		m_varreg(this),
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
	///		Call a member function of this Object.
	///	</summary>
	virtual std::string Call(
		const ObjectRegistry & objreg,
		const std::string & strFunctionName,
		const std::vector<std::string> & vecCommandLine,
		const std::vector<ObjectType> & vecCommandLineType,
		Object ** ppReturn
	);

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
	///		Get a reference to the VariableRegistry.
	///	</summary>
	VariableRegistry & GetVariableRegistry() {
		return m_varreg;
	}

	///	<summary>
	///		Get the Variable with the specified name.
	///	</summary>
	std::string GetVariable(
		const std::string & strVariableName,
		Variable ** ppvar
	);

	///	<summary>
	///		Add a new variable from a template.
	///	</summary>
	std::string AddVariableFromTemplate(
		const RecapConfigObject * pobjRecapConfig,
		const Variable * pvar,
		Variable ** ppvarOut
	);

	///	<summary>
	///		Add a new variable from a template and replace the vertical dimension.
	///	</summary>
	std::string AddVariableFromTemplateWithNewVerticalDim(
		const RecapConfigObject * pobjRecapConfig,
		const Variable * pvar,
		const std::string & strVerticalDimName,
		Variable ** ppvarOut
	);

public:
	///	<summary>
	///		When adding a lock to a RecapConfigObject, also
	///		lock the file_list and grid.
	///	</summary>
	virtual void AddLock() {
		FileListObject * pobjFileList = GetFileList();
		GridObject * pobjGrid = GetGrid();

		pobjFileList->AddLock();
		pobjGrid->AddLock();
		Object::AddLock();
	}

	///	<summary>
	///		When releasing a lock to a RecapConfigObject, also
	///		lock the file_list and grid.
	///	</summary>
	virtual void ReleaseLock() {
		FileListObject * pobjFileList = GetFileList();
		GridObject * pobjGrid = GetGrid();

		pobjFileList->ReleaseLock();
		pobjGrid->ReleaseLock();
		Object::ReleaseLock();
	}

protected:
	///	<summary>
	///		RecapConfigObject only supports three child objects:
	///		grid, data and lookup, each of which are singletons.
	///	</summary>
	virtual bool AddChild(Object * pChild);

protected:
	///	<summary>
	///		Access mode to data in this configuration.
	///	</summary>
	AccessMode m_eAccessMode;

	///	<summary>
	///		Variable registry used in this configuration.
	///	</summary>
	VariableRegistry m_varreg;
};

///////////////////////////////////////////////////////////////////////////////

#endif

