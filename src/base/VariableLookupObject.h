///////////////////////////////////////////////////////////////////////////////
///
///	\file    VariableLookupObject.h
///	\author  Paul Ullrich
///	\version March 14, 2017
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

#ifndef _VARIABLELOOKUPOBJECT_H_
#define _VARIABLELOOKUPOBJECT_H_

#include "Announce.h"
#include "Object.h"
#include "GlobalFunction.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A GlobalFunction that builds a new VariableLookupObject.
///	</summary>
class VariableLookupObjectConstructor : public GlobalFunction {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableLookupObjectConstructor(const std::string & strName) :
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
///		A data structure describing a lookup table for variable names.
///	</summary>
class VariableLookupObject : public Object {

public:
	///	<summary>
	///		A structure describing the content of the lookup table.
	///	</summary>
	struct LookupInfo {
		std::string strTargetName;
		std::string strUnits;
	};

	///	<summary>
	///		The lookup table map.
	///	</summary>
	typedef std::map<std::string, LookupInfo> LookupTableMap;

	///	<summary>
	///		The value_type of the lookup table map.
	///	</summary>
	typedef LookupTableMap::value_type LookupTableMapPair;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableLookupObject(
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
			new VariableLookupObject(strDuplicateName),
			objreg);
	}

public:
	///	<summary>
	///		Number of entries in lookup table.
	///	</summary>
	size_t LookupEntryCount() const {
		return m_mapLookupTable.size();	
	}

	///	<summary>
	///		Convert a standard variable string to a native variable string.
	///	</summary>
	std::string Convert(const std::string & strNative) {
		LookupTableMap::iterator iter = m_mapLookupTable.find(strNative);
		if (iter == m_mapLookupTable.end()) {
			return strNative;
		} else {
			return iter->second.strTargetName;
		}
	}

public:
	///	<summary>
	///		Populate the lookup table from a file.
	///	</summary>
	///	<returns>
	///		An error message if loading / parsing the lookup table fails.
	///	</returns>
	std::string PopulateFromFile(
		const std::string & strFilename
	);

protected:
	///	<summary>
	///		The lookup table.
	///	</summary>
	LookupTableMap m_mapLookupTable;
};

///////////////////////////////////////////////////////////////////////////////

#endif

