///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestRegridObject.h
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

#ifndef _TEMPESTREGRIDOBJECT_H_
#define _TEMPESTREGRIDOBJECT_H_

#include "Announce.h"
#include "Object.h"
#include "GridElements.h"
#include "OfflineMap.h"

///////////////////////////////////////////////////////////////////////////////

class RecapConfigObject;

///////////////////////////////////////////////////////////////////////////////

namespace HRegrid {

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A GlobalFunction that builds a new TempestRegridObject.
///	</summary>
class TempestRegridObjectConstructor : public GlobalFunction {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TempestRegridObjectConstructor(const std::string & strName) :
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
class TempestRegridObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TempestRegridObject(
		const std::string & strName
	) :
		Object(strName),
		m_pobjSourceConfig(NULL),
		m_pobjTargetConfig(NULL),
		m_pMeshOverlap(NULL)
	{ }

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~TempestRegridObject();

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new TempestRegridObject(strDuplicateName),
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

public:
	///	<summary>
	///		Initialize the Grid from command-line arguments.
	///	</summary>
	std::string Initialize(
		const ObjectRegistry & objreg,
		const std::vector<std::string> & vecFuncArguments,
		const std::vector<ObjectType> & vecFuncArgumentsType
	);

protected:
	///	<summary>
	///		Pointer to source configuration.
	///	</summary>
	RecapConfigObject * m_pobjSourceConfig;

	///	<summary>
	///		Pointer to target configuration.
	///	</summary>
	RecapConfigObject * m_pobjTargetConfig;

	///	<summary>
	///		The overlap Mesh used for building the offline map.
	///	</summary>
	Mesh * m_pMeshOverlap;

	///	<summary>
	///		The offline map.
	///	</summary>
	OfflineMap m_mapRemap;
};

///////////////////////////////////////////////////////////////////////////////

};

#endif

