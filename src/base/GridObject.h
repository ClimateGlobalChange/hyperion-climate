///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridObject.h
///	\author  Paul Ullrich
///	\version March 10, 2017
///
///	<remarks>
///		Copyright 2000- Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _GRIDOBJECT_H_
#define _GRIDOBJECT_H_

#include "Announce.h"
#include "Object.h"
#include "DataArray1D.h"
#include "GridElements.h"
#include "GlobalFunction.h"

#include <vector>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A GlobalFunction that builds a new GridObject.
///	</summary>
class GridObjectConstructor : public GlobalFunction {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridObjectConstructor(const std::string & strName) :
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
///		A data structure describing the grid, including coordinates of
///		each data point and graph connectivity of elements.
///	</summary>
class GridObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	GridObject(
		const std::string & strName
	) :
		Object(strName)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	GridObject(
		const std::string & strName,
		const std::string & strMeshFile
	) :
		Object(strName),
		m_mesh(strMeshFile)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		Announce("ERROR: Grid objects cannot be duplicated");
		return (NULL);
	}

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
	///		Initialize a RLL Grid.
	///	</summary>
	std::string InitializeRLLGrid(
		const ObjectRegistry & objreg,
		const Object * pobjParameters
	);

public:
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
	///		Get a reference to the Mesh.
	///	</summary>
	const Mesh & GetMesh() const {
		return m_mesh;
	}

	///	<summary>
	///		Get a reference to the Mesh.
	///	</summary>
	Mesh & GetMesh() {
		return m_mesh;
	}

protected:
	///	<summary>
	///		Mesh associated with this Grid.
	///	</summary>
	Mesh m_mesh;

public:
	///	<summary>
	///		Grid dimensions.
	///	</summary>
	std::vector<size_t> m_nGridDim;
};

///////////////////////////////////////////////////////////////////////////////

#endif

