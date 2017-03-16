///////////////////////////////////////////////////////////////////////////////
///
///	\file    Variable.h
///	\author  Paul Ullrich
///	\version November 18, 2015
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _VARIABLE_H_
#define _VARIABLE_H_

#include "netcdfcpp.h"

#include "DataArray1D.h"
#include "GridObject.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class RecapConfigObject;

class FileListObject;

class VariableInfo;

class Variable;

typedef std::vector<Variable> VariableVector;

typedef int VariableIndex;

typedef std::vector<VariableIndex> VariableIndexVector;

///////////////////////////////////////////////////////////////////////////////

class VariableRegistry {

public:
	///	<summary>
	///		Find or register a variable.
	///	</summary>
	Variable * FindOrRegister(
		const std::string & strVariableName,
		const FileListObject * pobjFileList
	);

	///	<summary>
	///		Register a variable.  Or return an index if the Variable already
	///		exists in the registry.
	///	</summary>
	bool FindOrRegister(
		const Variable & var
	);

	///	<summary>
	///		Get the variable with the specified name.
	///	</summary>
	Variable * Get(
		const std::string & strName
	);

	///	<summary>
	///		Unload all data.
	///	</summary>
	void UnloadAllGridData();

private:
	///	<summary>
	///		Set of Variables.
	///	</summary>
	std::map<std::string, Variable> m_mapVariables;

	///	<summary>
	///		Array of variables.
	///	</summary>
	VariableVector m_vecVariables;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a parsed variable name.
///	</summary>
class Variable {

public:
	///	<summary>
	///		Maximum number of arguments in variable.
	///	</summary>
	static const int MaxArguments = 4;

	///	<summary>
	///		A sTime value indicating that this variable is not yet loaded
	///	</summary>
	static const size_t InvalidTimeIndex = (-2);

	///	<summary>
	///		A sTime value indicating that this variable has no time dimension
	///		(so only load once).
	///	</summary>
	static const size_t SingleTimeIndex = (-1);

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Variable(
		VariableRegistry * pvarreg
	) :
		m_pvarreg(pvarreg),
		m_strName(),
		m_pvarinfo(NULL),
		m_fOp(false),
		m_strUnits(),
		m_nSpecifiedDim(0),
		m_fNoTimeInNcFile(false),
		m_sTime(InvalidTimeIndex)
	{
		memset(m_iDim, 0, MaxArguments * sizeof(int));
	}

	///	<summary>
	///		Primitive Variable constructor.
	///	</summary>
	Variable(
		VariableRegistry * pvarreg,
		const std::string & strName,
		const VariableInfo * pvarinfo
	);

public:
	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Variable & var);

	///	<summary>
	///		Comparator operator.
	///	</summary>
	bool operator<(const Variable & var) {
		return (m_strName < var.m_strName);
	}

public:
	///	<summary>
	///		Parse the variable information from a string.  Return the index
	///		of the first character after the variable information.
	///	</summary>
	int ParseFromString(
		const std::string & strIn
	);

	///	<summary>
	///		Get a string representation of this variable.
	///	</summary>
	std::string ToString() const;

	///	<summary>
	///		Load a data block from the NcFileVector.
	///	</summary>
	void LoadGridData(
		RecapConfigObject * pobjConfig,
		size_t sTime
	);

	///	<summary>
	///		Unload the current data block.
	///	</summary>
	void UnloadGridData();

	///	<summary>
	///		Get the data associated with this variable.
	///	</summary>
	const DataArray1D<float> & GetData() const {
		return m_data;
	}

public:
	///	<summary>
	///		Name accessor.
	///	</summary>
	const std::string & Name() const {
		return m_strName;
	}

	///	<summary>
	///		Units accessor.
	///	</summary>
	const std::string & Units() const {
		return m_strUnits;
	}

public:
	///	<summary>
	///		Pointer to the VariableRegistry.
	///	</summary>
	VariableRegistry * m_pvarreg;

	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strName;

	///	<summary>
	///		Pointer to the associated VariableInfo structure in
	///		FileListObject, or NULL if this is a combination of variables.
	///	</summary>
	const VariableInfo * m_pvarinfo;

	///	<summary>
	///		Flag indicating this is an operator.
	///	</summary>
	bool m_fOp;

	///	<summary>
	///		Variable units.
	///	</summary>
	std::string m_strUnits;

	///	<summary>
	///		Number of dimensions specified.
	///	</summary>
	int m_nSpecifiedDim;

	///	<summary>
	///		Specified dimension values.
	///	</summary>
	int m_iDim[MaxArguments];

	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	VariableIndexVector m_varArg;

	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	std::vector<std::string> m_vecArg;

public:
	///	<summary>
	///		Flag indicating this Variable has no time index in NetCDF file.
	///	</summary>
	bool m_fNoTimeInNcFile;

	///	<summary>
	///		Time index associated with data loaded in this Variable.
	///	</summary>
	size_t m_sTime;

	///	<summary>
	///		Data associated with this Variable.
	///	</summary>
	DataArray1D<float> m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

