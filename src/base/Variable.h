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

#include "AccessMode.h"
#include "DataArray1D.h"

#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class RecapConfigObject;

class VariableInfo;

class Variable;

typedef std::vector<Variable> VariableVector;

typedef int VariableIndex;

typedef std::vector<VariableIndex> VariableIndexVector;

///////////////////////////////////////////////////////////////////////////////

class VariableRegistry {

public:
	///	<summary>
	///		Constructor for VariableRegistry.
	///	</summary>
	VariableRegistry(
		RecapConfigObject * pobjConfig
	);

public:
	///	<summary>
	///		Get the pointer to the RecapConfigObject.
	///	</summary>
	RecapConfigObject * GetRecapConfigObject() {
		return m_pobjConfig;
	}

	///	<summary>
	///		Find or register a variable.  Automatically recursively identify
	///		arguments if this is an operator combination, unless fStandalone
	///		is set to true.
	///	</summary>
	std::string FindOrRegister(
		const std::string & strVariableName,
		Variable ** ppVariable
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

	///	<summary>
	///		Update time indices for all Variables.
	///	</summary>
	void UpdateTimeIndices(
		const std::map<size_t, size_t> & mapTimeIxToNewTimeIx
	);

private:
	///	<summary>
	///		Pointer to parent RecapConfigObject.
	///	</summary>
	RecapConfigObject * m_pobjConfig;

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

friend class VariableRegistry;

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
	///		Primitive Variable constructor.
	///	</summary>
	Variable(
		VariableRegistry * pvarreg,
		const std::string & strName,
		const VariableInfo * pvarinfo
	);

private:
	///	<summary>
	///		Initialize the Variable.
	///	</summary>
	std::string Initialize();

	///	<summary>
	///		Initialize auxiliary information associated with this Variable.
	///	</summary>
	std::string InitializeAuxiliary(
		const std::vector<Variable *> vecVarArguments
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
	///		Get a string representation of this variable.
	///	</summary>
	std::string ToString() const;

	///	<summary>
	///		Load a data block from disk.
	///	</summary>
	std::string LoadGridData(
		size_t sTime
	);

	///	<summary>
	///		Allocate an empty data block for the specified variable.
	///	</summary>
	void AllocateGridData(
		size_t sTime
	);

	///	<summary>
	///		Write a data block to disk.
	///	</summary>
	std::string WriteGridData(
		size_t sTime
	);

	///	<summary>
	///		Unload the current data block.
	///	</summary>
	void UnloadGridData();

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

	///	<summary>
	///		Dimension names accessor.
	///	</summary>
	const std::vector<std::string> & DimNames() const {
		return m_vecDimNames;
	}

	///	<summary>
	///		Check if this is an operator.
	///	</summary>
	bool IsOp() const {
		return m_fOp;
	}

	///	<summary>
	///		Check if this is a reduction operator.
	///	</summary>
	bool IsReductionOp() const {
		return m_fReductionOp;
	}

	///	<summary>
	///		Get the data associated with this variable.
	///	</summary>
	const DataArray1D<float> & GetData() const {
		return m_data;
	}

	///	<summary>
	///		Get the data associated with this variable.
	///	</summary>
	DataArray1D<float> & GetData() {
		return m_data;
	}

protected:
	///	<summary>
	///		Pointer to the VariableRegistry.
	///	</summary>
	VariableRegistry * m_pvarreg;

	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strName;

	///	<summary>
	///		Variable units.
	///	</summary>
	std::string m_strUnits;

	///	<summary>
	///		Auxiliary dimensions associated with this variable.
	///	</summary>
	std::vector<std::string> m_vecDimNames;

	///	<summary>
	///		Auxiliary indices associated with data loaded in this Variable.
	///	</summary>
	std::vector<long> m_vecAuxIndices;

	///	<summary>
	///		Index of the record dimension.
	///	</summary>
	int m_iTimeDimIx;

protected:
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
	///		Flag indicating this is a reduction operator.
	///	</summary>
	bool m_fReductionOp;

	///	<summary>
	///		Number of dimensions specified.
	///	</summary>
	int m_nSpecifiedDim;

	///	<summary>
	///		Specified dimension values.
	///	</summary>
	int m_iDim[MaxArguments];

	///	<summary>
	///		Name of the operator.
	///	</summary>
	std::string m_strOpName;
/*
	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	VariableIndexVector m_varArg;
*/
	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	std::vector<std::string> m_vecArg;

protected:
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

