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

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A map from the unique variable identifier to a pointer for that
///		variable.
///	</summary>
typedef std::map<std::string, Variable *> VariableMap;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An object holding dimension indices for a given Variable.
///	</summary>
typedef std::vector<long> VariableDimIndex;

///	<summary>
///		An object holding auxiliary indices for a given Variable.
///	</summary>
typedef VariableDimIndex VariableAuxIndex;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A registry holding information and data on active Variables.
///	</summary>
class VariableRegistry {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableRegistry(
		RecapConfigObject * pobjConfig
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	~VariableRegistry();

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
	VariableMap m_mapVariables;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class describing an iterator over the auxiliary indices of a
///		Variable.
///	</summary>
class VariableAuxIndexIterator {

friend class Variable;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableAuxIndexIterator() :
		m_fEnd(false)
	{ }

protected:
	///	<summary>
	///		Initializer, only accessible from Variable.
	///	</summary>
	void Initialize(
		const VariableAuxIndex & vecSize,
		bool fEnd
	) {
		m_fEnd = fEnd;
		m_vecSize = vecSize;
		m_vecValue.resize(m_vecSize.size());
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			if (m_vecSize[d] <= 0) {
				_EXCEPTIONT("Invalid auxiliary index size entry");
			}
		}
		if ((fEnd) && (vecSize.size() > 0)) {
			m_vecValue[0] = m_vecSize[0];
		}
	}

	///	<summary>
	///		Remove the dimension with specified index.
	///	</summary>
	void RemoveDim(
		size_t dim
	) {
		if (dim >= m_vecSize.size()) {
			_EXCEPTIONT("Logic errror");
		}
		m_vecSize.erase(m_vecSize.begin() + dim);
		m_vecValue.erase(m_vecValue.begin() + dim);
	}

public:
	///	<summary>
	///		Prefix incrementor.
	///	</summary>
	VariableAuxIndexIterator & operator++() {
		if (m_vecValue.size() == 0) {
			if (m_fEnd) {
				_EXCEPTIONT("Iterator exceeded bounds");
			} else {
				m_fEnd = true;
			}
			return (*this);
		}
		for (size_t d = 0; d < m_vecValue.size(); d++) {
			size_t dx = m_vecValue.size()-d-1;
			if (m_vecValue[dx] >= m_vecSize[dx]) {
				_EXCEPTIONT("Iterator exceeded bounds");
			} else if (m_vecValue[dx] == m_vecSize[dx]-1) {
				m_vecValue[dx] = 0;
			} else {
				m_vecValue[dx]++;
				return (*this);
			}
		}
		m_vecValue[0] = m_vecSize[0];
		m_fEnd = true;
		return (*this);
	}

	///	<summary>
	///		Postfix incrementor.
	///	</summary>
	VariableAuxIndexIterator operator++(int) {
		this->operator++();
		return (*this);
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator==(const VariableAuxIndexIterator & iter) const {
		if (m_vecSize.size() != iter.m_vecSize.size()) {
			_EXCEPTIONT("Invalid comparison");
		}
		if (m_vecSize.size() == 0) {
			if (m_fEnd == iter.m_fEnd) {
				return true;
			} else {
				return false;
			}
		}
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			if (m_vecValue[d] != iter.m_vecValue[d]) {
				return false;
			}
		}
		return true;
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator!=(const VariableAuxIndexIterator & iter) const {
		return !((*this)==(iter));
	}

	///	<summary>
	///		Cast to VariableAuxIndex.
	///	</summary>
	operator VariableAuxIndex() const {
		return m_vecValue;
	}

	///	<summary>
	///		Get the value of the iterator.
	///	</summary>
	const VariableAuxIndex & Value() const {
		return m_vecValue;
	}

	///	<summary>
	///		Convert to a std::string.
	///	</summary>
	std::string ToString() const {
		std::string strOut;
		if (m_vecSize.size() == 0) {
			if (m_fEnd) {
				strOut = std::string("(end)");
			} else {
				strOut = std::string("(begin)");
			}
			return strOut;
		}
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			strOut +=
				std::string("(")
				+ std::to_string(m_vecValue[d])
				+ std::string("/")
				+ std::to_string(m_vecSize[d])
				+ std::string(")");
		}
		return strOut;
	}

protected:
	///	<summary>
	///		The sizes of the auxiliary indices.
	///	</summary>
	VariableAuxIndex m_vecSize;

	///	<summary>
	///		The values of the auxiliary indices.
	///	</summary>
	VariableAuxIndex m_vecValue;

	///	<summary>
	///		At the end.
	///	</summary>
	bool m_fEnd;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a 2D slice of data.
///	</summary>
class Variable {

friend class VariableRegistry;

public:
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
	///		Load a data block from disk using the specified time index.
	///	</summary>
	std::string LoadGridData(
		size_t sTime
	);

	///	<summary>
	///		Load a data block from disk.
	///	</summary>
	std::string LoadGridData(
		const VariableAuxIndex & ixAux
	);

	///	<summary>
	///		Allocate an empty data block for the specified variable.
	///	</summary>
	void AllocateGridData();

	///	<summary>
	///		Load a data block from disk.
	///	</summary>
	std::string WriteGridData(
		const VariableAuxIndex & ixAux
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
	const std::vector<std::string> & AuxDimNames() const {
		return m_vecAuxDimNames;
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
	///		Get the VariableAuxIndexIterator corresponding to the first
	///		auxiliary index.
	///	</summary>
	const VariableAuxIndexIterator & GetAuxIndexBegin() const {
		return m_iterAuxBegin;
	}

	///	<summary>
	///		Get the VariableAuxIndexIterator corresponding to the end
	///		auxiliary index.
	///	</summary>
	const VariableAuxIndexIterator & GetAuxIndexEnd() const {
		return m_iterAuxEnd;
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
	///		Pointer to the associated VariableInfo structure in
	///		FileListObject, or NULL if this is a combination of variables.
	///	</summary>
	const VariableInfo * m_pvarinfo;

protected:
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
	std::vector<std::string> m_vecAuxDimNames;

	///	<summary>
	///		Index of the record dimension.
	///	</summary>
	int m_iTimeDimIx;

protected:
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
	///		Name of the operator.
	///	</summary>
	std::string m_strOpName;

	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	std::vector<std::string> m_vecArg;

protected:
	///	<summary>
	///		Auxiliary index iterator "begin".
	///	</summary>
	VariableAuxIndexIterator m_iterAuxBegin;

	///	<summary>
	///		Auxiliary index iterator "end".
	///	</summary>
	VariableAuxIndexIterator m_iterAuxEnd;

protected:
	///	<summary>
	///		Information on the slice of data currently stored in this Variable.
	///	</summary>
	//DataSliceTag m_tag;

	///	<summary>
	///		Flag indicating that this single time index Variable contains
	///		active data.
	///	</summary>
	bool m_fSingleTimeIndexAndHasData;

	///	<summary>
	///		Auxiliary indices associated with data loaded in this Variable.
	///	</summary>
	VariableAuxIndex m_vecAuxIndices;

	///	<summary>
	///		Data associated with this Variable.
	///	</summary>
	DataArray1D<float> m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

