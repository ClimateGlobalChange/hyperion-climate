///////////////////////////////////////////////////////////////////////////////
///
///	\file    MultiTypeDataArray2D.h
///	\author  Paul Ullrich
///	\version March 17, 2017
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

#ifndef _MULTITYPEDATAARRAY2D_H_
#define _MULTITYPEDATAARRAY2D_H_

#include "Announce.h"
#include "Object.h"
#include "DataArray2D.h"

#include <string>
#include <vector>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

class MultiTypeDataArray2D {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	MultiTypeDataArray2D() :
		m_nRows(0),
		m_nIntFields(0),
		m_nFloatFields(0),
		m_nDoubleFields(0)
	{ }

public:
	///	<summary>
	///		Set the number of fields used by the object.
	///	</summary>
	virtual void SetFieldCount(
		size_t nIntFields = 0,
		size_t nFloatFields = 0,
		size_t nDoubleFields = 0
	);

	///	<summary>
	///		Change the number of rows in the data.
	///	</summary>
	virtual void Resize(
		size_t nRows
	);

	///	<summary>
	///		Concatenate multiple PointDataObjects.
	///	</summary>
	virtual void Concatenate(
		const std::vector<MultiTypeDataArray2D *> & vecpobjPointData
	);

	///	<summary>
	///		Gather data to the head rank.
	///	</summary>
	virtual void Gather();

	///	<summary>
	///		Output as a comma-separated list to the given stream.
	///	</summary>
	virtual void OutputCSV(
		std::ostream & osOutput
	) const;

public:
	///	<summary>
	///		Get the number of rows in the data.
	///	</summary>
	size_t GetRows() const {
		return m_dDataInt.GetRows();
	}

	///	<summary>
	///		Get the number of integer fields in the data.
	///	</summary>
	size_t GetIntFieldCount() const {
		return m_nIntFields;
	}

	///	<summary>
	///		Get the number of float fields in the data.
	///	</summary>
	size_t GetFloatFieldCount() const {
		return m_nFloatFields;
	}

	///	<summary>
	///		Get the number of double fields in the data.
	///	</summary>
	size_t GetDoubleFieldCount() const {
		return m_nDoubleFields;
	}

	///	<summary>
	///		Data accessor to integer data.
	///	</summary>
	inline int & DataInt(size_t i, size_t j) {
		return m_dDataInt(i,j);
	}

	///	<summary>
	///		Data accessor to integer data.
	///	</summary>
	inline const int & DataInt(size_t i, size_t j) const {
		return m_dDataInt(i,j);
	}

	///	<summary>
	///		Data accessor to float data.
	///	</summary>
	inline float & DataFloat(size_t i, size_t j) {
		return m_dDataFloat(i,j);
	}

	///	<summary>
	///		Data accessor to float data.
	///	</summary>
	inline const float & DataFloat(size_t i, size_t j) const {
		return m_dDataFloat(i,j);
	}

	///	<summary>
	///		Data accessor to float data.
	///	</summary>
	inline double & DataDouble(size_t i, size_t j) {
		return m_dDataDouble(i,j);
	}

	///	<summary>
	///		Data accessor to float data.
	///	</summary>
	inline const double & DataDouble(size_t i, size_t j) const {
		return m_dDataDouble(i,j);
	}

protected:
	///	<summary>
	///		Number of rows in data.
	///	</summary>
	size_t m_nRows;

	///	<summary>
	///		Number of integer data fields required.
	///	</summary>
	size_t m_nIntFields;

	///	<summary>
	///		Integer data fields.
	///	</summary>
	DataArray2D<int> m_dDataInt;

	///	<summary>
	///		Number of float data fields required.
	///	</summary>
	size_t m_nFloatFields;

	///	<summary>
	///		Float data fields.
	///	</summary>
	DataArray2D<float> m_dDataFloat;

	///	<summary>
	///		Number of double data fields required.
	///	</summary>
	size_t m_nDoubleFields;

	///	<summary>
	///		Double data fields.
	///	</summary>
	DataArray2D<double> m_dDataDouble;

};

///////////////////////////////////////////////////////////////////////////////

#endif

