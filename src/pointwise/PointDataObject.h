///////////////////////////////////////////////////////////////////////////////
///
///	\file    PointDataObject.h
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

#ifndef _POINTDATAOBJECT_H_
#define _POINTDATAOBJECT_H_

#include "Announce.h"
#include "Object.h"

///////////////////////////////////////////////////////////////////////////////

class PointDataObject : public DistributedObject {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PointDataObject(const std::string & strName) :
		DistributedObject(strName)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new PointDataObject(strDuplicateName),
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
	) {
		if (strFunctionName == "output_csv") {
			if ((vecCommandLineType.size() != 1) ||
			    (vecCommandLineType[0] != ObjectType_String)
			) {
				return std::string("ERROR: Invalid parameters to function \"output_csv\"");
			}
			return std::string("");
		}
		return
			Object::Call(
				objreg,
				strFunctionName,
				vecCommandLine,
				vecCommandLineType,
				ppReturn);
	}

public:
	///	<summary>
	///		Set the number of fields used by the object.
	///	</summary>
	void SetFieldCount(
		size_t nIntFields = 0,
		size_t nFloatFields = 0,
		size_t nDoubleFields = 0
	) {
		m_nIntFields = nIntFields;
		m_nFloatFields = nFloatFields;
		m_nDoubleFields = nDoubleFields;

		for (size_t i = 0; i < m_vecDataInt.size(); i++) {
			m_vecDataInt[i].resize(m_nIntFields);
		}
		for (size_t i = 0; i < m_vecDataFloat.size(); i++) {
			m_vecDataFloat[i].resize(m_nFloatFields);
		}
		for (size_t i = 0; i < m_vecDataDouble.size(); i++) {
			m_vecDataDouble[i].resize(m_nDoubleFields);
		}
	}

	///	<summary>
	///		Change the number of rows in the data.
	///	</summary>
	void Resize(
		size_t nRows
	) {
		int nIntSize = m_vecDataInt.size();
		m_vecDataInt.resize(nRows);
		for (size_t i = nIntSize; i < m_vecDataInt.size(); i++) {
			m_vecDataInt[i].resize(m_nIntFields);
		}

		int nFloatSize = m_vecDataFloat.size();
		m_vecDataFloat.resize(nRows);
		for (size_t i = nFloatSize; i < m_vecDataFloat.size(); i++) {
			m_vecDataFloat[i].resize(m_nFloatFields);
		}

		int nDoubleSize = m_vecDataDouble.size();
		m_vecDataDouble.resize(nRows);
		for (size_t i = nDoubleSize; i < m_vecDataDouble.size(); i++) {
			m_vecDataDouble[i].resize(m_nDoubleFields);
		}
	}

	///	<summary>
	///		Get the number of rows in the data.
	///	</summary>
	size_t GetRows() const {
		return m_vecDataInt.size();
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
	int & DataInt(size_t i, size_t j) {
		return m_vecDataInt[i][j];
	}

	///	<summary>
	///		Data accessor to float data.
	///	</summary>
	float & DataFloat(size_t i, size_t j) {
		return m_vecDataFloat[i][j];
	}

	///	<summary>
	///		Data accessor to float data.
	///	</summary>
	double & DataDouble(size_t i, size_t j) {
		return m_vecDataDouble[i][j];
	}

	///	<summary>
	///		Concatenate multiple PointDataObjects.
	///	</summary>
	void Concatenate(
		const std::vector<PointDataObject *> vecpobjPointData
	) {
		if (vecpobjPointData.size() == 0) {
			return;
		}

		size_t sRows = 0;
		for (size_t s = 0; s < vecpobjPointData.size(); s++) {
			sRows += vecpobjPointData[s]->GetRows();
		}

		SetFieldCount(
			vecpobjPointData[0]->m_nIntFields,
			vecpobjPointData[0]->m_nFloatFields,
			vecpobjPointData[0]->m_nDoubleFields);

		Resize(sRows);

		sRows = 0;
		for (size_t s = 0; s < vecpobjPointData.size(); s++) {
			if (vecpobjPointData[s]->m_nIntFields != m_nIntFields) {
				_EXCEPTIONT("int field count mismatch");
			}
			if (vecpobjPointData[s]->m_nFloatFields != m_nFloatFields) {
				_EXCEPTIONT("float field count mismatch");
			}
			if (vecpobjPointData[s]->m_nDoubleFields != m_nDoubleFields) {
				_EXCEPTIONT("double field count mismatch");
			}

			memcpy(
				&(m_vecDataInt[sRows][0]),
				&(vecpobjPointData[s]->m_vecDataInt[0][0]),
				m_nIntFields * sizeof(int));

			memcpy(
				&(m_vecDataFloat[sRows][0]),
				&(vecpobjPointData[s]->m_vecDataFloat[0][0]),
				m_nFloatFields * sizeof(float));

			memcpy(
				&(m_vecDataDouble[sRows][0]),
				&(vecpobjPointData[s]->m_vecDataDouble[0][0]),
				m_nDoubleFields * sizeof(double));

			sRows += vecpobjPointData[s]->GetRows();
		}
	}

	///	<summary>
	///		Reduce data to the head rank.
	///	</summary>
	virtual void Reduce()
	{ }

protected:
	///	<summary>
	///		Number of integer data fields required.
	///	</summary>
	size_t m_nIntFields;

	///	<summary>
	///		Integer data fields.
	///	</summary>
	std::vector< std::vector<int> > m_vecDataInt;

	///	<summary>
	///		Number of float data fields required.
	///	</summary>
	size_t m_nFloatFields;

	///	<summary>
	///		Float data fields.
	///	</summary>
	std::vector< std::vector<float> > m_vecDataFloat;

	///	<summary>
	///		Number of double data fields required.
	///	</summary>
	size_t m_nDoubleFields;

	///	<summary>
	///		Double data fields.
	///	</summary>
	std::vector< std::vector<double> > m_vecDataDouble;
};

///////////////////////////////////////////////////////////////////////////////

#endif

