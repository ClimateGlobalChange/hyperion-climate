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
#include "MultiTypeDataArray2D.h"

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class PointDataObject :
	public Object,
	public MultiTypeDataArray2D
{

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PointDataObject(const std::string & strName) :
		Object(strName),
		MultiTypeDataArray2D()
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
	virtual void SetFieldCount(
		size_t nIntFields = 0,
		size_t nFloatFields = 0,
		size_t nDoubleFields = 0
	) {
		m_strIntFieldHeaders.resize(nIntFields);
		m_strFloatFieldHeaders.resize(nFloatFields);
		m_strDoubleFieldHeaders.resize(nDoubleFields);

		m_strFloatFieldUnits.resize(nFloatFields);
		m_strDoubleFieldUnits.resize(nDoubleFields);

		MultiTypeDataArray2D::SetFieldCount(
			nIntFields, nFloatFields, nDoubleFields);
	}

	///	<summary>
	///		Concatenate multiple PointDataObjects.
	///	</summary>
	virtual void Concatenate(
		const std::vector<PointDataObject *> & vecpobjPointData
	) {
		std::vector<MultiTypeDataArray2D *> vecpMultiTypeData;
		vecpMultiTypeData.resize(vecpobjPointData.size());
		for (size_t i = 0; i < vecpobjPointData.size(); i++) {
			vecpMultiTypeData[i] = vecpobjPointData[i];
		}
		MultiTypeDataArray2D::Concatenate(vecpMultiTypeData);
	}

public:
	///	<summary>
	///		Get the integer field header.
	///	</summary>
	const std::string & GetIntFieldHeader(size_t sFieldHeader) const {
		if (sFieldHeader > m_strIntFieldHeaders.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		return m_strIntFieldHeaders[sFieldHeader];
	}

	///	<summary>
	///		Set the integer field header.
	///	</summary>
	void SetIntFieldHeader(
		size_t sFieldHeader,
		const std::string & strIntFieldHeader
	) {
		if (sFieldHeader > m_strIntFieldHeaders.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		m_strIntFieldHeaders[sFieldHeader] = strIntFieldHeader;
	}

	///	<summary>
	///		Get the float field header.
	///	</summary>
	const std::string & GetFloatFieldHeader(size_t sFieldHeader) const {
		if (sFieldHeader > m_strFloatFieldHeaders.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		return m_strFloatFieldHeaders[sFieldHeader];
	}

	///	<summary>
	///		Set the float field header.
	///	</summary>
	void SetFloatFieldHeader(
		size_t sFieldHeader,
		const std::string & strFloatFieldHeader
	) {
		if (sFieldHeader > m_strFloatFieldHeaders.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		m_strFloatFieldHeaders[sFieldHeader] = strFloatFieldHeader;
	}

	///	<summary>
	///		Get the double field header.
	///	</summary>
	const std::string & GetDoubleFieldHeader(size_t sFieldHeader) const {
		if (sFieldHeader > m_strDoubleFieldHeaders.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		return m_strDoubleFieldHeaders[sFieldHeader];
	}

	///	<summary>
	///		Set the double field header.
	///	</summary>
	void SetDoubleFieldHeader(
		size_t sFieldHeader,
		const std::string & strDoubleFieldHeader
	) {
		if (sFieldHeader > m_strDoubleFieldHeaders.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		m_strDoubleFieldHeaders[sFieldHeader] = strDoubleFieldHeader;
	}

public:
	///	<summary>
	///		Get the float field header.
	///	</summary>
	const std::string & GetFloatFieldUnits(size_t sFieldUnits) const {
		if (sFieldUnits > m_strFloatFieldUnits.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		return m_strFloatFieldUnits[sFieldUnits];
	}

	///	<summary>
	///		Set the float field header.
	///	</summary>
	void SetFloatFieldUnits(
		size_t sFieldUnits,
		const std::string & strFloatFieldUnits
	) {
		if (sFieldUnits > m_strFloatFieldUnits.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		m_strFloatFieldUnits[sFieldUnits] = strFloatFieldUnits;
	}

	///	<summary>
	///		Get the double field header.
	///	</summary>
	const std::string & GetDoubleFieldUnits(size_t sFieldUnits) const {
		if (sFieldUnits > m_strDoubleFieldUnits.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		return m_strDoubleFieldUnits[sFieldUnits];
	}

	///	<summary>
	///		Set the double field header.
	///	</summary>
	void SetDoubleFieldUnits(
		size_t sFieldUnits,
		const std::string & strDoubleFieldUnits
	) {
		if (sFieldUnits > m_strDoubleFieldUnits.size()) {
			_EXCEPTIONT("Field header index out of range");
		}
		m_strDoubleFieldUnits[sFieldUnits] = strDoubleFieldUnits;
	}

public:
	///	<summary>
	///		Column headers for integer fields.
	///	</summary>
	std::vector<std::string> m_strIntFieldHeaders;

	///	<summary>
	///		Column headers for float fields.
	///	</summary>
	std::vector<std::string> m_strFloatFieldHeaders;

	///	<summary>
	///		Column headers for double fields.
	///	</summary>
	std::vector<std::string> m_strDoubleFieldHeaders;

	///	<summary>
	///		Column units for float fields.
	///	</summary>
	std::vector<std::string> m_strFloatFieldUnits;

	///	<summary>
	///		Column units for double fields.
	///	</summary>
	std::vector<std::string> m_strDoubleFieldUnits;
};

///////////////////////////////////////////////////////////////////////////////

#endif

