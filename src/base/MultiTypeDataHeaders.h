///////////////////////////////////////////////////////////////////////////////
///
///	\file    MultiTypeDataHeaders.h
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

#ifndef _MULTITYPEDATAHEADERS_H_
#define _MULTITYPEDATAHEADERS_H_

#include "Announce.h"
#include "Exception.h"

#include <string>
#include <vector>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

class MultiTypeDataHeaders {

protected:
	///	<summary>
	///		Constructor.
	///	</summary>
	MultiTypeDataHeaders()
	{ }

public:
	///	<summary>
	///		Copy over headers from another MultiTypeDataHeaders.
	///	</summary>
	void SetHeaders(const MultiTypeDataHeaders & headers) {
		m_strIntFieldHeaders = headers.m_strIntFieldHeaders;
		m_strFloatFieldHeaders = headers.m_strFloatFieldHeaders;
		m_strDoubleFieldHeaders = headers.m_strDoubleFieldHeaders;
		m_strFloatFieldUnits = headers.m_strFloatFieldUnits;
		m_strDoubleFieldUnits = headers.m_strDoubleFieldUnits;
	}

protected:
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
		m_strFloatFieldUnits.resize(nFloatFields);
		m_strDoubleFieldHeaders.resize(nDoubleFields);
		m_strDoubleFieldUnits.resize(nDoubleFields);
	}

	///	<summary>
	///		Output as a comma-separated list to the given stream.
	///	</summary>
	virtual void OutputCSV(
		std::ostream & osOutput
	) const;

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

protected:
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

