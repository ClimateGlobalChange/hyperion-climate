///////////////////////////////////////////////////////////////////////////////
///
///	\file    MultiTypeDataHeaders.cpp
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

#include "MultiTypeDataHeaders.h"

///////////////////////////////////////////////////////////////////////////////

void MultiTypeDataHeaders::OutputCSV(
	std::ostream & osOutput
) const {
	// Column headers
	bool fInitialComma = false;

	for (int j = 0; j < m_strIntFieldHeaders.size(); j++) {
		if (!fInitialComma) {
			fInitialComma = true;
		} else {
			osOutput << ",";
		}
		osOutput << m_strIntFieldHeaders[j];
	}
	for (int j = 0; j < m_strFloatFieldHeaders.size(); j++) {
		if (!fInitialComma) {
			fInitialComma = true;
		} else {
			osOutput << ",";
		}
		osOutput << m_strFloatFieldHeaders[j]
	   		<< "[" << m_strFloatFieldUnits[j] << "]";
	}
	for (int j = 0; j < m_strDoubleFieldHeaders.size(); j++) {
		if (!fInitialComma) {
			fInitialComma = true;
		} else {
			osOutput << ",";
		}
		osOutput << m_strDoubleFieldHeaders[j]
	   		<< "[" << m_strDoubleFieldUnits[j] << "]";
	}
	osOutput << std::endl;
}

///////////////////////////////////////////////////////////////////////////////

