///////////////////////////////////////////////////////////////////////////////
///
///	\file    PathDataObject.cpp
///	\author  Paul Ullrich
///	\version March 18, 2017
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

#include "PathDataObject.h"

#include <fstream>

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

std::string PathDataObject::OutputCSV(
	const std::string & strCSVOutputFilename
) {
#if defined(HYPERION_MPIOMP)
	// Only output on root thread
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	if (nRank != 0) {
		return std::string("");
	}
#endif

	// Open output file
	std::ofstream ofOutput(strCSVOutputFilename.c_str());
	if (!ofOutput.is_open()) {
		return std::string("Unable to open output file \"") + strCSVOutputFilename + "\"";
	}

	if (m_vecPathData.size() == 0) {
		return std::string("");
	}

	MultiTypeDataHeaders::OutputCSV(ofOutput);

	for (size_t i = 0; i < m_vecPathData.size(); i++) {
		ofOutput << "Path " << i << std::endl;
		m_vecPathData[i].OutputCSV(ofOutput);
	}
/*
	// Column headers
	bool fInitialComma = false;
	for (int j = 0; j < GetIntFieldCount(); j++) {
		if (!fInitialComma) {
			fInitialComma = true;
		} else {
			ofOutput << ",";
		}
		ofOutput << m_strIntFieldHeaders[j];
	}
	for (int j = 0; j < GetFloatFieldCount(); j++) {
		if (!fInitialComma) {
			fInitialComma = true;
		} else {
			ofOutput << ",";
		}
		ofOutput << m_strFloatFieldHeaders[j]
	   		<< "[" << m_strFloatFieldUnits[j] << "]";
	}
	for (int j = 0; j < GetDoubleFieldCount(); j++) {
		if (!fInitialComma) {
			fInitialComma = true;
		} else {
			ofOutput << ",";
		}
		ofOutput << m_strDoubleFieldHeaders[j]
	   		<< "[" << m_strDoubleFieldUnits[j] << "]";
	}
	ofOutput << std::endl;
	
	// Loop through all rows
	for (int i = 0; i < GetRows(); i++) {
		fInitialComma = false;

		// Output integer data
		for (int j = 0; j < GetIntFieldCount(); j++) {
			if (!fInitialComma) {
				fInitialComma = true;
			} else {
				ofOutput << ",";
			}
			ofOutput << DataInt(i,j);
		}
		for (int j = 0; j < GetFloatFieldCount(); j++) {
			if (!fInitialComma) {
				fInitialComma = true;
			} else {
				ofOutput << ",";
			}
			ofOutput << DataFloat(i,j);
		}
		for (int j = 0; j < GetDoubleFieldCount(); j++) {
			if (!fInitialComma) {
				fInitialComma = true;
			} else {
				ofOutput << ",";
			}
			ofOutput << DataDouble(i,j);
		}
		ofOutput << std::endl;
	}
*/
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

