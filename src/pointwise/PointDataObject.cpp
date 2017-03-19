///////////////////////////////////////////////////////////////////////////////
///
///	\file    PointDataObject.cpp
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

#include "PointDataObject.h"

#include <fstream>

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

std::string PointDataObject::OutputCSV(
	const std::string & strCSVOutputFilename
) {
	// Gather data to root node
	if (IsDistributed()) {
		Gather();
	}

#if defined(HYPERION_MPIOMP)
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

	// Output the headers
	MultiTypeDataHeaders::OutputCSV(ofOutput);
	MultiTypeDataArray2D::OutputCSV(ofOutput);

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

