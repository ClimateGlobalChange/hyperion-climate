///////////////////////////////////////////////////////////////////////////////
///
///	\file    recap.cpp
///	\author  Paul Ullrich
///	\version March 8, 2017
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

#include "netcdfcpp.h"

#include "Announce.h"
#include "Object.h"

#include <iostream>
#include <fstream>
#include <string>

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Output usage information
	if (argc < 2) {
		Announce("recap version 0.1, March 8th, 2017");
		Announce("Usage: %s <input file> [assignment list]", argv[0]);
		return (0);
	}

	// Load input file
	std::ifstream ifScript(argv[1]);
	if (!ifScript.is_open()) {
		Announce("ERROR: Unable to open input file \"%s\"", argv[1]);
		return (-1);
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}

///////////////////////////////////////////////////////////////////////////////

