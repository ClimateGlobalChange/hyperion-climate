///////////////////////////////////////////////////////////////////////////////
///
///	\file    Regrid.cpp
///	\author  Paul Ullrich
///	\version March 14, 2017
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

#include "Regrid.h"

#include "TempestRegridObject.h"

///////////////////////////////////////////////////////////////////////////////

void HRegrid::RegisterGlobalFunctions(
	GlobalFunctionRegistry & funcreg
) {
	funcreg.Assign(
		std::string("tempestregrid"),
		new TempestRegridObjectConstructor(std::string("tempestregrid")));
}

///////////////////////////////////////////////////////////////////////////////

