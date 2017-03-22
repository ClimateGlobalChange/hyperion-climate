///////////////////////////////////////////////////////////////////////////////
///
///	\file    Regrid.h
///	\author  Paul Ullrich
///	\version March 22, 2017
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

#ifndef _REGRID_H_
#define _REGRID_H_

#include "GlobalFunction.h"

namespace HRegrid {

///////////////////////////////////////////////////////////////////////////////

void RegisterGlobalFunctions(
	GlobalFunctionRegistry & funcreg
);

///////////////////////////////////////////////////////////////////////////////

};

#endif

