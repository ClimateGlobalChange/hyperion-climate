///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShpLoader.h
///	\author  Paul Ullrich
///	\version March 30, 2017
///
///	<remarks>
///		Copyright 2000- Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _SHPLOADER_H_
#define _SHPLOADER_H_

#include <string>

///////////////////////////////////////////////////////////////////////////////

class Mesh;

///////////////////////////////////////////////////////////////////////////////

std::string LoadShpFile(
	const std::string & strShpFilename,
	Mesh & mesh,
	bool fAppend = false
);

///////////////////////////////////////////////////////////////////////////////

#endif

