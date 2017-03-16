///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchNodes.h
///	\author  Paul Ullrich
///	\version March 16, 2017
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

#ifndef _STITCHNODES_H_
#define _STITCHNODES_H_

#include "GlobalFunction.h"

namespace HPointwise {

///////////////////////////////////////////////////////////////////////////////

class StitchNodesFunction : public GlobalFunction {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	StitchNodesFunction(const std::string & strName) :
		GlobalFunction(strName)
	{ }

public:
	///	<summary>
	///		Call a member function of this GlobalFunction.
	///	</summary>
	virtual std::string Call(
		const ObjectRegistry & objreg,
		const std::vector<std::string> & vecCommandLine,
		const std::vector<ObjectType> & vecCommandLineType,
		Object ** ppReturn
	);
};

///////////////////////////////////////////////////////////////////////////////

};

#endif

