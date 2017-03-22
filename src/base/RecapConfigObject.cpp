///////////////////////////////////////////////////////////////////////////////
///
///	\file    RecapConfigObject.cpp
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

#include "RecapConfigObject.h"

///////////////////////////////////////////////////////////////////////////////

std::string RecapConfigObjectConstructor::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	RecapConfigObject * pobjRecapConfig = new RecapConfigObject("");
	if (pobjRecapConfig == NULL) {
		_EXCEPTIONT("Unable to initialize RecapConfigObject");
	}

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjRecapConfig;
	} else {
		delete pobjRecapConfig;
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

