///////////////////////////////////////////////////////////////////////////////
///
///	\file    PathDataObject.h
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

#ifndef _PATHDATAOBJECT_H_
#define _PATHDATAOBJECT_H_

#include "Announce.h"
#include "Object.h"
#include "MultiTypeDataArray2D.h"

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class PathDataObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PathDataObject(const std::string & strName) :
		Object(strName)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new PathDataObject(strDuplicateName),
			objreg);
	}

	///	<summary>
	///		Call a member function of this Object.
	///	</summary>
	virtual std::string Call(
		const ObjectRegistry & objreg,
		const std::string & strFunctionName,
		const std::vector<std::string> & vecCommandLine,
		const std::vector<ObjectType> & vecCommandLineType,
		Object ** ppReturn
	) {
		if (strFunctionName == "output_csv") {
			if ((vecCommandLineType.size() != 1) ||
			    (vecCommandLineType[0] != ObjectType_String)
			) {
				return std::string("ERROR: Invalid parameters to function \"output_csv\"");
			}
			return std::string("");
		}
		return
			Object::Call(
				objreg,
				strFunctionName,
				vecCommandLine,
				vecCommandLineType,
				ppReturn);
	}

	///	<summary>
	///		Get the vector of path data information.
	///	</summary>
	std::vector< MultiTypeDataArray2D > & GetData() {
		return m_vecPathData;
	}

	///	<summary>
	///		Get the vector of path data information.
	///	</summary>
	const std::vector< MultiTypeDataArray2D > & GetData() const {
		return m_vecPathData;
	}

protected:
	///	<summary>
	///		A vector of path data information.
	///	</summary>
	std::vector< MultiTypeDataArray2D > m_vecPathData;
};

///////////////////////////////////////////////////////////////////////////////

#endif

