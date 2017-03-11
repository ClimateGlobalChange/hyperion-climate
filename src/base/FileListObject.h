///////////////////////////////////////////////////////////////////////////////
///
///	\file    FileListObject.h
///	\author  Paul Ullrich
///	\version March 10, 2017
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

#ifndef _FILELISTOBJECT_H_
#define _FILELISTOBJECT_H_

#include "Announce.h"
#include "Object.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing the grid, including coordinates of
///		each data point and graph connectivity of elements.
///	</summary>
class FileListObject : public ListObject {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FileListObject(
		const std::string & strName
	) :
		ListObject(strName)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new FileListObject(strDuplicateName),
			objreg);
	}

public:
	///	<summary>
	///		Populate from a search string.
	///	</summary>
	///	<returns>
	///		false if the directory specified by strSearchString does
	///		not exist.  true otherwise.
	///	</returns>
	bool PopulateFromSearchString(
		const std::string & strSearchString,
		ObjectRegistry & objreg
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

