///////////////////////////////////////////////////////////////////////////////
///
///	\file    Object.h
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

#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "Exception.h"

#include <string>
#include <set>

///////////////////////////////////////////////////////////////////////////////

class Object;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for registering Objects.
///	</summary>
class ObjectRegistry {

private:
	///	<summary>
	///		Set of Objects.
	///	</summary>
	std::set<Object> m_setObject;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for representing objects.
///	</summary>
class Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Object(const std::string & strName) :
		m_strName(strName),
		m_fOwnsData(false),
		m_pData(NULL)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~Object() {
		if (m_fOwnsData) {
			free(m_pData);
		}
	}

public:
	///	<summary>
	///		Inequality comparator.
	///	</summary>
	bool operator< (const Object & obj) {
		return (m_strName < obj.m_strName);
	}

	///	<summary>
	///		Equality comparator.
	///	</summary>
	bool operator== (const Object & obj) {
		return (m_strName == obj.m_strName);
	}

protected:
	///	<summary>
	///		Name of the object.
	///	</summary>
	std::string m_strName;

	///	<summary>
	///		List of attributes.
	///	</summary>
	std::set<Object> m_setAttributes;

	///	<summary>
	///		Flag indicating that this object owns its data.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		Data associatd with this object.
	///	</summary>
	void * m_pData;
};

///////////////////////////////////////////////////////////////////////////////

#endif

