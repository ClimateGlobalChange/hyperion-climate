///////////////////////////////////////////////////////////////////////////////
///
///	\file    Object.cpp
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

#include "Object.h"

///////////////////////////////////////////////////////////////////////////////
// ObjectRegistry
///////////////////////////////////////////////////////////////////////////////

ObjectRegistry::~ObjectRegistry() {
	ObjectMap::iterator iter = m_mapObjects.begin();
	for (; iter != m_mapObjects.end(); iter++) {
		delete(iter->second);
	}
}

///////////////////////////////////////////////////////////////////////////////

Object * ObjectRegistry::GetObject(
	const std::string & strName
) const {

	ObjectMap::const_iterator iter = m_mapObjects.find(strName);
	if (iter == m_mapObjects.end()) {
		return NULL;
	} else {
		return (iter->second);
	}
}

///////////////////////////////////////////////////////////////////////////////

bool ObjectRegistry::Insert(
	const std::string & strName,
	Object * pObject
) {
	//printf("INSERT %s\n", strName.c_str());

	// Check if this Object already exists
	ObjectMap::const_iterator iter = m_mapObjects.find(strName);
	if (iter != m_mapObjects.end()) {
		_EXCEPTION1("Duplicate Object \"%s\"", strName.c_str());
	}

	// Add the Object to its parent
	std::string strParent;
	for (int i = strName.length()-1; i >= 0; i--) {
		if (strName[i] == '.') {
			strParent = strName.substr(0,i);
			break;
		}
	}

	// Insert the Object into the parent's array
	if (strParent != "") {
		ObjectMap::const_iterator iterParent = m_mapObjects.find(strParent);
		if (iterParent == m_mapObjects.end()) {
			return false;
		}

		//printf("PARENT %s\n", strParent.c_str());
		iterParent->second->m_setChildren.insert(pObject);
	}

	// Insert the Object into the ObjectRegistry
	m_mapObjects.insert(
		ObjectMap::value_type(strName, pObject));

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Object
///////////////////////////////////////////////////////////////////////////////

Object * Object::_Duplicate(
	Object * pobjDuplicate,
	ObjectRegistry & objreg
) const {
	bool fSuccess = objreg.Insert(pobjDuplicate->m_strName, pobjDuplicate);
	if (!fSuccess) {
		_EXCEPTIONT("Failed to insert Object into registry");
	}

	const int nOriginalNameLength = m_strName.length();

	const int nNewNameLength = pobjDuplicate->m_strName.length();

	ObjectChildrenSet::const_iterator iter = m_setChildren.begin();
	for (; iter != m_setChildren.end(); iter++) {

		printf("ORIGINAL CHILD %s\n", (*iter)->m_strName.c_str());
		std::string strNewChildName =
			pobjDuplicate->m_strName
			+ (*iter)->m_strName.substr(
				nOriginalNameLength, std::string::npos);

		if (strNewChildName[nNewNameLength] != '.') {
			_EXCEPTIONT("Logic error: Invalid child");
		}

		(*iter)->Duplicate(strNewChildName, objreg);
	}

	return (pobjDuplicate);
}

///////////////////////////////////////////////////////////////////////////////

