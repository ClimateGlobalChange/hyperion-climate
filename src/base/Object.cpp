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

#include <cstdlib>

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

bool ObjectRegistry::Create(
	const ObjectType & objtype,
	const std::string & strName,
	const std::string & strValue
) {
	// String type on RHS
	if (objtype == ObjectType_String) {
		return Assign(
			strName,
			new StringObject(
				strName,
				strValue));

	// Integer type on RHS
	} else if (objtype == ObjectType_Integer) {
		return Assign(
			strName,
			new IntegerObject(
				strName,
				atoi(strValue.c_str())));

	// Floating point type on RHS
	} else if (objtype == ObjectType_FloatingPoint) {
		return Assign(
			strName,
			new FloatingPointObject(
				strName,
				atof(strValue.c_str())));

	// Object type on RHS
	} else if (objtype == ObjectType_Token) {
		Object * pobj = GetObject(strValue);
		if (pobj != NULL) {
			Object * pobjDuplicate = pobj->Duplicate(strName, *this);
			if (pobjDuplicate != NULL) {
				return true;
			} else {
				_EXCEPTIONT("Unknown error");
			}
		} else {
			return false;
		}
	}

	// Invalid type
	return false;
}

///////////////////////////////////////////////////////////////////////////////

void ObjectRegistry::Remove(
	const std::string & strName
) {
	// Check if this Object already exists
	ObjectMap::iterator iter = m_mapObjects.find(strName);
	if (iter == m_mapObjects.end()) {
		_EXCEPTION1("Object \"%s\" not found in registry", strName.c_str());
	}

	printf("REMOVE %s\n", strName.c_str());

	// Remove all children of this Object
	ObjectChildrenSet::const_iterator iterChildren
		= iter->second->m_setChildren.begin();
	for (; iterChildren != iter->second->m_setChildren.end(); iterChildren++) {
		Remove((*iterChildren)->m_strName);
	}

	// Remove this Object
	delete (iter->second);
	m_mapObjects.erase(iter);
}

///////////////////////////////////////////////////////////////////////////////

bool ObjectRegistry::Assign(
	const std::string & strName,
	Object * pObject
) {
	//printf("ASSIGN %s\n", strName.c_str());

	// Check if this Object already exists
	ObjectMap::const_iterator iter = m_mapObjects.find(strName);
	if (iter != m_mapObjects.end()) {
		Remove(strName);
	}

	// Add the Object to its parent
	std::string strParent;
	for (int i = strName.length()-1; i >= 0; i--) {
		if (strName[i] == '.') {
			strParent = strName.substr(0,i);
			break;
		}
	}

	// Assign the Object into the parent's array
	if (strParent != "") {
		ObjectMap::const_iterator iterParent = m_mapObjects.find(strParent);
		if (iterParent == m_mapObjects.end()) {
			return false;
		}

		//printf("PARENT %s\n", strParent.c_str());
		iterParent->second->m_setChildren.insert(pObject);
	}

	// Assign the Object into the ObjectRegistry
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
	bool fSuccess = objreg.Assign(pobjDuplicate->m_strName, pobjDuplicate);
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

