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
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class Object;

typedef std::set<Object *> ObjectChildrenSet;

typedef std::map<std::string, Object *> ObjectMap;

typedef int ObjectIndex;

///////////////////////////////////////////////////////////////////////////////

typedef int ObjectType;

static const ObjectType ObjectType_Op = 0;

static const ObjectType ObjectType_Token = 1;

static const ObjectType ObjectType_String = 2;

static const ObjectType ObjectType_Integer = 3;

static const ObjectType ObjectType_FloatingPoint = 4;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for registering Objects.
///	</summary>
class ObjectRegistry {

public:
	///	<summary>
	///		Destructor.
	///	</summary>
	~ObjectRegistry();

public:
	///	<summary>
	///		Get the Object with the specified name.
	///	</summary>
	Object * GetObject(
		const std::string & strName
	) const;

	///	<summary>
	///		Create an Object of the given type.
	///	</summary>
	bool Create(
		const ObjectType & objtype,
		const std::string & strName,
		const std::string & strValue
	);

	///	<summary>
	///		Remove the Object with the specified name.
	///	</summary>
	void Remove(
		const std::string & strName
	);

	///	<summary>
	///		Assign the Object with the specified name.
	///	</summary>
	///	<returns>
	///		true if insertion is successful.  false if parent Object could
	///		not be found in the ObjectRegistry.
	///	</returns>
	bool Assign(
		const std::string & strName,
		Object * pObject
	);

private:
	///	<summary>
	///		Map from Object name to Object instance.
	///	</summary>
	ObjectMap m_mapObjects;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for representing objects.
///	</summary>
class Object {

friend class ObjectRegistry;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Object(const std::string & strName) :
		m_strName(strName)
	{ }

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Object() {
	}

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(new Object(strDuplicateName), objreg);
	}

protected:
	///	<summary>
	///		Duplicate this object and its children.
	///	</summary>
	Object * _Duplicate(
		Object * pobjDuplicate,
		ObjectRegistry & objreg
	) const;

protected:
	///	<summary>
	///		Name of the Object.
	///	</summary>
	std::string m_strName;

	///	<summary>
	///		List of child Objects.
	///	</summary>
	ObjectChildrenSet m_setChildren;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for representing String objects.
///	</summary>
class StringObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	StringObject(
		const std::string & strName,
		const std::string & strValue
	) :
		Object(strName),
		m_strValue(strValue)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new StringObject(strDuplicateName, m_strValue),
			objreg);
	}

protected:
	///	<summary>
	///		String value of the Object.
	///	</summary>
	std::string m_strValue;

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for representing integer objects.
///	</summary>
class IntegerObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	IntegerObject(
		const std::string & strName,
		const int & iValue
	) :
		Object(strName),
		m_iValue(iValue)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new IntegerObject(strDuplicateName, m_iValue),
			objreg);
	}

protected:
	///	<summary>
	///		Integer value of the Object.
	///	</summary>
	int m_iValue;

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for representing floating point objects.
///	</summary>
class FloatingPointObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FloatingPointObject(
		const std::string & strName,
		const double & dValue
	) :
		Object(strName),
		m_dValue(dValue)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new FloatingPointObject(strDuplicateName, m_dValue),
			objreg);
	}

protected:
	///	<summary>
	///		Double value of the Object.
	///	</summary>
	double m_dValue;

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for representing Variable operations.
///	</summary>
class VariableObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableObject(
		const std::string & strName,
		const std::string & strValue
	) :
		Object(strName),
		m_strValue(strValue)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new VariableObject(strDuplicateName, m_strValue),
			objreg);
	}

protected:
	///	<summary>
	///		String describing the Variable operation.
	///	</summary>
	std::string m_strValue;

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for registering Objects.
///	</summary>
class ListObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ListObject(
		const std::string & strName
	) :
		Object(strName)
	{ }

	///	<summary>
	///		Add an Object to the ListObject.
	///	</summary>
	void PushBack(
		const std::string & strObject
	) {
		m_vecObjectNames.push_back(strObject);
	}

protected:
	///	<summary>
	///		Vector of Object names in this ListObject.
	///	</summary>
	std::vector<std::string> m_vecObjectNames;

};

///////////////////////////////////////////////////////////////////////////////

#endif

