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

class VariableRegistry;

class Object;

typedef std::vector<Object *> ObjectChildrenVector;

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
	///		Get the name of this Object.
	///	</summary>
	const std::string & Name() const {
		return m_strName;
	}

	///	<summary>
	///		Create a string containing this Object's name.
	///	</summary>
	std::string ChildName(
		const std::string & strChild
	) const {
		return (m_strName + "." + strChild);
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
		return std::string("ERROR: Unknown member function [")
			+ strFunctionName + std::string("]");
	}

	///	<summary>
	///		Number of children in Object.
	///	</summary>
	size_t ChildrenCount() const {
		return m_vecChildren.size();
	}

	///	<summary>
	///		Get child Object by index.
	///	</summary>
	Object * GetChild(size_t sChild) const {
		if (sChild >= m_vecChildren.size()) {
			_EXCEPTIONT("Children vector access out of range");
		}
		return m_vecChildren[sChild];
	}

	///	<summary>
	///		Get child Object by name.
	///	</summary>
	Object * GetChild(const std::string & strChildName) const {
		std::string strFullChildName = ChildName(strChildName);
		for (size_t i = 0; i < m_vecChildren.size(); i++) {
			if (m_vecChildren[i]->Name() == strFullChildName) {
				return m_vecChildren[i];
			}
		}
		return NULL;
	}

protected:
	///	<summary>
	///		Add a child Object.
	///	</summary>
	virtual bool AddChild(Object * pChild) {
		m_vecChildren.push_back(pChild);
		return true;
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
	ObjectChildrenVector m_vecChildren;
};

///	<summary>
///		A class for representing Objects that may be distributed
///		across MPI ranks.
///	</summary>
class DistributedObject : public Object {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DistributedObject(const std::string & strName) :
		Object(strName),
		m_fDistributed(true)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new DistributedObject(strDuplicateName),
			objreg);
	}

	///	<summary>
	///		Set the distributed flag to true.
	///	</summary>
	void SetDistributed() {
		m_fDistributed = true;
	}

	///	<summary>
	///		Set the distributed flag to false.
	///	</summary>
	void UnsetDistributed() {
		m_fDistributed = false;
	}

	///	<summary>
	///		Get the state of the distributed flag.
	///	</summary>
	bool IsDistributed() const {
		return m_fDistributed;
	}

protected:
	///	<summary>
	///		A flag indicating data remains distributed across ranks.
	///	</summary>
	bool m_fDistributed;
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

	///	<summary>
	///		Get the string.
	///	</summary>
	const std::string & Value() const {
		return m_strValue;
	}

public:
	///	<summary>
	///		Return the value of the string with the specified units.
	///	</summary>
	bool ToUnit(
		const std::string & strUnit,
		double * dValueOut,
		bool fIsDelta = false
	);

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

	///	<summary>
	///		Get the integer value.
	///	</summary>
	const int & Value() const {
		return m_iValue;
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
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new ListObject(strDuplicateName),
			objreg);
	}

public:
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

