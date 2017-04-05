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
#include "Announce.h"

#include <cstdlib>
#include <cmath>

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

std::string ObjectRegistry::Remove(
	const std::string & strName
) {
	// Check if this Object already exists
	ObjectMap::iterator iter = m_mapObjects.find(strName);
	if (iter == m_mapObjects.end()) {
		_EXCEPTION1("Object \"%s\" not found in registry", strName.c_str());
	}

	if (iter->second->m_nLocks != 0) {
		return std::string("ERROR: Cannot delete locked object \"")
			+ strName + std::string("\"");
	}

	Announce("REMOVE %s\n", strName.c_str());

	// Remove all children of this Object
	for (size_t i = 0; i < iter->second->m_vecChildren.size(); i++) {
		std::string strError =
			Remove(iter->second->m_vecChildren[i]->m_strName);

		if (strError != "") {
			return strError;
		}
	}

	// Remove this Object
	delete (iter->second);
	m_mapObjects.erase(iter);

	return std::string("");
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
		std::string strError = Remove(strName);
		if (strError != "") {
			Announce(strError.c_str());
			return false;
		}
	}

	// Assign the Object its name
	if (pObject->m_strName == "") {
		pObject->m_strName = strName;
	}

	// Add the Object to its parent
	std::string strParent;
	for (int i = strName.length()-1; i >= 0; i--) {
		if (strName[i] == '.') {
			strParent = strName.substr(0,i);
			break;
		}
	}

	// Put a pointer to the Object into the parent's array
	if (strParent != "") {
		ObjectMap::const_iterator iterParent = m_mapObjects.find(strParent);
		if (iterParent == m_mapObjects.end()) {
			return false;
		}

		//printf("PARENT %s\n", strParent.c_str());
		bool fSuccess = iterParent->second->AddChild(pObject);
		if (!fSuccess) {
			return false;
		}
	}

	// Assign the Object into the ObjectRegistry
	m_mapObjects.insert(
		ObjectMap::value_type(strName, pObject));

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Object
///////////////////////////////////////////////////////////////////////////////

std::string ObjectConstructor::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	if (ppReturn != NULL) {
		Object * pobj = new Object("");
		if (pobj == NULL) {
			_EXCEPTIONT("Unable to initialize GridObject");
		}

		// Set the return value
		(*ppReturn) = pobj;
	}

	return std::string("");
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

	for (size_t i = 0; i < m_vecChildren.size(); i++) {

		//printf("ORIGINAL CHILD %s\n", (*iter)->m_strName.c_str());
		std::string strNewChildName =
			pobjDuplicate->m_strName
			+ m_vecChildren[i]->m_strName.substr(
				nOriginalNameLength, std::string::npos);

		if (strNewChildName[nNewNameLength] != '.') {
			_EXCEPTIONT("Logic error: Invalid child");
		}

		m_vecChildren[i]->Duplicate(strNewChildName, objreg);
	}

	return (pobjDuplicate);
}

///////////////////////////////////////////////////////////////////////////////
// StringObject
///////////////////////////////////////////////////////////////////////////////

bool StringObject::ToUnit(
	const std::string & strTargetUnit,
	double * dValueOut,
	bool fIsDelta
) {
	// Check argument
	if (dValueOut == NULL) {
		_EXCEPTIONT("Invalid pointer to double on return");
	}

	// Extract the value and unit from this String
	enum ParseMode {
		ParseMode_WS,
		ParseMode_Number,
		ParseMode_Unit
	};

	ParseMode mode = ParseMode_WS;
	ParseMode modeNext = ParseMode_Number;

	bool fHasPeriod = false;

	std::string strNumber;
	std::string strUnit;

	int iPos = 0;
	for (;;) {
		if (iPos >= m_strValue.length()) {
			break;
		}

		// Whitespace
		if (mode == ParseMode_WS) {
			if ((m_strValue[iPos] == ' ') ||
			    (m_strValue[iPos] == '\t')
			) {
				iPos++;
			} else {
				mode = modeNext;
			}

		// Number
		} else if (mode == ParseMode_Number) {
			if ((m_strValue[iPos] >= '0') && (m_strValue[iPos] <= '9')) {
				strNumber += m_strValue[iPos];
				iPos++;

			} else if (m_strValue[iPos] == '.') {
				if (fHasPeriod) {
					return false;
				} else {
					strNumber += m_strValue[iPos];
					fHasPeriod = true;
					iPos++;
				}

			} else if (m_strValue[iPos] == '-') {
				if (strNumber.length() != 0) {
					return false;
				} else {
					strNumber += m_strValue[iPos];
					iPos++;
				}

			} else if (
				((m_strValue[iPos] >= 'a') && (m_strValue[iPos] <= 'z')) ||
				((m_strValue[iPos] >= 'A') && (m_strValue[iPos] <= 'Z'))
			) {
				mode = ParseMode_Unit;

			} else if (
				(m_strValue[iPos] == ' ') ||
			    (m_strValue[iPos] == '\t')
			) {
				mode = ParseMode_WS;
				modeNext = ParseMode_Unit;

			} else {
				return false;
			}

		// Unit
		} else if (mode == ParseMode_Unit) {
			if (((m_strValue[iPos] >= 'a') && (m_strValue[iPos] <= 'z')) ||
				((m_strValue[iPos] >= 'A') && (m_strValue[iPos] <= 'Z'))
			) {
				strUnit += m_strValue[iPos];
				iPos++;

			} else if (
				(m_strValue[iPos] == ' ') ||
			    (m_strValue[iPos] == '\t')
			) {
				break;

			} else {
				return false;
			}

		// Invalid mode
		} else {
			_EXCEPTIONT("Invalid mode");
		}
	}

	if (strNumber.length() == 0) {
		return false;
	}

	// Value
	const double dValue = atof(strNumber.c_str());

	// Unit is equal to TargetUnit
	if (strUnit == strTargetUnit) {
		(*dValueOut) = dValue;

	// Perform unit conversion from great circle distance (degrees)
	} else if ((strUnit == "deg") || (strUnit == "degrees_north") || (strUnit == "degrees_east")) {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			(*dValueOut) = dValue;

		} else if (strTargetUnit == "rad") {
			(*dValueOut) = dValue * M_PI / 180.0;

		} else if (strTargetUnit == "m") {
			(*dValueOut) = 6.37122e6 * dValue * M_PI / 180.0;

		} else if (strTargetUnit == "km") {
			(*dValueOut) = 6.37122e3 * dValue * M_PI / 180.0;

		} else {
			return false;
		}

	// Perform unit conversion from great circle distance (radians)
	} else if (strUnit == "rad") {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			(*dValueOut) = 180.0 / M_PI * dValue;

		} else if (strTargetUnit == "m") {
			(*dValueOut) = 6.37122e6 * dValue;

		} else if (strTargetUnit == "km") {
			(*dValueOut) = 6.37122e3 * dValue;

		} else {
			return false;
		}

	// Perform unit conversion from great circle distance (meters)
	// or altitude (meters)
	} else if (strUnit == "m") {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			(*dValueOut) = 180.0 / M_PI * dValue / 6.37122e6;

		} else if (strTargetUnit == "rad") {
			(*dValueOut) = dValue / 6.37122e6;

		} else if (strTargetUnit == "km") {
			(*dValueOut) = dValue / 1000.0;

		} else if (strTargetUnit == "m2/s2") {
			(*dValueOut) = dValue * 9.80616;

		} else {
			return false;
		}

	// Perform unit conversion from great circle distance (kilometers)
	} else if (strUnit == "km") {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			(*dValueOut) = 180.0 / M_PI * dValue / 6.37122e3;

		} else if (strTargetUnit == "rad") {
			(*dValueOut) = dValue / 6.37122e3;

		} else if (strTargetUnit == "m") {
			(*dValueOut) = dValue * 1000.0;

		} else if (strTargetUnit == "m2/s2") {
			(*dValueOut) = dValue * 1000.0 * 9.80616;

		} else {
			return false;
		}

	// Perform unit conversion from temperature (K)
	} else if (strUnit == "K") {
		if (strTargetUnit == "degC") {
			if (fIsDelta) {
				(*dValueOut) = dValue;
			} else {
				(*dValueOut) = dValue - 273.15;
			}

		} else {
			return false;
		}

	// Perform unit conversion from temperature (degC)
	} else if (strUnit == "degC") {
		if (strTargetUnit == "K") {
			if (fIsDelta) {
				(*dValueOut) = dValue;
			} else {
				(*dValueOut) = dValue + 273.15;
			}

		} else {
			return false;
		}

	// Perform unit conversion from pressure (Pa)
	} else if (strUnit == "Pa") {
		if ((strTargetUnit == "hPa") ||
		    (strTargetUnit == "mb") ||
		    (strTargetUnit == "mbar")
		) {
			(*dValueOut) = dValue / 100.0;

		} else if (strTargetUnit == "atm") {
			(*dValueOut) = dValue / 101325.0;

		} else {
			return false;
		}

	// Perform unit conversion from pressure (hPa,mb,mbar)
	} else if ((strUnit == "hPa") || (strUnit == "mb") || (strUnit == "mbar")) {
		if (strTargetUnit == "Pa") {
			(*dValueOut) = dValue * 100.0;

		} else if (
		    (strTargetUnit == "hPa") ||
		    (strTargetUnit == "mb") ||
		    (strTargetUnit == "mbar")
		) {
			(*dValueOut) = dValue;

		} else if (strTargetUnit == "atm") {
			(*dValueOut) = dValue / 1013.25;

		} else {
			return false;
		}

	// Perform unit conversion from pressure (atm)
	} else if (strUnit == "atm") {
		if (strTargetUnit == "Pa") {
			(*dValueOut) = dValue * 101325.0;

		} else if (
		    (strTargetUnit == "hPa") ||
		    (strTargetUnit == "mb") ||
		    (strTargetUnit == "mbar")
		) {
			(*dValueOut) = dValue * 1013.25;

		} else {
			return false;
		}

	// Perform unit conversion from geopotential (m2/s2)
	} else if (strUnit == "m2/s2") {
		if (strTargetUnit == "m") {
			(*dValueOut) = dValue / 9.80616;

		} else if (strUnit == "km") {
			(*dValueOut) = dValue / 9.80616 / 1000.0;

		} else {
			return false;
		}

	} else {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

