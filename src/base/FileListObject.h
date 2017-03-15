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
#include "TimeObj.h"

///////////////////////////////////////////////////////////////////////////////

class VariableInfo {

public:
	///	<summary>
	///		A local (file,time) pair.
	///	</summary>
	typedef std::pair<size_t, int> LocalFileTimePair;

	///	<summary>
	///		A map from time indices to local (file,time) pairs.
	///	</summary>
	typedef std::map<size_t, LocalFileTimePair> VariableTimeFileMap;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableInfo(
		const std::string & strVariableName
	) :
		m_strVariableName(strVariableName),
		m_iTimeDimIx(-1)
	{ } 

public:
	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strVariableName;

	///	<summary>
	///		Units for the Variable.
	///	</summary>
	std::string m_strUnits;

	///	<summary>
	///		Dimension of time variable.
	///	</summary>
	std::vector<int> m_vecDimSizes;

	///	<summary>
	///		Index of time dimension or (-1) if time dimension doesn't exist.
	///	</summary>
	int m_iTimeDimIx;

	///	<summary>
	///		Map from Times to filename index and time index.
	///	</summary>
	VariableTimeFileMap m_mapTimeFile;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing a list of files.
///	</summary>
class FileListObject : public Object {

public:
	///	<summary>
	///		Invalid Time index.
	///	</summary>
	static const size_t InvalidTimeIx;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FileListObject(
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
			new FileListObject(strDuplicateName),
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
			return OutputTimeVariableIndexCSV(vecCommandLine[0]);
		}
		return
			Object::Call(
				objreg,
				strFunctionName,
				vecCommandLine,
				vecCommandLineType,
				ppReturn);
	}

public:
	///	<summary>
	///		Get the VariableInfo associated with a given variable name.
	///	</summary>
	const VariableInfo * GetVariableInfo(
		const std::string & strVariableName
	) const {
		for (size_t i = 0; i < m_vecVariableInfo.size(); i++) {
			if (m_vecVariableInfo[i].m_strVariableName == strVariableName) {
				return &(m_vecVariableInfo[i]);
			}
		}
		return NULL;
	}

	///	</summary>
	///	<summary>
	///		Populate from a search string.
	///	</summary>
	///	<returns>
	///		An error message if an error occurred.  A blank string otherwise.
	///	</returns>
	std::string PopulateFromSearchString(
		const std::string & strSearchString,
		ObjectRegistry & objreg
	);

protected:
	///	<summary>
	///		Index variable data.
	///	</summary>
	std::string IndexVariableData();

protected:
	///	<summary>
	///		Output the time-variable index as a CSV.
	///	</summary>
	std::string OutputTimeVariableIndexCSV(
		const std::string & strCSVOutput
	);

protected:
	///	<summary>
	///		The list of filenames.
	///	</summary>
	std::vector<std::string> m_vecFilenames;

	///	<summary>
	///		The list of Times that appear in the FileList.
	///	</summary>
	std::vector<Time> m_vecTimes;

	///	<summary>
	///		The list of variable names that appear in the FileList.
	///	</summary>
	std::vector<VariableInfo> m_vecVariableInfo;

};

///////////////////////////////////////////////////////////////////////////////

#endif

