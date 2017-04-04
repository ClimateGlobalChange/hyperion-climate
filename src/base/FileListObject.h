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
#include "DataArray1D.h"
#include "GlobalFunction.h"

///////////////////////////////////////////////////////////////////////////////

class RecapConfigObject;

class GridObject;

class Variable;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A local (file,time) pair.
///	</summary>
typedef std::pair<size_t, int> LocalFileTimePair;

///	<summary>
///		A map from time indices to local (file,time) pairs.
///	</summary>
typedef std::map<size_t, LocalFileTimePair> VariableTimeFileMap;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class that describes primitive variable information from a FileList.
///	</summary>
class VariableInfo {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableInfo(
		const std::string & strVariableName
	) :
		m_strVariableName(strVariableName),
		m_iTimeDimIx(-1)
		//m_nTimeSliceDims(1)
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
	///		Dimension names.
	///	</summary>
	std::vector<std::string> m_vecDimNames;

	///	<summary>
	///		Index of time dimension or (-1) if time dimension doesn't exist.
	///	</summary>
	int m_iTimeDimIx;
/*
	///	<summary>
	///		Dimensionality of one time slice of data.
	///	</summary>
	int m_nTimeSliceDims;
*/
	///	<summary>
	///		Map from Times to filename index and time index.
	///	</summary>
	VariableTimeFileMap m_mapTimeFile;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A GlobalFunction that builds a new FileListObject.
///	</summary>
class FileListObjectConstructor : public GlobalFunction {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FileListObjectConstructor(const std::string & strName) :
		GlobalFunction(strName)
	{ }

public:
	///	<summary>
	///		Call a member function of this GlobalFunction.
	///	</summary>
	virtual std::string Call(
		const ObjectRegistry & objreg,
		const std::vector<std::string> & vecCommandLine,
		const std::vector<ObjectType> & vecCommandLineType,
		Object ** ppReturn
	);
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing a list of files.
///	</summary>
class FileListObject : public Object {

public:
	///	<summary>
	///		Invalid File index.
	///	</summary>
	static const size_t InvalidFileIx;

	///	<summary>
	///		Invalid Time index.
	///	</summary>
	static const size_t InvalidTimeIx;

	///	<summary>
	///		A value to denote that a dimension has inconsistent sizes
	///		across files.
	///	</summary>
	static const long InconsistentDimensionSizes;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FileListObject(
		const std::string & strName
	) :
		Object(strName),
		m_pobjRecapConfig(NULL),
		m_strRecordDimName("time"),
		m_sReduceTargetIx(InvalidFileIx)
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		FileListObject * pobjDuplicate =
			new FileListObject(strDuplicateName);
		if (pobjDuplicate == NULL) {
			_EXCEPTIONT("Out of memory");
		}

		pobjDuplicate->m_strRecordDimName = m_strRecordDimName;
		pobjDuplicate->m_vecFilenames = m_vecFilenames;
		pobjDuplicate->m_vecTimes = m_vecTimes;
		pobjDuplicate->m_vecVariableInfo = m_vecVariableInfo;
		pobjDuplicate->m_mapDimNameSize = m_mapDimNameSize;
		pobjDuplicate->m_sReduceTargetIx = m_sReduceTargetIx;
		pobjDuplicate->m_mapOutputTimeFile = m_mapOutputTimeFile;

		return _Duplicate(pobjDuplicate, objreg);
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
	);

	///	<summary>
	///		Set the RecapConfigObject pointer for this FileListObject.
	///	</summary>
	void SetRecapConfigObject(
		RecapConfigObject * pobjRecapConfig
	) {
		m_pobjRecapConfig = pobjRecapConfig;
	}

public:
	///	<summary>
	///		Get the count of filenames.
	///	</summary>
	size_t GetFilenameCount() const {
		return m_vecFilenames.size();
	}

	///	<summary>
	///		Get the vector of filenames.
	///	</summary>
	const std::string & GetFilename(size_t f) const {
		if (f >= m_vecFilenames.size()) {
			_EXCEPTIONT("Index out of range");
		}
		return m_vecFilenames[f];
	}

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

	///	<summary>
	///		Populate from a search string.
	///	</summary>
	std::string PopulateFromSearchString(
		const std::string & strSearchString
	);

	///	<summary>
	///		Add a single timeslice file with the given filename.
	///	</summary>
	std::string CreateFileNoTime(
		const std::string & strFilename,
		const GridObject * pobjGrid
	);

	///	<summary>
	///		Set the reduce target by filename.
	///	</summary>
	std::string SetReduceTarget(
		const std::string & strTargetFilename
	);

	///	<summary>
	///		Check if the FileList has a reduce target.
	///	</summary>
	bool HasReduceTarget() const {
		if (m_sReduceTargetIx != InvalidFileIx) {
			return true;
		}
		return false;
	}

public:
	///	<summary>
	///		Get the record dimension name.
	///	</summary>
	const std::string & GetRecordDimName() const {
		return m_strRecordDimName;
	}

	///	<summary>
	///		Get the number of time indices in the file list.
	///	</summary>
	size_t GetTimeCount() const {
		return m_vecTimes.size();
	}

	///	<summary>
	///		Get the Time with the specified index.
	///	</summary>
	const Time & GetTime(int iTime) const {
		if (iTime >= m_vecTimes.size()) {
			_EXCEPTIONT("Out of range");
		}
		return m_vecTimes[iTime];
	}

	///	<summary>
	///		Get the vector of Times associated with the FileList.
	///	</summary>
	const std::vector<Time> & GetTimes() const {
		return m_vecTimes;
	}

	///	<summary>
	///		Get the size of the specified dimension.
	///	</summary>
	long GetDimSize(const std::string & strDimName) const {
		std::map<std::string, long>::const_iterator iter =
			m_mapDimNameSize.find(strDimName);

		if (iter != m_mapDimNameSize.end()) {
			return iter->second;
		}

		_EXCEPTIONT("Invalid dimension name");
	}

	///	<summary>
	///		Check if another FileListObject has a compatible set of
	///		Times indices.
	///	</summary>
	bool IsCompatible(
		const FileListObject * pobjFileList
	);

	///	<summary>
	///		Distribute available time indices across MPI ranks.
	///	</summary>
	void GetOnRankTimeIndices(
		std::vector<size_t> & vecTimeIndices
	);

	///	<summary>
	///		Load the data from a particular variable into the given array.
	///	</summary>
	std::string LoadData_float(
		const std::string & strVariableName,
		const std::vector<long> & vecAuxIndices,
		size_t sTime,
		DataArray1D<float> & data
	);

	///	<summary>
	///		Write the data from the given array to disk.
	///	</summary>
	std::string WriteData_float(
		const std::string & strVariableName,
		const std::vector<long> & vecAuxIndices,
		size_t sTime,
		const DataArray1D<float> & data
	);

	///	<summary>
	///		Add a new variable from a template.
	///	</summary>
	std::string AddVariableFromTemplate(
		const FileListObject * pobjSourceFileList,
		const Variable * pvar,
		VariableInfo ** ppvarinfo
	);

public:
	///	<summary>
	///		Add the given dimension to this FileListObject.
	///	</summary>
	std::string AddDimension(
		const std::string & strDimName,
		long lDimSize,
		bool fGridDim = false
	);

	///	<summary>
	///		Get the size of the specified dimension.
	///	</summary>
	long GetDimensionSize(
		const std::string & strDimName
	) const;

protected:
	///	<summary>
	///		Index variable data.
	///	</summary>
	std::string IndexVariableData(
		size_t sFileIx = InvalidFileIx
	);

protected:
	///	<summary>
	///		Output the time-variable index as a CSV.
	///	</summary>
	std::string OutputTimeVariableIndexCSV(
		const std::string & strCSVOutput
	);

protected:
	///	<summary>
	///		Pointer to the associated RecapConfigObject.
	///	</summary>
	RecapConfigObject * m_pobjRecapConfig;

	///	<summary>
	///		The name of the record dimension (default "time")
	///	</summary>
	std::string m_strRecordDimName;

	///	<summary>
	///		The list of filenames.
	///	</summary>
	std::vector<std::string> m_vecFilenames;

	///	<summary>
	///		The list of Times that appear in the FileList.
	///	</summary>
	std::vector<Time> m_vecTimes;

	///	<summary>
	///		A map from Time to m_vecTimes vector index
	///	</summary>
	std::map<Time, size_t> m_mapTimeToIndex;

	///	<summary>
	///		The list of variable names that appear in the FileList.
	///	</summary>
	std::vector<VariableInfo> m_vecVariableInfo;

	///	<summary>
	///		A map between the variable name and index in m_vecVariableInfo.
	///	</summary>
	std::map<std::string, size_t> m_mapVariableNameToIndex;

	///	<summary>
	///		A map between dimension name and size in the FileList.
	///	</summary>
	std::map<std::string, long> m_mapDimNameSize;

	///	<summary>
	///		Names of grid dimensions for this FileList.
	///	</summary>
	std::vector<std::string> m_vecGridDimNames;

	///	<summary>
	///		Filename index that is the target of reductions (output mode).
	///	</summary>
	size_t m_sReduceTargetIx;

	///	<summary>
	///		Filename index for each of the time indices (output mode).
	///	</summary>
	std::map<size_t, LocalFileTimePair> m_mapOutputTimeFile;
};

///////////////////////////////////////////////////////////////////////////////

#endif

