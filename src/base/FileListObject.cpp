///////////////////////////////////////////////////////////////////////////////
///
///	\file    FileListObject.cpp
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

#include "FileListObject.h"
#include "RecapConfigObject.h"
#include "STLStringHelper.h"
#include "GridElements.h"
#include "GridObject.h"
#include "Variable.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "netcdfcpp.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <fstream>

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////
// FileListObjectConstructor
///////////////////////////////////////////////////////////////////////////////

std::string FileListObjectConstructor::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	FileListObject * pobjFileList = new FileListObject("");
	if (pobjFileList == NULL) {
		_EXCEPTIONT("Unable to initialize FileListObject");
	}

	// Constructor accepts 0 or 1 arguments
	if (vecCommandLine.size() > 1) {
		return std::string("ERROR: Invalid arguments to file_list()");
	}
	if (vecCommandLine.size() == 1) {
		if (vecCommandLineType[0] != ObjectType_String) {
			return std::string("ERROR: Invalid first argument [filename] "
				"in file_list() call");
		}

		// Initialize the FileList with the given search string
		std::string strError =
			pobjFileList->PopulateFromSearchString(
				vecCommandLine[0]);

		if (strError != "") {
			return strError;
		}
	}

	// Set the return value
	if (ppReturn != NULL) {
		(*ppReturn) = pobjFileList;
	} else {
		delete pobjFileList;
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////
// FileListObject
///////////////////////////////////////////////////////////////////////////////

const size_t FileListObject::InvalidFileIx = (-1);

///////////////////////////////////////////////////////////////////////////////

const size_t FileListObject::InvalidTimeIx = (-1);

///////////////////////////////////////////////////////////////////////////////

const long FileListObject::InconsistentDimensionSizes = (-1);

///////////////////////////////////////////////////////////////////////////////

FileListObject::~FileListObject() {
	for (int v = 0; v < m_vecVariableInfo.size(); v++) {
		delete m_vecVariableInfo[v];
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::Call(
	const ObjectRegistry & objreg,
	const std::string & strFunctionName,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	// Output information about the FileList to a CSV file
	if (strFunctionName == "output_csv") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"output_csv\"");
		}
		return OutputTimeVariableIndexCSV(vecCommandLine[0]);
	}

	// Create a copy of the FileList with modified filenames
	if (strFunctionName == "duplicate_for_writing") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"duplicate_for_writing\"");
		}

		FileListObject * pobjNewFileList = new FileListObject("");
		if (pobjNewFileList == NULL) {
			return std::string("ERROR: Unable to allocate new FileListObject");
		}

		pobjNewFileList->m_strBaseDir = vecCommandLine[0];
		if (vecCommandLine[0][vecCommandLine[0].length()-1] != '/') {
			pobjNewFileList->m_strBaseDir += "/";
		}

		// Create the directory if it doesn't exist
		DIR * pDir = opendir(pobjNewFileList->m_strBaseDir.c_str());
		if (pDir == NULL) {
			int iError =
				mkdir(
					pobjNewFileList->m_strBaseDir.c_str(),
					S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

			if (iError != 0) {
				return std::string("Unable to create directory \"")
					+ pobjNewFileList->m_strBaseDir + std::string("\"");
			}

		} else {
			closedir(pDir);
		}

		// Copy times
		pobjNewFileList->m_vecTimes = m_vecTimes;
		pobjNewFileList->m_mapTimeToIndex = m_mapTimeToIndex;
		pobjNewFileList->m_strTimeUnits = m_strTimeUnits;

/*
		// Add file extension
		for (size_t f = 0; f < pobjNewFileList->m_vecFilenames.size(); f++) {
			int nLength = pobjNewFileList->m_vecFilenames[f].length();
			int iExt = nLength-1;
			for (; iExt >= 0; iExt--) {
				if (pobjNewFileList->m_vecFilenames[f][iExt] == '.') {
					break;
				}
			}
			if (iExt == (-1)) {
				pobjNewFileList->m_vecFilenames[f] =
					m_vecFilenames[f]
					+ vecCommandLine[0];
			} else {
				pobjNewFileList->m_vecFilenames[f] =
					m_vecFilenames[f].substr(0,iExt)
					+ vecCommandLine[0]
					+ m_vecFilenames[f].substr(iExt);
			}
		}
*/
		if (ppReturn != NULL) {
			(*ppReturn) = pobjNewFileList;
		} else {
			delete pobjNewFileList;
		}
		return std::string("");
	}

	// Append files to the file_list
	if (strFunctionName == "append") {
		if ((vecCommandLineType.size() != 1) ||
		    (vecCommandLineType[0] != ObjectType_String)
		) {
			return std::string("ERROR: Invalid parameters to function \"add_file\"");
		}

		if (IsLocked()) {
			return std::string("ERROR: Cannot append files to a locked file_list");
		}

		return PopulateFromSearchString(m_strBaseDir + vecCommandLine[0]);
	}

	return
		Object::Call(
			objreg,
			strFunctionName,
			vecCommandLine,
			vecCommandLineType,
			ppReturn);
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::PopulateFromSearchString(
	const std::string & strSearchString
) {
	// Check if already initialized
	if (m_vecFilenames.size() != 0) {
		_EXCEPTIONT("FileListObject has already been initialized");
	}

	// File the directory in the search string
	std::string strFileSearchString;
	for (int i = strSearchString.length(); i >= 0; i--) {
		if (strSearchString[i] == '/') {
			m_strBaseDir = strSearchString.substr(0,i+1);
			strFileSearchString =
				strSearchString.substr(i+1, std::string::npos);
			break;
		}
	}
	if ((m_strBaseDir == "") && (strFileSearchString == "")) {
		strFileSearchString = strSearchString;
		m_strBaseDir = "./";
	}

	// Open the directory
	DIR * pDir = opendir(m_strBaseDir.c_str());
	if (pDir == NULL) {
		return std::string("Unable to open directory \"")
			+ m_strBaseDir + std::string("\"");
	}

	// Search all files in the directory for match to search string
	size_t iFileBegin = m_vecFilenames.size();
	struct dirent * pDirent;
	while ((pDirent = readdir(pDir)) != NULL) {
		std::string strFilename = pDirent->d_name;
		if (STLStringHelper::WildcardMatch(
				strFileSearchString.c_str(),
				strFilename.c_str())
		) {
			// File found, insert into list of filenames
			m_vecFilenames.push_back(strFilename);
		}
	}
	closedir(pDir);

	// Index the variable data
	return IndexVariableData(iFileBegin, m_vecFilenames.size());
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::CreateFilesFromTemplate(
	const std::string & strFilenameTemplate,
	const GridObject * pobjGrid,
	int nTimesPerFile
) {
	if (pobjGrid == NULL) {
		_EXCEPTION();
	}

	// Search for time specifier in template
	int iTemplatePos = (-1);
	for (int i = 0; i < strFilenameTemplate.length()-1; i++) {
		if ((strFilenameTemplate[i] == '%') &&
			(strFilenameTemplate[i+1] == 'T')
		) {
			if (iTemplatePos == (-1)) {
				iTemplatePos = i;
			} else {
				return std::string("ERROR: Time specifier appears more than once in filename template");
			}
		}
	}
	if (iTemplatePos == (-1)) {
		return std::string("ERROR: Time specifier missing from filename template");
	}

	// Loop through all times
	for (int t = 0; t < m_vecTimes.size(); t += nTimesPerFile) {
		std::string strFilename =
			strFilenameTemplate.substr(0,iTemplatePos)
			+ m_vecTimes[t].ToShortString()
			+ strFilenameTemplate.substr(iTemplatePos+2, std::string::npos);

		// Create the file from the Grid
		std::string strError =
			CreateFileNoTime(
				strFilename,
				pobjGrid);

		if (strError != "") {
			return strError;
		}

		// Add record dimension
		std::string strFullFilename = m_strBaseDir + strFilename;
		NcFile ncFile(strFullFilename.c_str(), NcFile::Write);
		if (!ncFile.is_valid()) {
			return std::string("ERROR: Unable to open \"")
				+ strFilename + std::string("\" for writing");
		}

		if (m_strRecordDimName == "") {
			return std::string("ERROR: NetCDF file template has no record dimension");
		}

		NcDim * dimRecord = ncFile.add_dim(m_strRecordDimName.c_str(), 0);
		if (dimRecord == NULL) {
			return std::string("ERROR: Unable to create record dimension \"")
				+ m_strRecordDimName + std::string("\" in file ")
				+ strFilename;
		}

		NcVar * varRecord = ncFile.add_var(m_strRecordDimName.c_str(), ncDouble, dimRecord);
		if (varRecord == NULL) {
			return std::string("ERROR: Unable to create record variable \"")
				+ m_strRecordDimName + std::string("\" in file ")
				+ strFilename;
		}

		NcBool fAttUnits = varRecord->add_att("units", m_strTimeUnits.c_str());
		if (!fAttUnits) {
			return std::string("ERROR: Unable to add attributes to variable \"")
				+ m_strRecordDimName + std::string("\" in file ")
				+ strFilename;
		}

		// Output times
		int nTimes = nTimesPerFile;
		if (t + nTimes >= m_vecTimes.size()) {
			nTimes = m_vecTimes.size() - t;
		}

		std::vector<double> dTimes;
		dTimes.resize(nTimes);
		for (int s = 0; s < nTimes; s++) {
			dTimes[s] = m_vecTimes[t+s].GetCFCompliantUnitsOffsetDouble(m_strTimeUnits);
		}
		varRecord->put(&(dTimes[0]), nTimes);

		// Populate m_mapOutputTimeFile
		for (int s = 0; s < nTimes; s++) {
			m_mapOutputTimeFile.insert(
				std::pair<size_t, LocalFileTimePair>(
					t+s, LocalFileTimePair(t/nTimesPerFile, s)));
		}
	}
	
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::CreateFileNoTime(
	const std::string & strFilename,
	const GridObject * pobjGrid
) {
	if (pobjGrid == NULL) {
		_EXCEPTION();
	}

	// Check if file already exists
	for (size_t f = 0; f < m_vecFilenames.size(); f++) {
		if (m_vecFilenames[f] == strFilename) {
			return std::string("ERROR: File \"") + strFilename
				+ std::string("\" already exists in file_list");
		}
	}

	// Add the file to the registry
	size_t sNewFileIx = m_vecFilenames.size();
	m_vecFilenames.push_back(strFilename);

	std::string strFullFilename = m_strBaseDir + strFilename;

	NcFile ncFile(strFullFilename.c_str(), NcFile::Replace);
	if (!ncFile.is_valid()) {
		return std::string("ERROR: Unable to open \"")
			+ strFullFilename + std::string("\" for writing");
	}

	// Add grid dimensions to file
	const Mesh & mesh = pobjGrid->GetMesh();
	for (size_t i = 0; i < mesh.vecDimNames.size(); i++) {
		NcDim * dim =
			ncFile.add_dim(
				mesh.vecDimNames[i].c_str(),
				mesh.vecDimSizes[i]);

		if (dim == NULL) {
			_EXCEPTION1("Unable to create dimension \"%s\" in file",
				mesh.vecDimNames[i].c_str());
		}

		if (mesh.vecDimValues.size() > i) {
			if (mesh.vecDimValues[i].size() != mesh.vecDimSizes[i]) {
				_EXCEPTIONT("Mesh DimValues / DimSizes mismatch");
			}

			NcVar * var =
				ncFile.add_var(
					mesh.vecDimNames[i].c_str(),
					ncDouble,
					dim);

			if (var == NULL) {
				_EXCEPTION1("Unable to create variable \"%s\" in file",
					mesh.vecDimNames[i].c_str());
			}

			var->set_cur((long)0);
			var->put(&(mesh.vecDimValues[i][0]), mesh.vecDimSizes[i]);
		}
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::SetReduceTarget(
	const std::string & strTargetFilename
) {

	// Check if file exists
	for (size_t f = 0; f < m_vecFilenames.size(); f++) {
		if (m_vecFilenames[f] == strTargetFilename) {
			m_sReduceTargetIx = f;
			return std::string("");
		}
	}

	return std::string("ERROR: Reduce target not in file_list");
}

///////////////////////////////////////////////////////////////////////////////

bool FileListObject::IsCompatible(
	const FileListObject * pobjFileList
) {
	if (m_vecTimes.size() != pobjFileList->m_vecTimes.size()) {
		return false;
	}
	for (size_t t = 0; t < m_vecTimes.size(); t++) {
		if (m_vecTimes[t] != pobjFileList->m_vecTimes[t]) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

void FileListObject::GetOnRankTimeIndices(
	std::vector<size_t> & vecTimeIndices,
	size_t sTimeStride
) {
	if ((sTimeStride == 0) || (sTimeStride > 1000)) {
		_EXCEPTIONT("timestride out of range");
	}

#if defined(HYPERION_MPIOMP)
	size_t sTimeCount = m_vecTimes.size() / sTimeStride;

	int nCommRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nCommRank);

	int nCommSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);

	int nMaxTimesPerRank = static_cast<int>(sTimeCount) / nCommSize;
	if (sTimeCount % nCommSize != 0) {
		nMaxTimesPerRank++;
	}

	int iBegin = nMaxTimesPerRank * nCommRank;
	int iEnd = nMaxTimesPerRank * (nCommRank + 1);
	if (iEnd > sTimeCount) {
		iEnd = sTimeCount;
	}

	for (size_t i = iBegin; i < iEnd; i++) {
		vecTimeIndices.push_back(i * sTimeStride);
	}
#else
	for (size_t i = 0; i < m_vecTimes.size(); i += sTimeStride) {
		vecTimeIndices.push_back(i);
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::LoadData_float(
	const std::string & strVariableName,
	const std::vector<long> & vecAuxIndices,
	//size_t sTime,
	DataArray1D<float> & data
) {
	// Find the VariableInfo structure for this Variable
	size_t iVarInfo = 0;
	for (; iVarInfo < m_vecVariableInfo.size(); iVarInfo++) {
		if (strVariableName == m_vecVariableInfo[iVarInfo]->m_strVariableName) {
			break;
		}
	}
	if (iVarInfo == m_vecVariableInfo.size()) {
		_EXCEPTION1("Variable \"%s\" not found in file_list index",
			strVariableName.c_str());
	}

	// Find the local file/time pair associated with this global time index
	VariableInfo & varinfo = *(m_vecVariableInfo[iVarInfo]);

	// Get the time index
	if ((varinfo.m_iTimeDimIx != (-1)) &&
		(varinfo.m_iTimeDimIx >= vecAuxIndices.size())
	) {
		_EXCEPTIONT("time index exceeds auxiliary index size");
	}

	// Extract the global time index
	size_t sTime = (-1);
	if (varinfo.m_iTimeDimIx != (-1)) {
		sTime = vecAuxIndices[varinfo.m_iTimeDimIx];
	}

	// Find local file/time index
	VariableTimeFileMap::const_iterator iter =
		varinfo.m_mapTimeFile.find(sTime);

	if (iter == varinfo.m_mapTimeFile.end()) {
		_EXCEPTION2("sTime (%s) (%lu) not found", strVariableName.c_str(), sTime);
	}

	size_t sFile = iter->second.first;
	int iTime = iter->second.second;

	{
		std::string strLoading =
			std::string("Loading data [") + strVariableName + std::string("]");

		for (size_t d = 0; d < vecAuxIndices.size(); d++) {
			if (d == varinfo.m_iTimeDimIx) {
				strLoading +=
					std::string(" [")
					+ m_vecTimes[sTime].ToString()
					+ std::string("]");
			} else {
				strLoading +=
					std::string(" [")
					+ varinfo.m_vecAuxDimNames[d]
					+ std::string(": ")
					+ std::to_string(vecAuxIndices[d])
					+ std::string("]");
			}
		}
		Announce(strLoading.c_str());
	}

	// Open the correct NetCDF file
	std::string strFullFilename = m_strBaseDir + m_vecFilenames[sFile];
	NcFile ncfile(strFullFilename.c_str());
	if (!ncfile.is_valid()) {
		_EXCEPTION1("Cannot open file \"%s\"", strFullFilename.c_str());
	}

	// Get the correct variable from the file
	NcVar * var = ncfile.get_var(strVariableName.c_str());
	const long nDims = var->num_dims() - m_vecGridDimNames.size();

	if (var == NULL) {
		_EXCEPTION1("Variable \"%s\" no longer found in file",
			strVariableName.c_str());
	}
	if ((var->type() != ncFloat) && (var->type() != ncDouble)) {
		return std::string("Variable \"")
			+ strVariableName
			+ std::string("\" is not of type float or double");
	}
	if (nDims != vecAuxIndices.size()) {
		_EXCEPTION2("Auxiliary index array size mismatch (%li / %lu)",
			nDims, vecAuxIndices.size());
	}

	// Set the data position and size
	long lTotalSize = 1;
	std::vector<long> vecPos = vecAuxIndices;
	std::vector<long> vecSize = vecAuxIndices;

	for (size_t d = 0; d < vecAuxIndices.size(); d++) {
		if (d == varinfo.m_iTimeDimIx) {
			vecPos[d] = iTime;
		}
		vecSize[d] = 1;
	}
	for (size_t d = 0; d < m_vecGridDimNames.size(); d++) {
		NcDim * dimGrid = var->get_dim(vecPos.size());
		vecPos.push_back(0);
		vecSize.push_back(dimGrid->size());
		lTotalSize *= dimGrid->size();
	}

	if (data.GetRows() != lTotalSize) {
		//for (int d = 0; d < vecAuxIndices.size(); d++) {
		//	printf("%s %i\n", varinfo.m_vecDimNames[d].c_str(), vecAuxIndices[d]);
		//}
		//for (int d = 0; d < nDims; d++) {
		//	printf("%i %i\n", vecPos[d], vecSize[d]);
		//}
		_EXCEPTION2("Data size mismatch (%i/%lu)", data.GetRows(), lTotalSize);
	}

	// Set the position
	var->set_cur(&(vecPos[0]));

	// Load the data
/*
	std::cout << varinfo.m_strVariableName << std::endl;
		for (int d = 0; d < vecAuxIndices.size(); d++) {
			printf("%s %i\n", varinfo.m_vecDimNames[d].c_str(), vecAuxIndices[d]);
		}
		for (int d = 0; d < nDims; d++) {
			printf("%i %i\n", vecPos[d], vecSize[d]);
		}
*/
	if (var->type() == ncDouble) {
		DataArray1D<double> data_dbl(data.GetRows());
		var->get(&(data_dbl[0]), &(vecSize[0]));
		for (int i = 0; i < data.GetRows(); i++) {
			data[i] = static_cast<float>(data_dbl[i]);
		}
	} else {
		var->get(&(data[0]), &(vecSize[0]));
	}

	NcError err;
	if (err.get_err() != NC_NOERR) {
		_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
	}

	// Cleanup
	ncfile.close();

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::WriteData_float(
	const std::string & strVariableName,
	const std::vector<long> & vecAuxIndices,
	const DataArray1D<float> & data
) {
	// Find the VariableInfo structure for this Variable
	size_t iVarInfo = 0;
	for (; iVarInfo < m_vecVariableInfo.size(); iVarInfo++) {
		if (strVariableName == m_vecVariableInfo[iVarInfo]->m_strVariableName) {
			break;
		}
	}

	// File index for write
	size_t sFile = (-1);
	int iLocalTime = (-1);

	// Not found
	if (iVarInfo == m_vecVariableInfo.size()) {
		_EXCEPTION();
	}

	// Get file index from varinfo
	VariableInfo & varinfo = *(m_vecVariableInfo[iVarInfo]);

	// Extract the time index
	size_t sTime = (-1);
	if (varinfo.m_iTimeDimIx != (-1)) {
		sTime = vecAuxIndices[varinfo.m_iTimeDimIx];
	}

	// Write message
	{
		std::string strWriting =
			std::string("Writing data [") + strVariableName + std::string("]");

		for (size_t d = 0; d < vecAuxIndices.size(); d++) {
			if (d == varinfo.m_iTimeDimIx) {
				strWriting +=
					std::string(" [")
					+ m_vecTimes[sTime].ToString()
					+ std::string("]");
			} else {
				strWriting +=
					std::string(" [")
					+ varinfo.m_vecAuxDimNames[d]
					+ std::string(": ")
					+ std::to_string(vecAuxIndices[d])
					+ std::string("]");
			}
		}
		Announce(strWriting.c_str());
	}

	// Check consistency with indexing procedure
	if ((sTime == (-1)) && (varinfo.m_iTimeDimIx != (-1))) {
		_EXCEPTIONT("Attempting to write single time index "
			"for multi-time-index variable");
	}
	if ((sTime != (-1)) && (varinfo.m_iTimeDimIx == (-1))) {
		_EXCEPTIONT("Attempting to write single time indexed "
			"variable as multi-indexed variable");
	}

	VariableTimeFileMap::const_iterator iterTimeFile =
		varinfo.m_mapTimeFile.find(sTime);

	if (iterTimeFile == varinfo.m_mapTimeFile.end()) {
		if (sTime == InvalidTimeIx) {
			sFile = m_sReduceTargetIx;
			iLocalTime = (-1);

		} else {
			std::map<size_t, LocalFileTimePair>::const_iterator iterFileIx =
				m_mapOutputTimeFile.find(sTime);

			if (iterFileIx == m_mapOutputTimeFile.end()) {
				_EXCEPTIONT("Unable to determine output file");
			}

			sFile = iterFileIx->second.first;
			iLocalTime = iterFileIx->second.second;

			varinfo.m_mapTimeFile.insert(
				VariableTimeFileMap::value_type(
					sTime,
					LocalFileTimePair(sFile, sTime)));
		}

	} else {
		sFile = iterTimeFile->second.first;
		iLocalTime = iterTimeFile->second.second;
	}

	if (sFile == (-1)) {
		_EXCEPTIONT("Logic error");
	}

	// Set the data position and size
	long lTotalSize = 1;
	std::vector<long> vecPos = vecAuxIndices;
	std::vector<long> vecSize = vecAuxIndices;

	for (size_t d = 0; d < vecAuxIndices.size(); d++) {
		if (d == varinfo.m_iTimeDimIx) {
			vecPos[d] = iLocalTime;
		}
		vecSize[d] = 1;
	}
	for (size_t d = 0; d < m_vecGridDimNames.size(); d++) {
		std::map<std::string, long>::const_iterator iter =
			m_mapDimNameSize.find(m_vecGridDimNames[d]);

		if (iter == m_mapDimNameSize.end()) {
			_EXCEPTIONT("Dimension not found in map");
		}

		vecPos.push_back(0);
		vecSize.push_back(iter->second);
		lTotalSize *= iter->second;
	}

	if (data.GetRows() != lTotalSize) {
		_EXCEPTION2("Data size mismatch (%i/%lu)", data.GetRows(), lTotalSize);
	}

	// Write data
	std::string strFullFilename = m_strBaseDir + m_vecFilenames[sFile];
	NcFile ncout(strFullFilename.c_str(), NcFile::Write);
	if (!ncout.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\"", strFullFilename.c_str());
	}

	// Get dimensions
	std::vector<NcDim *> vecDims;
	vecDims.resize(varinfo.m_vecDimNames.size());
	for (int d = 0; d < vecDims.size(); d++) {
		long lDimSize = GetDimensionSize(varinfo.m_vecDimNames[d]);
		vecDims[d] = ncout.get_dim(varinfo.m_vecDimNames[d].c_str());
		if (vecDims[d] != NULL) {
			if (vecDims[d]->size() != lDimSize) {
				_EXCEPTION1("Dimension %s mismatch",
					varinfo.m_vecDimNames[d].c_str());
			}
		} else {
			vecDims[d] =
				ncout.add_dim(
					varinfo.m_vecDimNames[d].c_str(),
					lDimSize);

			if (vecDims[d] == NULL) {
				_EXCEPTION1("Cannot add dimension %s",
					varinfo.m_vecDimNames[d].c_str());
			}
		}
	}

	// Create variable
	NcVar * var = ncout.get_var(strVariableName.c_str());
	if (var == NULL) {
		var = ncout.add_var(
			strVariableName.c_str(),
			ncFloat,
			vecDims.size(),
			(const NcDim **)(&(vecDims[0])));

		if (var == NULL) {
			_EXCEPTION1("Unable to create variable \"%s\"",
				strVariableName.c_str());
		}

		var->add_att("units", varinfo.m_strUnits.c_str());
	}

	// Set current position
	var->set_cur(&(vecPos[0]));

	// Write data
	var->put(&(data[0]), &(vecSize[0]));

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::AddVariableFromTemplate(
	const FileListObject * pobjSourceFileList,
	const Variable * pvar,
	VariableInfo ** ppvarinfo
) {
	// Check arguments
	if (ppvarinfo == NULL) {
		_EXCEPTIONT("Invalid value for \"ppvarinfo\"");
	}

	// Get the FileList object from the other RecapConfig
	if (GetRecordDimName() != pobjSourceFileList->GetRecordDimName()) {
		_EXCEPTIONT("Record dim name mismatch");
	}

	// Check if variable already exists
	for (int v = 0; v < m_vecVariableInfo.size(); v++) {
		if (m_vecVariableInfo[v]->m_strVariableName == pvar->Name()) {
			return std::string("ERROR: Variable already exists in file_list");
		}
	}

	// VariableInfo
	VariableInfo * pvarinfo = new VariableInfo(pvar->Name());
	if (pvarinfo == NULL) {
		_EXCEPTIONT("Unable to allocate VariableInfo");
	}

	pvarinfo->m_strUnits = pvar->Units();

	// Add dimensions from variable to output FileList
	const std::vector<std::string> & vecDimNames = pvar->AuxDimNames();
	for (int d = 0; d < vecDimNames.size(); d++) {
		bool fGridDim = false;
		for (int e = 0; e < pobjSourceFileList->m_vecGridDimNames.size(); e++) {
			if (vecDimNames[d] == pobjSourceFileList->m_vecGridDimNames[e]) {
				fGridDim = true;
				break;
			}
		}
		if (fGridDim) {
			continue;
		}

		long lDimSize = pobjSourceFileList->GetDimensionSize(vecDimNames[d]);
		if (lDimSize == (-1)) {
			_EXCEPTIONT("Logic error");
		}

		pvarinfo->m_vecDimNames.push_back(vecDimNames[d]);
		pvarinfo->m_vecDimSizes.push_back(lDimSize);
		pvarinfo->m_vecAuxDimNames.push_back(vecDimNames[d]);
		pvarinfo->m_vecAuxDimSizes.push_back(lDimSize);

		AddDimension(
			vecDimNames[d], 
			lDimSize);
	}

	// Add grid dimensions
	for (int d = 0; d < m_vecGridDimNames.size(); d++) {
		pvarinfo->m_vecDimNames.push_back(m_vecGridDimNames[d]);
	}

	// Identify the time dimension
	pvarinfo->m_iTimeDimIx = (-1);
	for (int d = 0; d < pvarinfo->m_vecDimNames.size(); d++) {
		if (pvarinfo->m_vecDimNames[d] == m_strRecordDimName) {
			pvarinfo->m_iTimeDimIx = d;
			break;
		}
	}

	// Add this VariableInfo to the vector of VariableInfos
	m_vecVariableInfo.push_back(pvarinfo);

	(*ppvarinfo) = pvarinfo;

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::AddDimension(
	const std::string & strDimName,
	long lDimSize,
	bool fGridDim
) {
	std::map<std::string, long>::const_iterator iter =
		m_mapDimNameSize.find(strDimName);

	if (iter != m_mapDimNameSize.end()) {
		if (iter->second != lDimSize) {
			_EXCEPTION3("Dimension size mismatch: %s, %li, %li",
				strDimName.c_str(), lDimSize, iter->second);
		}

	} else {
		m_mapDimNameSize.insert(
			std::pair<std::string, long>(
				strDimName, lDimSize));
	}

	// Add a grid dimension
	if (fGridDim) {
		m_vecGridDimNames.push_back(strDimName);
		for (size_t v = 0; v < m_vecVariableInfo.size(); v++) {
			if (m_vecVariableInfo[v]->m_vecDimNames.size() !=
				m_vecVariableInfo[v]->m_vecDimSizes.size()
			) {
				_EXCEPTIONT("Logic error");
			}
			m_vecVariableInfo[v]->m_vecAuxDimNames.clear();
			m_vecVariableInfo[v]->m_vecAuxDimSizes.clear();
			for (size_t d = 0; d < m_vecVariableInfo[v]->m_vecDimSizes.size(); d++) {
				if (m_vecVariableInfo[v]->m_vecDimNames[d] != strDimName) {
					m_vecVariableInfo[v]->m_vecAuxDimNames.push_back(
						m_vecVariableInfo[v]->m_vecDimNames[d]);
					m_vecVariableInfo[v]->m_vecAuxDimSizes.push_back(
						m_vecVariableInfo[v]->m_vecDimSizes[d]);
				}
			}
		}
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

long FileListObject::GetDimensionSize(
	const std::string & strDimName
) const {
	std::map<std::string, long>::const_iterator iter =
		m_mapDimNameSize.find(strDimName);

	if (iter != m_mapDimNameSize.end()) {
		return iter->second;
	} else {
		return (-1);
	}
}

///////////////////////////////////////////////////////////////////////////////

void FileListObject::SortTimeArray() {

	if (m_vecTimes.size() != m_mapTimeToIndex.size()) {
		_EXCEPTIONT("vecTimes / mapTimeToIndex mismatch");
	}

	// Check if the array needs sorting, and map from old indices to new
	bool fSorted = true;
	std::map<size_t, size_t>  mapTimeIxToNewTimeIx;
	std::map<Time, size_t>::iterator iterTime = m_mapTimeToIndex.begin();
	for (size_t i = 0; iterTime != m_mapTimeToIndex.end(); iterTime++, i++) {
		if ((fSorted) && (i != iterTime->second)) {
			fSorted = false;
		}
		mapTimeIxToNewTimeIx.insert(
			std::pair<size_t, size_t>(iterTime->second, i));
	}

	if (fSorted) {
		return;
	}

	// Sort m_vecTimes and m_mapTimeToIndex
	iterTime = m_mapTimeToIndex.begin();
	for (size_t i = 0; iterTime != m_mapTimeToIndex.end(); iterTime++, i++) {
		m_vecTimes[i] = iterTime->first;
		iterTime->second = i;
	}

	// Rebuild VariableInfo VariableTimeFileMap with new time indices
	for (size_t i = 0; i < m_vecVariableInfo.size(); i++) {	
		VariableTimeFileMap mapTimeFileBak = m_vecVariableInfo[i]->m_mapTimeFile;
		m_vecVariableInfo[i]->m_mapTimeFile.clear();

		VariableTimeFileMap::const_iterator iterFileMap = mapTimeFileBak.begin();
		for (; iterFileMap != mapTimeFileBak.end(); iterFileMap++) {
			if (iterFileMap->first == InvalidTimeIx) {
				m_vecVariableInfo[i]->m_mapTimeFile.insert(
					VariableTimeFileMap::value_type(
						InvalidTimeIx,
						iterFileMap->second));
			} else {
				m_vecVariableInfo[i]->m_mapTimeFile.insert(
					VariableTimeFileMap::value_type(
						mapTimeIxToNewTimeIx[iterFileMap->first],
						iterFileMap->second));
			}
		}
	}

	// Get the VariableRegistry
	if (m_pobjRecapConfig != NULL) {
		VariableRegistry & varreg =
			m_pobjRecapConfig->GetVariableRegistry();

		varreg.UpdateTimeIndices(mapTimeIxToNewTimeIx);
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::IndexVariableData(
	size_t sFileIxBegin,
	size_t sFileIxEnd
) {
	// Open all files
	if (sFileIxBegin == InvalidFileIx) {
		sFileIxBegin = 0;
	}
	if (sFileIxEnd == InvalidFileIx) {
		sFileIxEnd = m_vecFilenames.size();
	}

	for (size_t f = sFileIxBegin; f < sFileIxEnd; f++) {

		// Open the NetCDF file
		std::string strFullFilename = m_strBaseDir + m_vecFilenames[f];
		NcFile ncFile(strFullFilename.c_str(), NcFile::ReadOnly);
		if (!ncFile.is_valid()) {
			return std::string("Unable to open data file \"")
				+ strFullFilename + std::string("\" for reading");
		}

		printf("Indexing %s\n", strFullFilename.c_str());

		// time indices stored in this file
		std::vector<size_t> vecFileTimeIndices;

		// Find the time variable, if it exists
		NcVar * varTime = ncFile.get_var(m_strRecordDimName.c_str());
		if (varTime != NULL) {
			if (varTime->num_dims() != 1) {
				return std::string("\"")
					+ m_strRecordDimName
					+ std::string("\" variable must contain exactly one dimension in \"")
					+ m_vecFilenames[f] + std::string("\"");
			}
			if ((varTime->type() != ncInt) && (varTime->type() != ncDouble)) {
				return std::string("\"")
					+ m_strRecordDimName
					+ std::string("\" variable must be ncInt or ncDouble in \"")
					+ m_vecFilenames[f] + std::string("\"");
			}

			Time::CalendarType timecal;

			NcDim * dimTime = varTime->get_dim(0);
			if (dimTime == NULL) {
				_EXCEPTION1("Malformed NetCDF file \"%s\"",
					m_vecFilenames[f].c_str());
			}

			// Get calendar
			NcAtt * attTimeCalendar = varTime->get_att("calendar");
			if (attTimeCalendar == NULL) {
				timecal = Time::CalendarStandard;
			} else {
				std::string strTimeCalendar = attTimeCalendar->as_string(0);
				timecal = Time::CalendarTypeFromString(strTimeCalendar);
				if (timecal == Time::CalendarUnknown) {
					return std::string("Unknown calendar \"") + strTimeCalendar
						+ std::string("\" in \"")
						+ m_vecFilenames[f] + std::string("\"");
				}
			}

			// Get units attribute
			NcAtt * attTimeUnits = varTime->get_att("units");
			if (attTimeUnits != NULL) {
				m_strTimeUnits = attTimeUnits->as_string(0);
			}
			if (m_strTimeUnits == "") {
				return std::string("Unknown units for \"")
					+ m_strRecordDimName
					+ std::string("\" in \"")
					+ m_vecFilenames[f] + std::string("\"");
			}

			// Add Times to master array and store corresponding indices
			// in vecFileTimeIndices.
			DataArray1D<int> nTimes(dimTime->size());
			if (varTime->type() == ncInt) {
				varTime->set_cur((long)0);
				varTime->get(&(nTimes[0]), dimTime->size());
			}

			DataArray1D<double> dTimes(dimTime->size());
			if (varTime->type() == ncDouble) {
				varTime->set_cur((long)0);
				varTime->get(&(dTimes[0]), dimTime->size());
			}

			for (int t = 0; t < dimTime->size(); t++) {
				Time time(timecal);
				if (m_strTimeUnits != "") {
					if (varTime->type() == ncInt) {
						time.FromCFCompliantUnitsOffsetInt(
							m_strTimeUnits,
							nTimes[t]);

					} else if (varTime->type() == ncDouble) {
						time.FromCFCompliantUnitsOffsetDouble(
							m_strTimeUnits,
							dTimes[t]);
					}
				}

				std::map<Time, size_t>::const_iterator iterTime =
					m_mapTimeToIndex.find(time);

				if (iterTime == m_mapTimeToIndex.end()) {
					size_t sNewIndex = m_vecTimes.size();
					m_vecTimes.push_back(time);
					vecFileTimeIndices.push_back(sNewIndex);
					m_mapTimeToIndex.insert(
						std::pair<Time, size_t>(time, sNewIndex));
				} else {
					vecFileTimeIndices.push_back(iterTime->second);
				}
			}
		}

		printf("File contains %lu times\n", vecFileTimeIndices.size());

		// Loop over all Dimensions
		const int nDims = ncFile.num_dims();
		for (int d = 0; d < nDims; d++) {
			NcDim * dim = ncFile.get_dim(d);
			std::string strDimName(dim->name());
			std::map<std::string, long>::iterator iterDim =
				m_mapDimNameSize.find(strDimName);

			if (iterDim == m_mapDimNameSize.end()) {
				m_mapDimNameSize.insert(
					std::pair<std::string, long>(
						strDimName, dim->size()));

			} else if (iterDim->second != dim->size()) {
				iterDim->second = InconsistentDimensionSizes;
			}
		}

		// Loop over all Variables
		const int nVariables = ncFile.num_vars();
		for (int v = 0; v < nVariables; v++) {
			NcVar * var = ncFile.get_var(v);
			if (var == NULL) {
				_EXCEPTION1("Malformed NetCDF file \"%s\"",
					m_vecFilenames[f].c_str());
			}

			std::string strVariableName = var->name();

			if (strVariableName == m_strRecordDimName) {
				continue;
			}

			//printf("Variable %s\n", strVariableName.c_str());

			// Find the corresponding VariableInfo structure
			std::map<std::string, size_t>::const_iterator iter =
				m_mapVariableNameToIndex.find(strVariableName);

			size_t sVarIndex;

			if (iter == m_mapVariableNameToIndex.end()) {
				sVarIndex = m_vecVariableInfo.size();

				m_vecVariableInfo.push_back(
					new VariableInfo(strVariableName));

				iter = m_mapVariableNameToIndex.insert(
					std::pair<std::string, size_t>(
						strVariableName, sVarIndex)).first;

			} else {
				sVarIndex = iter->second;
			}

			VariableInfo & info = *(m_vecVariableInfo[sVarIndex]);

			// Get units, if available
			NcAtt * attUnits = var->get_att("units");
			if (attUnits != NULL) {
				std::string strUnits = attUnits->as_string(0);
				if (strUnits != "") {
					if ((info.m_strUnits != "") && (strUnits != info.m_strUnits)) {
						return std::string("ERROR: Variable \"") + strVariableName
							+ std::string("\" has inconsistent units across files");
					}
					if (info.m_strUnits == "") {
						info.m_strUnits = strUnits;
					}
				}
			}

			// Load dimension information
			const int nDims = var->num_dims();
/*
			if (info.m_vecDimSizes.size() != 0) {
				if (info.m_vecDimSizes.size() != nDims) {
					return std::string("Variable \"") + strVariableName
						+ std::string("\" has inconsistent dimensionality across files");
				}
			}
*/
			info.m_vecDimNames.resize(nDims);
			info.m_vecDimSizes.resize(nDims);
			for (int d = 0; d < nDims; d++) {
				info.m_vecDimNames[d] = var->get_dim(d)->name();

				if (info.m_vecDimNames[d] == m_strRecordDimName) {
					if (info.m_iTimeDimIx == (-1)) {
						info.m_iTimeDimIx = d;
					} else if (info.m_iTimeDimIx != d) {
						return std::string("ERROR: Variable \"") + strVariableName
							+ std::string("\" has inconsistent \"time\" dimension across files");
					}
					info.m_vecDimSizes[d] = (-1);
				} else {
					info.m_vecDimSizes[d] = var->get_dim(d)->size();
				}
			}

			// No time information on this Variable
			if (info.m_iTimeDimIx == (-1)) {
				if (info.m_mapTimeFile.size() == 0) {
					info.m_mapTimeFile.insert(
						std::pair<size_t, LocalFileTimePair>(
							InvalidTimeIx,
							LocalFileTimePair(f, 0)));

				} else if (info.m_mapTimeFile.size() == 1) {
					VariableTimeFileMap::const_iterator iterTimeFile =
						info.m_mapTimeFile.begin();

					if (iterTimeFile->first != InvalidTimeIx) {
						return std::string("Variable \"") + strVariableName
							+ std::string("\" has inconsistent \"time\" dimension across files");
					}

				} else {
					return std::string("Variable \"") + strVariableName
						+ std::string("\" has inconsistent \"time\" dimension across files");
				}

			// Add file and time indices to VariableInfo
			} else {
				for (int t = 0; t < vecFileTimeIndices.size(); t++) {
					VariableTimeFileMap::const_iterator iterTimeFile =
						info.m_mapTimeFile.find(vecFileTimeIndices[t]);

					if (iterTimeFile == info.m_mapTimeFile.end()) {
						info.m_mapTimeFile.insert(
							std::pair<size_t, LocalFileTimePair>(
								vecFileTimeIndices[t],
								LocalFileTimePair(f, t)));

					} else {
						return std::string("Variable \"") + strVariableName
							+ std::string("\" has repeated time across files:\n")
							+ std::string("Time: ") + m_vecTimes[vecFileTimeIndices[t]].ToString() + std::string("\n")
							+ std::string("File1: ") + m_vecFilenames[iterTimeFile->second.first] + std::string("\n")
							+ std::string("File2: ") + m_vecFilenames[f];
					}
				}
			}
		}
	}

	// Sort the Time array
	SortTimeArray();

	// Update the time dimension size
	for (int v = 0; v < m_vecVariableInfo.size(); v++) {
		int iTimeDimIx = m_vecVariableInfo[v]->m_iTimeDimIx;
		if (iTimeDimIx != (-1)) {
			if (m_vecVariableInfo[v]->m_vecDimSizes.size() < iTimeDimIx) {
				_EXCEPTIONT("Logic error");
			}
			m_vecVariableInfo[v]->m_vecDimSizes[iTimeDimIx] =
				m_vecVariableInfo[v]->m_mapTimeFile.size();
		}
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::OutputTimeVariableIndexCSV(
	const std::string & strCSVOutputFilename
) {
#if defined(HYPERION_MPIOMP)
	// Only output on root thread
	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	if (nRank != 0) {
		return std::string("");
	}
#endif

	if (m_vecTimes.size() != m_mapTimeToIndex.size()) {
		_EXCEPTIONT("vecTimes / mapTimeToIndex mismatch");
	}

	std::vector< std::pair<size_t,int> > iTimeVariableIndex;
	iTimeVariableIndex.resize(m_vecVariableInfo.size());

	// Open output file
	std::ofstream ofOutput(strCSVOutputFilename.c_str());
	if (!ofOutput.is_open()) {
		return std::string("Unable to open output file \"") + strCSVOutputFilename + "\"";
	}

	// Output variables across header
	ofOutput << "time";
	for (int v = 0; v < m_vecVariableInfo.size(); v++) {
		ofOutput << "," << m_vecVariableInfo[v]->m_strVariableName;
	}
	ofOutput << std::endl;

	// Output variables with no time dimension
	ofOutput << "NONE";
	for (size_t v = 0; v < m_vecVariableInfo.size(); v++) {
		if (m_vecVariableInfo[v]->m_iTimeDimIx == (-1)) {
			ofOutput << ",X";
		} else {
			ofOutput << ",";
		}
	}
	ofOutput << std::endl;

	// Output variables with time dimension
	for (size_t t = 0; t < m_vecTimes.size(); t++) {
		ofOutput << m_vecTimes[t].ToString();

		for (size_t v = 0; v < m_vecVariableInfo.size(); v++) {

			VariableTimeFileMap::const_iterator iterTimeFile =
				m_vecVariableInfo[v]->m_mapTimeFile.find(t);

			if (iterTimeFile == m_vecVariableInfo[v]->m_mapTimeFile.end()) {
				ofOutput << ",";
			} else {
				ofOutput << "," << iterTimeFile->second.first
					<< ":" << iterTimeFile->second.second;
			}
		}
		ofOutput << std::endl;
	}

	// Output file names
	ofOutput << std::endl << std::endl;

	ofOutput << "file_ix,filename" << std::endl;
	for (size_t f = 0; f < m_vecFilenames.size(); f++) {
		ofOutput << f << ",\"" << m_strBaseDir << m_vecFilenames[f] << "\"" << std::endl;
	}

	return ("");
}

///////////////////////////////////////////////////////////////////////////////

