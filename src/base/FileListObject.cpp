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
#include "STLStringHelper.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "netcdfcpp.h"

#include <sys/types.h>
#include <dirent.h>
#include <fstream>

///////////////////////////////////////////////////////////////////////////////

const size_t FileListObject::InvalidTimeIx = (-1);

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::PopulateFromSearchString(
	const std::string & strSearchString,
	ObjectRegistry & objreg
) {
	// File the directory in the search string
	std::string strDir;
	std::string strFileSearchString;
	for (int i = strSearchString.length(); i >= 0; i--) {
		if (strSearchString[i] == '/') {
			strDir = strSearchString.substr(0,i);
			strFileSearchString =
				strSearchString.substr(i+1, std::string::npos);
			break;
		}
	}
	if ((strDir == "") && (strFileSearchString == "")) {
		strFileSearchString = strSearchString;
		strDir = "./";
	}

	// Open the directory
	DIR * pDir = opendir(strDir.c_str());
	if (pDir == NULL) {
		return std::string("Unable to open directory \"")
			+ strDir + std::string("\"");
	}

	// Search all files in the directory for match to search string
	//int iFile = 0;
	struct dirent * pDirent;
	while ((pDirent = readdir(pDir)) != NULL) {
		std::string strFilename = pDirent->d_name;
		if (STLStringHelper::WildcardMatch(
				strFileSearchString.c_str(),
				strFilename.c_str())
		) {
			//std::string strChild = m_strName + "._" + std::to_string(iFile);
			std::string strFullFilename = strDir + strFilename;
			//printf("%s %s\n", strChild.c_str(), strFilename.c_str());

			// File found, insert into list of filenames
			m_vecFilenames.push_back(strFullFilename);
			//bool fSuccess =
			//	objreg.Assign(
			//		strChild,
			//		new StringObject(
			//			strChild,
			//			strFullFilename));

			//if (!fSuccess) {
			//	_EXCEPTIONT("Failed to register Object");
			//}

			//iFile++;
			//printf ("[%s]\n", pDirent->d_name);
		}
	}
	closedir(pDir);

	return IndexVariableData();;
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::IndexVariableData() {

	// A map from Variable name to m_vecVariableInfo vector index
	std::map<std::string, size_t> mapVariableNameToIndex;

	// A map from Time to m_vecTimes vector index
	std::map<Time, size_t> mapTimeToIndex;

	// Open all files
	for (size_t f = 0; f < m_vecFilenames.size(); f++) {

		// Open the NetCDF file
		NcFile ncFile(m_vecFilenames[f].c_str(), NcFile::ReadOnly);
		if (!ncFile.is_valid()) {
			return std::string("Unable to open data file \"")
				+ m_vecFilenames[f] + std::string("\" for reading");
		}

		printf("Indexing %s\n", m_vecFilenames[f].c_str());

		// time indices stored in this file
		std::vector<size_t> vecFileTimeIndices;

		// Find the time variable, if it exists
		NcVar * varTime = ncFile.get_var("time");
		if (varTime != NULL) {

			if (varTime->num_dims() != 1) {
				return std::string("\"time\" variable must contain exactly one dimension in \"")
					+ m_vecFilenames[f] + std::string("\"");
			}
			if ((varTime->type() != ncInt) && (varTime->type() != ncDouble)) {
				return std::string("\"time\" variable must be ncInt or ncDouble in \"")
					+ m_vecFilenames[f] + std::string("\"");
			}

			Time::CalendarType timecal;
			std::string strTimeUnits;

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
				strTimeUnits = attTimeUnits->as_string(0);
			}
			if (strTimeUnits == "") {
				return std::string("Unknown units for \"time\" in \"")
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
				if (strTimeUnits != "") {
					if (varTime->type() == ncInt) {
						time.FromCFCompliantUnitsOffsetInt(
							strTimeUnits,
							nTimes[t]);

					} else if (varTime->type() == ncDouble) {
						time.FromCFCompliantUnitsOffsetDouble(
							strTimeUnits,
							dTimes[t]);
					}
				}

				std::map<Time, size_t>::const_iterator iterTime =
					mapTimeToIndex.find(time);

				if (iterTime == mapTimeToIndex.end()) {
					size_t sNewIndex = m_vecTimes.size();
					m_vecTimes.push_back(time);
					vecFileTimeIndices.push_back(sNewIndex);
					mapTimeToIndex.insert(
						std::pair<Time, size_t>(time, sNewIndex));
				} else {
					vecFileTimeIndices.push_back(iterTime->second);
				}
			}
		}

		printf("File contains %lu times\n", vecFileTimeIndices.size());

		// Loop over all Variables
		const int nVariables = ncFile.num_vars();
		for (int v = 0; v < nVariables; v++) {
			NcVar * var = ncFile.get_var(v);
			if (var == NULL) {
				_EXCEPTION1("Malformed NetCDF file \"%s\"",
					m_vecFilenames[f].c_str());
			}

			std::string strVariableName = var->name();

			// Ignore the "time" variable
			if (strVariableName == "time") {
				continue;
			}

			//printf("Variable %s\n", strVariableName.c_str());

			// Find the corresponding VariableInfo structure
			std::map<std::string, size_t>::const_iterator iter =
				mapVariableNameToIndex.find(strVariableName);

			size_t sVarIndex;

			if (iter == mapVariableNameToIndex.end()) {
				sVarIndex = m_vecVariableInfo.size();

				m_vecVariableInfo.push_back(
					VariableInfo(strVariableName));

				iter = mapVariableNameToIndex.insert(
					std::pair<std::string, size_t>(
						strVariableName, sVarIndex)).first;

			} else {
				sVarIndex = iter->second;
			}

			VariableInfo & info = m_vecVariableInfo[sVarIndex];

			// Get units, if available
			NcAtt * attUnits = var->get_att("units");
			if (attUnits != NULL) {
				std::string strUnits = attUnits->as_string(0);
				if (strUnits != "") {
					if ((info.m_strUnits != "") && (strUnits != info.m_strUnits)) {
						return std::string("Variable \"") + strVariableName
							+ std::string("\" has inconsistent units across files");
					}
					if (info.m_strUnits == "") {
						info.m_strUnits = strUnits;
					}
				}
			}

			// Load dimension sizes
			const int nDims = var->num_dims();
			if (info.m_vecDimSizes.size() != 0) {
				if (info.m_vecDimSizes.size() != nDims) {
					return std::string("Variable \"") + strVariableName
						+ std::string("\" has inconsistent dimensionality across files");
				}
			}

			info.m_vecDimSizes.resize(nDims);
			for (int d = 0; d < nDims; d++) {
				info.m_vecDimSizes[d] = var->get_dim(d)->size();
			}

			// Check for time dimesion
			for (int d = 0; d < nDims; d++) {
				if (0 == strcmp(var->get_dim(d)->name(), "time")) {
					if (info.m_iTimeDimIx == (-1)) {
						info.m_iTimeDimIx = d;
					} else if (info.m_iTimeDimIx != d) {
						return std::string("Variable \"") + strVariableName
							+ std::string("\" has inconsistent \"time\" dimension across files");
					}
					break;
				}
			}

			// No time information on this Variable
			if (info.m_iTimeDimIx == (-1)) {
				if (info.m_mapTimeFile.size() == 0) {
					info.m_mapTimeFile.insert(
						std::pair<size_t, VariableInfo::LocalFileTimePair>(
							InvalidTimeIx,
							VariableInfo::LocalFileTimePair(f, 0)));

				} else if (info.m_mapTimeFile.size() == 1) {
					VariableInfo::VariableTimeFileMap::const_iterator iterTimeFile =
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
					VariableInfo::VariableTimeFileMap::const_iterator iterTimeFile =
						info.m_mapTimeFile.find(vecFileTimeIndices[t]);

					if (iterTimeFile == info.m_mapTimeFile.end()) {
						info.m_mapTimeFile.insert(
							std::pair<size_t, VariableInfo::LocalFileTimePair>(
								vecFileTimeIndices[t],
								VariableInfo::LocalFileTimePair(f, t)));

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

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string FileListObject::OutputTimeVariableIndexCSV(
	const std::string & strCSVOutputFilename
) {
	std::vector< std::pair<size_t,int> > iTimeVariableIndex;
	iTimeVariableIndex.resize(m_vecVariableInfo.size());

	// Load input file
	std::ofstream ofOutput(strCSVOutputFilename.c_str());
	if (!ofOutput.is_open()) {
		return std::string("Unable to open output file \"") + strCSVOutputFilename + "\"";
	}

	// Output variables across header
	for (int v = 0; v < m_vecVariableInfo.size(); v++) {
		ofOutput << "," << m_vecVariableInfo[v].m_strVariableName;
	}
	ofOutput << std::endl;

	// Output variables with no time dimension
	ofOutput << "NONE";
	for (size_t v = 0; v < m_vecVariableInfo.size(); v++) {
		if (m_vecVariableInfo[v].m_iTimeDimIx == (-1)) {
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

			VariableInfo::VariableTimeFileMap::const_iterator iterTimeFile =
				m_vecVariableInfo[v].m_mapTimeFile.find(t);

			if (iterTimeFile == m_vecVariableInfo[v].m_mapTimeFile.end()) {
				ofOutput << ",";
			} else {
				ofOutput << "," << iterTimeFile->second.first
					<< ":" << iterTimeFile->second.second;
			}
		}
		ofOutput << std::endl;
	}

	return ("");
}

///////////////////////////////////////////////////////////////////////////////

