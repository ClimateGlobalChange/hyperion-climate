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

#include <sys/types.h>
#include <dirent.h>

///////////////////////////////////////////////////////////////////////////////

bool FileListObject::PopulateFromSearchString(
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
		strDir = ".";
	}

	// Open the directory
	DIR * pDir = opendir(strDir.c_str());
	if (pDir == NULL) {
		return false;
	}

	// Search all files in the directory for match to search string
	int iFile = 0;
	struct dirent * pDirent;
	while ((pDirent = readdir(pDir)) != NULL) {
		std::string strFilename = pDirent->d_name;
		if (STLStringHelper::WildcardMatch(
				strFileSearchString.c_str(),
				strFilename.c_str())
		) {
			std::string strChild = m_strName + "._" + std::to_string(iFile);
			std::string strFullFilename = strDir + strFilename;
			//printf("%s %s\n", strChild.c_str(), strFilename.c_str());

			// File found, insert as a member of the ListObject
			bool fSuccess =
				objreg.Assign(
					strChild,
					new StringObject(
						strChild,
						strFullFilename));

			if (!fSuccess) {
				_EXCEPTIONT("Failed to register Object");
			}

			iFile++;
			//printf ("[%s]\n", pDirent->d_name);
		}
	}
	closedir(pDir);

	return true;
}

///////////////////////////////////////////////////////////////////////////////

