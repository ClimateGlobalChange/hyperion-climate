///////////////////////////////////////////////////////////////////////////////
///
///	\file    VariableLookupObject.cpp
///	\author  Paul Ullrich
///	\version March 14, 2017
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

#include "VariableLookupObject.h"

#include <fstream>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

std::string VariableLookupObject::PopulateFromFile(
	const std::string & strFilename
) {

	// Load input file
	std::ifstream ifScript(strFilename);
	if (!ifScript.is_open()) {
		return std::string("Unable to open input file \"") + strFilename + "\"";
	}

	// Parse the input file
	int iLine = 0;
	for (std::string strLine; std::getline(ifScript, strLine); ) {

		iLine++;

		int iPos = 0;

		// Loop to find first non-whitespace character
		for (; iPos < strLine.length(); iPos++) {
			if ((strLine[iPos] != ' ') ||
			    (strLine[iPos] != '\t')
			) {
				break;
			}
		}
		if (iPos == strLine.length()) {
			continue;
		}
		if (strLine[iPos] == '#') {
			continue;
		}

		// Parse mode
		enum ParseMode {
			ParseMode_WS,
			ParseMode_CFVariable,
			ParseMode_TargetName,
			ParseMode_Units,
			ParseMode_Done
		};

		ParseMode mode = ParseMode_WS;
		ParseMode modeNext = ParseMode_CFVariable;

		// Begin parsing
		std::string strCFVariable;
		LookupInfo info;

		for (;;) {
			if (iPos == strLine.length()) {
				break;
			}
			if (strLine[iPos] == '#') {
				break;
			}

			// Ignore initial whitespace
			if (mode == ParseMode_WS) {
				if ((strLine[iPos] == ' ') ||
				    (strLine[iPos] == '\t')
				) {
					iPos++;
				} else {
					mode = modeNext;
				}

			// CF variable name
			} else if (mode == ParseMode_CFVariable) {
				if ((strLine[iPos] == ' ') ||
				    (strLine[iPos] == '\t') ||
					(strLine[iPos] == ',')
				) {
					mode = ParseMode_WS;
					modeNext = ParseMode_TargetName;
					iPos++;

				} else {
					strCFVariable += strLine[iPos];
					iPos++;
				}

			// Target name
			} else if (mode == ParseMode_TargetName) {
				if ((strLine[iPos] == ' ') ||
				    (strLine[iPos] == '\t') ||
					(strLine[iPos] == ',')
				) {
					mode = ParseMode_WS;
					modeNext = ParseMode_Units;
					iPos++;

				} else {
					info.strTargetName += strLine[iPos];
					iPos++;
				}

			// Units
			} else if (mode == ParseMode_Units) {
				if ((strLine[iPos] == ' ') ||
				    (strLine[iPos] == '\t') ||
					(strLine[iPos] == ',')
				) {
					mode = ParseMode_WS;
					modeNext = ParseMode_Done;
					iPos++;

				} else {
					info.strUnits += strLine[iPos];
					iPos++;
				}

			// Done
			} else if (mode == ParseMode_Done) {
				return std::string("Too many entries in lookup table on line ")
					+ std::to_string(iLine);
			}
		}

		// Create lookup table entry
		if (strCFVariable != "") {
			LookupTableMap::const_iterator iter =
				m_mapLookupTable.find(strCFVariable);

			if (iter != m_mapLookupTable.end()) {
				return std::string("Repeated lookup table entry on line ")
					+ std::to_string(iLine);
			}

			m_mapLookupTable.insert(
				LookupTableMapPair(strCFVariable, info));
		}
	}

	return ("");
}

///////////////////////////////////////////////////////////////////////////////

