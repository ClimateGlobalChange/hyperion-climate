///////////////////////////////////////////////////////////////////////////////
///
///	\file    recap.cpp
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

#include "netcdfcpp.h"

#include "Announce.h"
#include "Object.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(HYPERION_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Output usage information
	if (argc < 2) {
		Announce("recap version 0.1, March 8th, 2017");
		Announce("Usage: %s <input file> [assignment list]", argv[0]);
		return (0);
	}

	// Load input file
	std::ifstream ifScript(argv[1]);
	if (!ifScript.is_open()) {
		Announce("ERROR: Unable to open input file \"%s\"", argv[1]);
		return (-1);
	}

	// Load command-line variables


	// Begin parsing
	int iLine = 0;

	for (std::string strLine; std::getline(ifScript, strLine); ) {

		// Increment line number 
		iLine++;

		// Parser state
		enum ParserState {
			ParserState_WS,
			ParserState_Token,
			ParserState_String,
			ParserState_Op
		};

		// Parsed command line
		std::vector<std::string> vecCommandLine;
		std::vector<ParserState> vecCommandLineType;

		ParserState parse_state = ParserState_WS;

		std::stack<char> stParentheses;

		// Begin parsing
		bool fError = false;
		int iPos = 0;
		for (int i = 0; i < strLine.length();) {

			// Whitespace
			if (parse_state == ParserState_WS) {
				if ((strLine[i] == ' ') || (strLine[i] == '\t')) {
					i++;
					iPos = i;
					continue;

				} else if (strLine[i] == '#') {
					break;

				} else if (((strLine[i] >= 'A') && (strLine[i] <= 'Z')) ||
				    ((strLine[i] >= 'a') && (strLine[i] <= 'z')) ||
					(strLine[i] == '_')
				) {
					parse_state = ParserState_Token;
					continue;

				} else if (
					(strLine[i] == '=') ||
					(strLine[i] == '(') ||
					(strLine[i] == ')') ||
					(strLine[i] == '[') ||
					(strLine[i] == ']') ||
					(strLine[i] == ',') ||
					(strLine[i] == '\"') ||
					(strLine[i] == '.')
				) {
					if (vecCommandLine.size() == 0) {
						Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
					parse_state = ParserState_Op;
					continue;

				} else {
					Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
					fError = true;
					break;
				}

			// Token
			} else if (parse_state == ParserState_Token) {
				if (((strLine[i] >= 'A') && (strLine[i] <= 'Z')) ||
				    ((strLine[i] >= 'a') && (strLine[i] <= 'z')) ||
					(strLine[i] == '_')
				) {
					i++;
					continue;

				} else {
					if (iPos != i) {
						vecCommandLine.push_back(
							strLine.substr(iPos, i-iPos));
						vecCommandLineType.push_back(
							ParserState_Token);
						iPos = i;
					}

					if (
						(strLine[i] == '=') ||
						(strLine[i] == '(') ||
						(strLine[i] == ')') ||
						(strLine[i] == '[') ||
						(strLine[i] == ']') ||
						(strLine[i] == ',') ||
						(strLine[i] == '\"') ||
						(strLine[i] == '.')
					) {
						parse_state = ParserState_Op;
						continue;

					} else if (strLine[i] == ' ') {
						parse_state = ParserState_WS;
						continue;
					}
				}

			// String
			} else if (parse_state == ParserState_String) {
				if (strLine[i] == '\"') {
					vecCommandLine.push_back(
						strLine.substr(iPos, i-iPos));
					vecCommandLineType.push_back(
						ParserState_String);
					parse_state = ParserState_Op;
					continue;
				}
				i++;

			// Op
			} else if (parse_state == ParserState_Op) {
				if (strLine[i] == ' ') {
					parse_state = ParserState_WS;
					continue;
				}

				if (strLine[i] == '(') {
					stParentheses.push('(');
					vecCommandLine.push_back("(");
					vecCommandLineType.push_back(ParserState_Op);
					parse_state = ParserState_Token;

				} else if (strLine[i] == '[') {
					stParentheses.push('[');
					vecCommandLine.push_back("[");
					vecCommandLineType.push_back(ParserState_Op);
					parse_state = ParserState_Token;

				} else if (strLine[i] == '\"') {
					if ((stParentheses.size() != 0) &&
					    (stParentheses.top() == '\"')
					) {
						stParentheses.pop();

					} else {
						stParentheses.push('\"');
						parse_state = ParserState_String;
					}

				} else if (strLine[i] == ')') {
					if (stParentheses.size() == 0) {
						Announce("ERROR: Unbalanced parenthesis on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
					if (stParentheses.top() != '(') {
						Announce("ERROR: Unbalanced parenthesis on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
					stParentheses.pop();
					vecCommandLine.push_back(")");
					vecCommandLineType.push_back(ParserState_Op);

				} else if (strLine[i] == ']') {
					if (stParentheses.size() == 0) {
						Announce("ERROR: Unbalanced bracket on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
					if (stParentheses.top() != '[') {
						Announce("ERROR: Unbalanced bracket on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
					stParentheses.pop();
					vecCommandLine.push_back("]");
					vecCommandLineType.push_back(ParserState_Op);

				} else if (strLine[i] == '=') {
					vecCommandLine.push_back("=");
					vecCommandLineType.push_back(ParserState_Op);
					parse_state = ParserState_Token;

				} else if (strLine[i] == ',') {
					vecCommandLine.push_back(",");
					vecCommandLineType.push_back(ParserState_Op);
					parse_state = ParserState_Token;

				} else if (strLine[i] == '.') {
					vecCommandLine.push_back(".");
					vecCommandLineType.push_back(ParserState_Op);
					parse_state = ParserState_Token;

				} else {
					Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
					fError = true;
					break;
				}

				i++;
				iPos = i;
			}
		}
		if (vecCommandLine.size() != vecCommandLineType.size()) {
			_EXCEPTIONT("CommandLine/CommandLineType mismatch");
		}

		if (vecCommandLine.size() == 0) {
			continue;
		}

		if (parse_state == ParserState_String) {
			Announce("ERROR: Unbalanced quotation marks on line %i", iLine);
			fError = true;
		}

		if (stParentheses.size() != 0) {
			Announce("ERROR: Unbalanced parenthesis or bracket on line %i", iLine);
			std::cout << stParentheses.top() << std::endl;
			fError = true;
		}

		if (fError) {
			return (-1);
		}

		for (int i = 0; i < vecCommandLine.size(); i++) {
			printf("%s:", vecCommandLine[i].c_str());
		}
		printf("\n");
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(HYPERION_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}

///////////////////////////////////////////////////////////////////////////////

