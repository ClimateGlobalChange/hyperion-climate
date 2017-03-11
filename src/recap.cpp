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
#include "Variable.h"
#include "GridObject.h"
#include "FileListObject.h"

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


	// Create Object and Variable registry
	ObjectRegistry objreg;

	VariableRegistry varreg;

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
			ParserState_Number,
			ParserState_Op
		};

		// Parsed command line
		std::vector<std::string> vecCommandLine;
		std::vector<ObjectType> vecCommandLineType;

		ObjectType objtype;
		ParserState parse_state = ParserState_WS;

		std::stack<char> stParentheses;

		// Begin parsing
		bool fError = false;
		int iPos = 0;
		for (int i = 0; i <= strLine.length();) {

			// Whitespace
			if (parse_state == ParserState_WS) {
				if ((strLine[i] == ' ') || (strLine[i] == '\t')) {
					i++;
					iPos = i;
					continue;

				} else if (((strLine[i] >= 'A') && (strLine[i] <= 'Z')) ||
				    ((strLine[i] >= 'a') && (strLine[i] <= 'z')) ||
					(strLine[i] == '_')
				) {
					parse_state = ParserState_Token;
					continue;

				} else if (((strLine[i] >= '0') && (strLine[i] <= '9')) ||
					(strLine[i] == '.')
				) {
					objtype = ObjectType_Integer;
					parse_state = ParserState_Number;
					continue;

				} else if (
					(strLine[i] == '=') ||
					(strLine[i] == '(') ||
					(strLine[i] == ')') ||
					(strLine[i] == '[') ||
					(strLine[i] == ']') ||
					(strLine[i] == ',') ||
					(strLine[i] == '\"')
				) {
					if (vecCommandLine.size() == 0) {
						Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
					parse_state = ParserState_Op;
					continue;

				} else if (
					(strLine[i] == '#') ||
					(strLine[i] == '\0')
				) {
					break;

				} else {
					Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
					fError = true;
					break;
				}

			// Token
			} else if (parse_state == ParserState_Token) {
				if (((strLine[i] >= 'A') && (strLine[i] <= 'Z')) ||
				    ((strLine[i] >= 'a') && (strLine[i] <= 'z')) ||
					(strLine[i] == '_') ||
					(strLine[i] == '.')
				) {
					i++;
					continue;

				} else if (
					(i != iPos) &&
					((strLine[i] >= '0') && (strLine[i] <= '9'))
				) {
					i++;
					continue;

				} else {
					if (iPos != i) {
						vecCommandLine.push_back(
							strLine.substr(iPos, i-iPos));
						vecCommandLineType.push_back(
							ObjectType_Token);
						iPos = i;
					}

					if ((strLine[i] == '#') ||
					    (strLine[i] == '\0')
					) {
						break;
					}

					if (strLine[i] == '\"') {
						Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
						fError = true;
						break;

					} else if (
						(strLine[i] == '=') ||
						(strLine[i] == '(') ||
						(strLine[i] == ')') ||
						(strLine[i] == '[') ||
						(strLine[i] == ']') ||
						(strLine[i] == ',')
					) {
						parse_state = ParserState_Op;
						continue;

					} else if (strLine[i] == ' ') {
						parse_state = ParserState_WS;
						continue;

					} else {
						Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
				}

			// String
			} else if (parse_state == ParserState_String) {

				if (strLine[i] == '\0') {
					Announce("ERROR: Incomplete string on line %i", iLine);
					fError = true;
					break;
				}

				if (strLine[i] == '\"') {
					vecCommandLine.push_back(
						strLine.substr(iPos, i-iPos));
					vecCommandLineType.push_back(
						ObjectType_String);
					parse_state = ParserState_Op;
					continue;
				}
				i++;

			// Number
			} else if (parse_state == ParserState_Number) {
				if (((strLine[i] >= '0') && (strLine[i] <= '9')) ||
					(strLine[i] == '.')
				) {
					if ((strLine[i] == '.') && (objtype == ObjectType_Integer)) {
						objtype = ObjectType_FloatingPoint;
					}
					i++;
					continue;

				} else {
					if (iPos != i) {
						vecCommandLine.push_back(
							strLine.substr(iPos, i-iPos));
						vecCommandLineType.push_back(
							objtype);
						iPos = i;
					}

					if ((strLine[i] == '#') ||
					    (strLine[i] == '\0')
					) {
						break;
					}

					if ((strLine[i] == '=') ||
						(strLine[i] == '(') ||
						(strLine[i] == '[') ||
						(strLine[i] == '\"')
					) {
						Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
						fError = true;
						break;

					} else if (
						(strLine[i] == ',') ||
						(strLine[i] == ')') ||
						(strLine[i] == ']')
					)  {
						parse_state = ParserState_Op;
						continue;

					} else if (strLine[i] == ' ') {
						parse_state = ParserState_WS;
						continue;

					} else {
						Announce("ERROR: Unexpected character on line %i (%i)", iLine, i);
						fError = true;
						break;
					}
				}

			// Op
			} else if (parse_state == ParserState_Op) {
				if (strLine[i] == ' ') {
					parse_state = ParserState_WS;
					continue;
				}

				if ((strLine[i] == '#') ||
				    (strLine[i] == '\0')
				) {
					break;

				} else if (strLine[i] == '(') {
					stParentheses.push('(');
					vecCommandLine.push_back("(");
					vecCommandLineType.push_back(ObjectType_Op);
					parse_state = ParserState_WS;

				} else if (strLine[i] == '[') {
					stParentheses.push('[');
					vecCommandLine.push_back("[");
					vecCommandLineType.push_back(ObjectType_Op);
					parse_state = ParserState_WS;

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
					vecCommandLineType.push_back(ObjectType_Op);
					parse_state = ParserState_WS;

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
					vecCommandLineType.push_back(ObjectType_Op);
					parse_state = ParserState_WS;

				} else if (strLine[i] == '=') {
					vecCommandLine.push_back("=");
					vecCommandLineType.push_back(ObjectType_Op);
					parse_state = ParserState_WS;

				} else if (strLine[i] == ',') {
					vecCommandLine.push_back(",");
					vecCommandLineType.push_back(ObjectType_Op);
					parse_state = ParserState_WS;

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

		if (fError) {
			return (-1);
		}

		if (parse_state == ParserState_String) {
			Announce("ERROR: Unbalanced quotation marks on line %i", iLine);
			return (-1);

		} else if (stParentheses.size() != 0) {
			Announce("ERROR: Unbalanced parenthesis or bracket on line %i", iLine);
			return (-1);
		}
/*
		for (int i = 0; i < vecCommandLine.size(); i++) {
			printf("%s:", vecCommandLine[i].c_str());
		}
		printf("\n");
*/
		// Blank command line (do nothing)
		if (vecCommandLine.size() == 0) {
			continue;

		// No valid command lines with one token
		} else if (vecCommandLine.size() == 1) {
			Announce("ERROR: Syntax error on line %i", iLine);
			return (-1);

		// Assignment or Evaluate op
		} else {
			int iAssignmentOp = (-1);
			int iEvaluateOp = (-1);
			for (int i = 0; i < vecCommandLine.size(); i++) {
				if ((vecCommandLineType[i] == ObjectType_Op) &&
				    (vecCommandLine[i] == "=")
				) {
					if (iEvaluateOp != (-1)) {
						Announce("ERROR: Syntax error on line %i", iLine);
						return (-1);
					}
					if (i != 1) {
						Announce("ERROR: Syntax error on line %i", iLine);
						return (-1);
					}
					iAssignmentOp = i;
				}
				if ((vecCommandLineType[i] == ObjectType_Op) &&
				    (vecCommandLine[i] == "(")
				) {
					iEvaluateOp = i;
				}
			}

			// Assignment but insufficient tokens in line
			if ((iAssignmentOp != (-1)) && (vecCommandLine.size() < 3)) {
				Announce("ERROR: Syntax error in assignment on line %i", iLine);
				return (-1);
			}

			// No assignment or function evaluation: syntax error
			if ((iAssignmentOp == (-1)) && (iEvaluateOp == (-1))) {
				Announce("ERROR: Syntax error on line %i", iLine);
				return (-1);
			}

			// Assignment of primitive type to LHS
			if ((iAssignmentOp != (-1)) && (iEvaluateOp == (-1))) {
				bool fSuccess = true;

				// List type on RHS
				if ((vecCommandLineType[2] == ObjectType_Op) &&
					(vecCommandLine[2] == "[")
				) {
					printf("LIST: %s\n", vecCommandLine[0].c_str());

					ListObject * pobjList = new ListObject(vecCommandLine[0]);

					objreg.Assign(
						vecCommandLine[0],
						pobjList);

					int iListEntry = 0;
					for (int i = 3; i < vecCommandLine.size(); i++) {

						if (vecCommandLine[i] == "]") {
							break;
						}
						if (vecCommandLine[i] == ",") {
							if (iListEntry == 0) {
								Announce("ERROR: Malformed list after entry 0 on line %i",
									iLine);
								return (-1);

							} else {
								i++;
							}

						} else if (iListEntry != 0) {
							Announce("ERROR: Malformed list after entry %i on line %i",
								iListEntry, iLine);
							return (-1);
						}

						std::string strChildName =
							vecCommandLine[0] + "._" + std::to_string(iListEntry);
/*
						printf("%i %i %i %s %s\n",
							i, iListEntry,
							vecCommandLineType[i],
							strChildName.c_str(),
							vecCommandLine[i].c_str());
*/
						fSuccess =
							objreg.Create(
								vecCommandLineType[i],
								strChildName,
								vecCommandLine[i]);

						if (!fSuccess) {
							Announce("ERROR: Invalid list entry %i on line %i",
								iListEntry, iLine);
							return (-1);
						}

						pobjList->PushBack(strChildName);
						iListEntry++;
					}

				// Try to create the object as a primitive
				} else {
					fSuccess =
						objreg.Create(
							vecCommandLineType[2],
							vecCommandLine[0],
							vecCommandLine[2]);
				}

				if (!fSuccess) {
					Announce("ERROR: Invalid RHS in assignment on line %i", iLine);
					return (-1);
				}
			}

			// Assignment of function to LHS
			if ((iAssignmentOp != (-1)) && (iEvaluateOp != (-1))) {

				// parameter_list() type
				if (vecCommandLine[2] == "parameter_list") {
					//printf("PLI: %s\n", vecCommandLine[0].c_str()); 
					objreg.Assign(
						vecCommandLine[0],
						new Object(vecCommandLine[0]));

				// variable(STRING) type
				} else if (vecCommandLine[2] == "variable") {
					if (vecCommandLine.size() != 6) {
						Announce("ERROR: Variable op declaration missing"
							" operation on line %i", iLine);
						return (-1);
					}
					if (vecCommandLine[3] != "(") {
						Announce("ERROR: Syntax error on line %i", iLine);
						return (-1);
					}
					if (vecCommandLineType[4] != ObjectType_String) {
						Announce("ERROR: Invalid variable op declaration on line %i", iLine);
						return (-1);
					}
					//printf("VAR: %s %s\n", vecCommandLine[0].c_str(), vecCommandLine[4].c_str());
					objreg.Assign(
						vecCommandLine[0],
						new VariableObject(
							vecCommandLine[0], vecCommandLine[4]));

				// grid(STRING) type
				} else if (vecCommandLine[2] == "grid_file") {
					if (vecCommandLine.size() != 6) {
						Announce("ERROR: grid_file filename argument missing"
							" on line %i", iLine);
						return (-1);
					}
					if (vecCommandLine[3] != "(") {
						Announce("ERROR: Syntax error on line %i", iLine);
						return (-1);
					}
					if (vecCommandLineType[4] != ObjectType_String) {
						Announce("ERROR: Invalid grid_file declaration on line %i", iLine);
						return (-1);
					}

					printf("GRID: %s %s\n", vecCommandLine[0].c_str(), vecCommandLine[4].c_str());
					objreg.Assign(
						vecCommandLine[0],
						new GridObject(
							vecCommandLine[0], vecCommandLine[4]));

				// file_list(STRING) type
				} else if (vecCommandLine[2] == "file_list") {
					if (vecCommandLine.size() != 6) {
						Announce("ERROR: file_list search string argument missing"
							" on line %i", iLine);
						return (-1);
					}
					if (vecCommandLine[3] != "(") {
						Announce("ERROR: Syntax error on line %i", iLine);
						return (-1);
					}
					if (vecCommandLineType[4] != ObjectType_String) {
						Announce("ERROR: Invalid file_list declaration on line %i", iLine);
						return (-1);
					}

					printf("FILELIST: %s %s\n", vecCommandLine[0].c_str(), vecCommandLine[4].c_str());
					FileListObject * pObj =
						new FileListObject(
							vecCommandLine[0]);

					objreg.Assign(vecCommandLine[0], pObj);

					pObj->PopulateFromSearchString(
							vecCommandLine[4],
							objreg);

					Announce("file_list %s contains %i entries",
						vecCommandLine[0].c_str(),
						pObj->ChildrenCount());

				// Unknown function (keep as a WARNING for now)
				} else {
					Announce("WARNING: Unknown function \"%s\" on line %i",
						vecCommandLine[2].c_str(), iLine);
				}
			}
		}
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

