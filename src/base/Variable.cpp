///////////////////////////////////////////////////////////////////////////////
///
///	\file    Variable.cpp
///	\author  Paul Ullrich
///	\version November 18, 2015
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "Variable.h"

#include "RecapConfigObject.h"
#include "FileListObject.h"

#include <set>

///////////////////////////////////////////////////////////////////////////////
// VariableRegistry
///////////////////////////////////////////////////////////////////////////////

Variable * VariableRegistry::FindOrRegister(
	const std::string & strVariableName,
	const FileListObject * pobjFileList
) {
	std::map<std::string, Variable>::iterator iterVar =
		m_mapVariables.find(strVariableName);

	// Variable exists in registry
	if (iterVar != m_mapVariables.end()) {
		return &(iterVar->second);
	}

	// Check if variable is in the FileList, and add to registry if present
	const VariableInfo * pvarinfo =
		pobjFileList->GetVariableInfo(strVariableName);
	if (pvarinfo != NULL) {
		iterVar =
			m_mapVariables.insert(
				std::pair<std::string, Variable>(
					strVariableName,
					Variable(this, strVariableName, pvarinfo))).first;

		return &(iterVar->second);
	}

	// Combinations not yet supported
	return NULL;
}

///////////////////////////////////////////////////////////////////////////////

bool VariableRegistry::FindOrRegister(
	const Variable & var
) {
	std::map<std::string, Variable>::const_iterator iter =
		m_mapVariables.find(var.m_strName);

	if (iter == m_mapVariables.end()) {
		m_mapVariables.insert(
			std::pair<std::string, Variable>(var.m_strName, var));

		return true;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

Variable * VariableRegistry::Get(
	const std::string & strName
) {
	std::map<std::string, Variable>::iterator iter =
		m_mapVariables.find(strName);

	if (iter == m_mapVariables.end()) {
		return NULL;
	} else {
		return &(iter->second);
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::UnloadAllGridData() {
	std::map<std::string, Variable>::iterator iter =
		m_mapVariables.begin();
	for (; iter != m_mapVariables.end(); iter++) {
		iter->second.UnloadGridData();
	}
}

///////////////////////////////////////////////////////////////////////////////
// Variable
///////////////////////////////////////////////////////////////////////////////

Variable::Variable(
	VariableRegistry * pvarreg,
	const std::string & strName,
	const VariableInfo * pvarinfo
) :
	m_pvarreg(pvarreg),
	m_strName(strName),
	m_pvarinfo(pvarinfo),
	m_fOp(false),
	m_strUnits(),
	m_nSpecifiedDim(0),
	m_fNoTimeInNcFile(false),
	m_sTime(-2)
{
	if (pvarinfo != NULL) {
		m_strUnits = pvarinfo->m_strUnits;
	}
	memset(m_iDim, 0, MaxArguments * sizeof(int));
}

///////////////////////////////////////////////////////////////////////////////

bool Variable::operator==(
	const Variable & var
) {
	if (m_fOp != var.m_fOp) {
		return false;
	}
	if (m_strName != var.m_strName) {
		return false;
	}
	if (m_nSpecifiedDim != var.m_nSpecifiedDim) {
		return false;
	}
	for (int i = 0; i < m_nSpecifiedDim; i++) {
		if (m_iDim[i] != var.m_iDim[i]) {
			return false;
		}
	}
	for (int i = 0; i < m_varArg.size(); i++) {
		if (m_varArg[i] != var.m_varArg[i]) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

int Variable::ParseFromString(
	const std::string & strIn
) {
	// Verify initialization
	if (m_pvarreg == NULL) {
		_EXCEPTIONT("Invalid VariableRegistry");
	}

	m_fOp = false;
	m_strName = "";
	m_nSpecifiedDim = 0;
	m_varArg.clear();

	bool fDimMode = false;
	std::string strDim;

	if (strIn.length() >= 1) {
		if (strIn[0] == '_') {
			m_fOp = true;
		}
	}

	for (int n = 0; n <= strIn.length(); n++) {
		// Reading the variable name
		if (!fDimMode) {
			if (n == strIn.length()) {
				m_strName = strIn;
				return n;
			}
			if (strIn[n] == ',') {
				m_strName = strIn.substr(0, n);
				return n;
			}
			if (strIn[n] == '(') {
				m_strName = strIn.substr(0, n);
				fDimMode = true;
				continue;
			}
			if (strIn[n] == ')') {
				m_strName = strIn.substr(0, n);
				return n;
			}

		// Reading in dimensions
		} else if (!m_fOp) {
			if (m_nSpecifiedDim == 4) {
				_EXCEPTIONT("Only 4 dimensions / arguments may "
					"be specified");
			}
			if (n == strIn.length()) {
				_EXCEPTION1("Variable dimension list must be terminated"
					" with ): %s", strIn.c_str());
			}
			if (strIn[n] == ',') {
				if (strDim.length() == 0) {
					_EXCEPTIONT("Invalid dimension index in variable");
				}
				m_iDim[m_nSpecifiedDim] = atoi(strDim.c_str());
				m_nSpecifiedDim++;
				strDim = "";

			} else if (strIn[n] == ')') {
				if (strDim.length() == 0) {
					_EXCEPTIONT("Invalid dimension index in variable");
				}
				m_iDim[m_nSpecifiedDim] = atoi(strDim.c_str());
				m_nSpecifiedDim++;
				return (n+1);

			} else {
				strDim += strIn[n];
			}

		// Reading in arguments
		} else {
			if (m_nSpecifiedDim == 4) {
				_EXCEPTIONT("Only 4 dimensions / arguments may "
					"be specified");
			}
			if (n == strIn.length()) {
				_EXCEPTION1("Op argument list must be terminated"
					" with ): %s", strIn.c_str());
			}

			// No arguments
			if (strIn[n] == ')') {
				return (n+1);
			}

			// Parse arguments
			m_varArg.resize(m_nSpecifiedDim+1);

			Variable var(m_pvarreg);
			n += var.ParseFromString(strIn.substr(n));

			m_pvarreg->FindOrRegister(var);
			m_vecArg.push_back(var.m_strName);

			m_nSpecifiedDim++;

			if (strIn[n] == ')') {
				return (n+1);
			}
		}
	}

	_EXCEPTION1("Malformed variable string \"%s\"", strIn.c_str());
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::ToString() const {

	// Verify initialization
	if (m_pvarreg == NULL) {
		_EXCEPTIONT("Invalid VariableRegistry");
	}

	// Convert to string
	char szBuffer[20];
	std::string strOut = m_strName;
	if (m_nSpecifiedDim == 0) {
		return strOut;
	}
	strOut += "(";
	for (int d = 0; d < m_nSpecifiedDim; d++) {
		if (m_fOp) {
			Variable * var = m_pvarreg->Get(m_vecArg[d]);
			if (var == NULL) {
				_EXCEPTION();
			} else {
				strOut += var->ToString();
			}
		} else {
			sprintf(szBuffer, "%i", m_iDim[d]);
			strOut += szBuffer;
		}
		if (d != m_nSpecifiedDim-1) {
			strOut += ",";
		} else {
			strOut += ")";
		}
	}
	return strOut;
}

///////////////////////////////////////////////////////////////////////////////

void Variable::LoadGridData(
	RecapConfigObject * pobjConfig,
	size_t sTime
) {
	// Verify initialization
	if (m_pvarreg == NULL) {
		_EXCEPTIONT("Invalid VariableRegistry");
	}

	// Verify argument
	if (pobjConfig == NULL) {
		_EXCEPTIONT("Invalid pobjConfig argument");
	}

	// Get Mesh
	GridObject * pobjGrid = pobjConfig->GetGrid();
	if (pobjGrid == NULL) {
		_EXCEPTIONT("Invalid configuration: Missing grid");
	}

	const Mesh & mesh = pobjGrid->GetMesh();

	// Get FileList
	FileListObject * pobjFileList = pobjConfig->GetFileList();
	if (pobjFileList == NULL) {
		_EXCEPTIONT("Invalid configuration: Missing file_list");
	}

	// Check if data already loaded
	if (sTime == m_sTime) {
		if (m_data.GetRows() != mesh.sDOFCount) {
			_EXCEPTIONT("Logic error");
		}
		return;
	}

	// Allocate data
	m_data.Allocate(mesh.sDOFCount);
	m_sTime = sTime;

	// Get the data directly from a variable
	if (!m_fOp) {
		pobjFileList->LoadData_float(
			m_strName,
			sTime,
			m_data);

		return;
	}
/*
	// Evaluate the vector magnitude operator
	if (m_strName == "_VECMAG") {
		if (m_varArg.size() != 2) {
			_EXCEPTION1("_VECMAG expects two arguments: %i given",
				m_varArg.size());
		}
		Variable & varLeft = varreg.Get(m_varArg[0]);
		Variable & varRight = varreg.Get(m_varArg[1]);

		varLeft.LoadGridData(varreg, vecFiles, grid, sTime);
		varRight.LoadGridData(varreg, vecFiles, grid, sTime);

		const DataArray1D<float> & dataLeft  = varLeft.GetData();
		const DataArray1D<float> & dataRight = varRight.GetData();

		for (int i = 0; i < m_data.GetRows(); i++) {
			m_data[i] =
				sqrt(dataLeft[i] * dataLeft[i]
					+ dataRight[i] * dataRight[i]);
		}

	// Evaluate the absolute value operator
	} else if (m_strName == "_ABS") {
		if (m_varArg.size() != 1) {
			_EXCEPTION1("_ABS expects one argument: %i given",
				m_varArg.size());
		}

		m_data.Zero();

		Variable & varParam = varreg.Get(m_varArg[0]);
		varParam.LoadGridData(varreg, vecFiles, grid, sTime);
		const DataArray1D<float> & dataParam  = varParam.m_data;

		for (int i = 0; i < m_data.GetRows(); i++) {
			m_data[i] = fabs(dataParam[i]);
		}

	// Evaluate the average operator
	} else if (m_strName == "_AVG") {
		if (m_varArg.size() <= 1) {
			_EXCEPTION1("_AVG expects at least two arguments: %i given",
				m_varArg.size());
		}

		m_data.Zero();
		for (int v = 0; v < m_varArg.size(); v++) {
			Variable & varParam = varreg.Get(m_varArg[v]);
			varParam.LoadGridData(varreg, vecFiles, grid, sTime);
			const DataArray1D<float> & dataParam  = varParam.m_data;

			for (int i = 0; i < m_data.GetRows(); i++) {
				m_data[i] += dataParam[i];
			}
		}

		for (int i = 0; i < m_data.GetRows(); i++) {
			m_data[i] /= static_cast<double>(m_varArg.size());
		}

	// Evaluate the minus operator
	} else if (m_strName == "_DIFF") {
		if (m_varArg.size() != 2) {
			_EXCEPTION1("_VECMAG expects two arguments: %i given",
				m_varArg.size());
		}
		Variable & varLeft = varreg.Get(m_varArg[0]);
		Variable & varRight = varreg.Get(m_varArg[1]);

		varLeft.LoadGridData(varreg, vecFiles, grid, sTime);
		varRight.LoadGridData(varreg, vecFiles, grid, sTime);

		const DataArray1D<float> & dataLeft  = varLeft.m_data;
		const DataArray1D<float> & dataRight = varRight.m_data;

		for (int i = 0; i < m_data.GetRows(); i++) {
			m_data[i] = dataLeft[i] - dataRight[i];
		}

	// Evaluate the Coriolis parameter operator
	} else if (m_strName == "_F") {
		if (m_varArg.size() != 0) {
			_EXCEPTION1("_F expects zero arguments: %i given",
				m_varArg.size());
		}

		for (int i = 0; i < m_data.GetRows(); i++) {
			m_data[i] = 2.0 * 7.2921e-5 * sin(mesh.dLat[i]);
		}

	// Evaluate the mean operator
	} else if (m_strName == "_MEAN") {
		if (m_varArg.size() != 2) {
			_EXCEPTION1("_MEAN expects two arguments: %i given",
				m_varArg.size());
		}

		// Check for adjacency list
		if (mesh.adjlist.size() == 0) {
			_EXCEPTIONT("Adjacency list not available for Grid");
		}

		// Obtain field and distance
		Variable & varField = varreg.Get(m_varArg[0]);
		Variable & varDist = varreg.Get(m_varArg[1]);

		varField.LoadGridData(varreg, vecFiles, grid, sTime);

		// Load distance (in degrees) and convert to radians
		double dDist = atof(varDist.m_strName.c_str());

		if ((dDist < 0.0) || (dDist > 360.0)) {
			_EXCEPTION1("Distance argument in _MEAN out of range\n"
				"Expected [0,360], found %1.3e", dDist);
		}

		// Calculate mean of field
		m_data.Zero();

		if (mesh.adjlist.size() != m_data.GetRows()) {
			_EXCEPTIONT("Invalid grid connectivity array");
		}

		for (int i = 0; i < m_data.GetRows(); i++) {
			std::set<int> setNodesVisited;
			std::set<int> setNodesToVisit;
			setNodesToVisit.insert(i);

			double dLat0 = mesh.dLat[i];
			double dLon0 = mesh.dLon[i];

			while (setNodesToVisit.size() != 0) {

				// Next node to explore
				int j = *(setNodesToVisit.begin());

				setNodesToVisit.erase(setNodesToVisit.begin());
				setNodesVisited.insert(j);

				// Update the mean
				m_data[i] += varField.m_data[j];

				// Find additional neighbors to explore
				std::set<int>::const_iterator iter = mesh.adjlist[j].begin();
				for (; iter != mesh.adjlist[j].end(); iter++) {
					int l = (*iter);

					// Check if already visited
					if (setNodesVisited.find(l) != setNodesVisited.end()) {
						continue;
					}

					// Test that this node satisfies the distance criteria
					double dLat1 = mesh.dLat[l];
					double dLon1 = mesh.dLon[l];

					double dR =
						sin(dLat0) * sin(dLat1)
						+ cos(dLat0) * cos(dLat1) * cos(dLon1 - dLon0);

					if (dR >= 1.0) {
						dR = 0.0;
					} else if (dR <= -1.0) {
						dR = 180.0;
					} else {
						dR = 180.0 / M_PI * acos(dR);
					}
					if (dR != dR) {
						_EXCEPTIONT("NaN value detected");
					}

					if (dR > dDist) {
						continue;
					}

					// Add node to visit
					setNodesToVisit.insert(l);
				}
			}

			// Average data to obtain mean
			if (setNodesVisited.size() == 0) {
				_EXCEPTIONT("Logic error");
			}

			m_data[i] /= static_cast<float>(setNodesVisited.size());
		}

	} else {
		_EXCEPTION1("Unexpected operator \"%s\"", m_strName.c_str());
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void Variable::UnloadGridData() {

	// Force data to be loaded within this structure
	m_sTime = (-2);
}

///////////////////////////////////////////////////////////////////////////////

