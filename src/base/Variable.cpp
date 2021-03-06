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
#include "GridObject.h"

///////////////////////////////////////////////////////////////////////////////
// VariableRegistry
///////////////////////////////////////////////////////////////////////////////

VariableRegistry::VariableRegistry(
	RecapConfigObject * pobjConfig
) :
	m_pobjConfig(pobjConfig)
{
	if (pobjConfig == NULL) {
		_EXCEPTIONT("Invalid RecapConfigObject pointer");
	}
}

///////////////////////////////////////////////////////////////////////////////

VariableRegistry::~VariableRegistry() {
	VariableMap::iterator iter = m_mapVariables.begin();
	for (; iter != m_mapVariables.end(); iter++) {
		delete iter->second;
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string VariableRegistry::FindOrRegister(
	const std::string & strVariableName,
	Variable ** ppVariable
) {
	if (strVariableName.length() == 0) {
		_EXCEPTIONT("Attempting to find/register zero length Variable name");
	}
	if (ppVariable == NULL) {
		_EXCEPTIONT("Invalid ppVariable argument");
	}

	// Chck if Variable exists in registry
	VariableMap::iterator iterVar = m_mapVariables.find(strVariableName);
	if (iterVar != m_mapVariables.end()) {
		(*ppVariable) = iterVar->second;

		return std::string("");
	}

	std::cout << "Registering " << m_pobjConfig->Name() << " : " << strVariableName << std::endl;

	// Check if Variable is in the FileList, and add to registry if present
	FileListObject * pobjFileList =
		m_pobjConfig->GetFileList();
	const VariableInfo * pvarinfo =
		pobjFileList->GetVariableInfo(strVariableName);
	if (pvarinfo != NULL) {
		Variable * pvar = new Variable(this, strVariableName, pvarinfo);

		std::pair<VariableMap::iterator, bool> prSuccess =
			m_mapVariables.insert(
				std::pair<std::string, Variable *>(strVariableName, pvar));

		if (!prSuccess.second) {
			_EXCEPTIONT("Repeated Variable found in Variable map");
		}

		iterVar = prSuccess.first;

		(*ppVariable) = iterVar->second;

		return std::string("");
	}

	// Check if Variable is an operator combination
	if (strVariableName[0] == '_') {

		// Register the operator combination
		Variable * pvarNew = new Variable(this, strVariableName, NULL);

		std::string strError = pvarNew->Initialize();
		if (strError != "") {
			return strError;
		}

		// Recursively register arguments to operator combination
		std::vector<Variable *> vecArguments;
		vecArguments.resize(pvarNew->m_vecArg.size());

		for (int a = 0; a < pvarNew->m_vecArg.size(); a++) {
			strError =
				m_pobjConfig->GetVariable(
					pvarNew->m_vecArg[a],
					&(vecArguments[a]));

			if (strError != "") {
				return strError;
			}
		}

		// Verify compatibility and inherit auxiliary data from arguments
		strError = pvarNew->InitializeAuxiliary(vecArguments);

		if (strError != "") {
			return strError;
		}

		// Insert new variable into registry
		std::pair<VariableMap::iterator, bool> prSuccess =
			m_mapVariables.insert(
				std::pair<std::string, Variable *>(
					strVariableName, pvarNew));

		if (!prSuccess.second) {
			_EXCEPTIONT("Repeated Variable found in Variable map");
		}

		iterVar = prSuccess.first;

		(*ppVariable) = iterVar->second;

		return std::string("");
	}

	// Check for special variable names
	if ((strVariableName == "lon") ||
	    (strVariableName == "lat")
	) {
		// Register the operator combination
		Variable * pvarNew = new Variable(this, strVariableName, NULL);
		pvarNew->m_strUnits = "deg";

		// Insert new variable into registry
		std::pair<VariableMap::iterator, bool> prSuccess =
			m_mapVariables.insert(
				std::pair<std::string, Variable *>(
					strVariableName, pvarNew));

		if (!prSuccess.second) {
			_EXCEPTIONT("Repeated Variable found in Variable map");
		}

		iterVar = prSuccess.first;

		(*ppVariable) = iterVar->second;

		return std::string("");
	}

	// Invalid variable name
	return std::string("Invalid variable name \"")
		+ strVariableName + std::string("\"");
}

///////////////////////////////////////////////////////////////////////////////

bool VariableRegistry::FindOrRegister(
	const Variable & var
) {
	_EXCEPTIONT("Deprecated -- don't call this routine");
/*
	std::map<std::string, Variable>::const_iterator iter =
		m_mapVariables.find(var.m_strName);

	if (iter == m_mapVariables.end()) {
		m_mapVariables.insert(
			std::pair<std::string, Variable>(var.m_strName, var));

		return true;
	}
*/
	return false;
}

///////////////////////////////////////////////////////////////////////////////

Variable * VariableRegistry::Get(
	const std::string & strName
) {
	std::map<std::string, Variable *>::iterator iter =
		m_mapVariables.find(strName);

	if (iter == m_mapVariables.end()) {
		return NULL;
	} else {
		return iter->second;
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::UnloadAllGridData() {
	std::map<std::string, Variable *>::iterator iter =
		m_mapVariables.begin();
	for (; iter != m_mapVariables.end(); iter++) {
		iter->second->UnloadAllGridData();
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::UpdateTimeIndices(
	const std::map<size_t, size_t> & mapTimeIxToNewTimeIx
) {
	_EXCEPTIONT("FIX");
/*
	for (size_t v = 0; v < m_vecVariables.size(); v++) {
		if (m_vecVariables[v].m_iTimeDimIx == (-1)) {
			continue;
		}

		size_t sTime =
			static_cast<size_t>(
				m_vecVariables[v].m_vecAuxIndices[m_vecVariables[v].m_iTimeDimIx]);

		std::map<size_t, size_t>::const_iterator iterTimeIxMap =
			mapTimeIxToNewTimeIx.find(sTime);

		if (iterTimeIxMap == mapTimeIxToNewTimeIx.end()) {
			_EXCEPTIONT("Invalid Variable time index");
		}

		m_vecVariables[v].m_vecAuxIndices[m_vecVariables[v].m_iTimeDimIx] =
			static_cast<long>(iterTimeIxMap->second);
	}
*/
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
	m_fReductionOp(false),
	m_strUnits(),
	m_iTimeDimIx(-1),
	m_iVerticalDimIx(-1),
	m_nVerticalDimOrder(+1),
	m_nSpecifiedDim(0),
	m_fLimitToSingleDataInstance(true)
{
	if (strName.length() == 0) {
		_EXCEPTIONT("Variable name must be given");
	}
	if (pvarreg == NULL) {
		_EXCEPTIONT("pvarreg must be specified in Variable() constuctor");
	}
	if (pvarinfo != NULL) {
		m_strUnits = pvarinfo->m_strUnits;
		m_vecAuxDimNames = pvarinfo->m_vecAuxDimNames;
		m_iTimeDimIx = pvarinfo->m_iTimeDimIx;
		m_iVerticalDimIx = pvarinfo->m_iVerticalDimIx;
		m_nVerticalDimOrder = pvarinfo->m_nVerticalDimOrder;

		//std::cout << "TEST " << m_strName << std::endl;
		//for (int d = 0; d < pvarinfo->m_vecAuxDimSizes.size(); d++) {
		//	std::cout << "  " << pvarinfo->m_vecAuxDimNames[d] << " " << pvarinfo->m_vecAuxDimSizes[d] << std::endl;
		//}
		m_iterAuxBegin.Initialize(pvarinfo->m_vecAuxDimSizes, false);
		m_iterAuxEnd.Initialize(pvarinfo->m_vecAuxDimSizes, true);
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::Initialize() {

	// Operator
	if (m_strName[0] == '_') {

		// Set operator flag
		m_fOp = true;

		// Parse arguments
		int iArgBegin = 0;
		int iArgEnd = m_strName.length()-1;

		for (; iArgBegin < m_strName.length(); iArgBegin++) {
			if (m_strName[iArgBegin] == '(') {
				break;
			}
		}
		if (iArgBegin == m_strName.length()) {
			return
				std::string("ERROR: Variable operator")
				+ m_strName + std::string(" missing open parenthesis");
		}
	
		for (; iArgEnd >= 0; iArgEnd--) {
			if (m_strName[iArgEnd] == ')') {
				break;
			}
		}
		if (iArgEnd < 0) {
			return
				std::string("ERROR: Variable operator")
				+ m_strName + std::string(" missing close parenthesis");
		}
	
		// Operator name
		m_strOpName = m_strName.substr(0, iArgBegin);
	
		// Parse argument list
		int iPos = iArgBegin+1;
		int iParenthesisNest = 0;
		for (int i = iArgBegin+1; i <= iArgEnd; i++) {
			if (m_strName[i] == '(') {
				iParenthesisNest++;
			}
			if (((m_strName[i] == ',') && (iParenthesisNest == 0)) ||
			    ((m_strName[i] == ')') && (i == iArgEnd) && (iParenthesisNest == 0))
			) {
				// Add argument to Variable's argument list
				std::string strArg = m_strName.substr(iPos, i-iPos);
				m_vecArg.push_back(strArg);
	
				iPos = i+1;
	
				if (i == iArgEnd) {
					break;
				}
			}
	
			if (m_strName[i] == ')') {
				if (iParenthesisNest > 0) {
					iParenthesisNest--;
				} else {
					return
						std::string("ERROR: Variable operator")
						+ m_strName + std::string(" has unbalanced parentheses");
				}
			}
		}

		// Check operator properties
		if (m_strOpName == "_CLIMMEAN") {
			m_fReductionOp = true;
			if (m_vecArg.size() != 1) {
				return std::string("ERROR: _CLIMMEAN expects exactly one argument");
			}
		}
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::InitializeAuxiliary(
	const std::vector<Variable *> vecVarArguments
) {
	// Inherit attributes from argument by default
	if (vecVarArguments.size() != 0) {

		m_strUnits = vecVarArguments[0]->m_strUnits;
		m_vecAuxDimNames = vecVarArguments[0]->m_vecAuxDimNames;
		m_iTimeDimIx = vecVarArguments[0]->m_iTimeDimIx;
		m_iVerticalDimIx = vecVarArguments[0]->m_iVerticalDimIx;
		m_nVerticalDimOrder = vecVarArguments[0]->m_nVerticalDimOrder;

		m_iterAuxBegin = vecVarArguments[0]->m_iterAuxBegin;
		m_iterAuxEnd = vecVarArguments[0]->m_iterAuxEnd;

		if (m_iterAuxBegin.size() != m_vecAuxDimNames.size()) {
			_EXCEPTIONT("Logic error");
		}
	}

	// Reduction operator
	if (m_strOpName == "_CLIMMEAN") {
		if (vecVarArguments.size() != 1) {
			return std::string("_CLIMMEAN expectes one argument");
		}

		// Remove time dimension
		m_vecAuxDimNames.erase(m_vecAuxDimNames.begin() + m_iTimeDimIx);

		m_iterAuxBegin.RemoveDim(m_iTimeDimIx);
		m_iterAuxEnd.RemoveDim(m_iTimeDimIx);

		m_iTimeDimIx = (-1);

		return std::string("");
	}

	// Vector magnitude
	if (m_strOpName == "_VECMAG") {
		if (vecVarArguments.size() != 2) {
			return std::string("_VECMAG expectes two arguments");
		}

		return std::string("");
	}

	// Vector difference
	if (m_strOpName == "_DIFF") {
		if (vecVarArguments.size() != 2) {
			return std::string("_DIFF expectes two arguments");
		}

		return std::string("");
	}

	// Vector product
	if (m_strOpName == "_PROD") {
		if (vecVarArguments.size() != 2) {
			return std::string("_PROD expectes two arguments");
		}

		return std::string("");
	}

	_EXCEPTIONT("Missing auxiliary specification for op");
}

///////////////////////////////////////////////////////////////////////////////

bool Variable::operator==(
	const Variable & var
) {
	_EXCEPTIONT("Not available");
/*
	if (m_strName != var.m_strName) {
		return false;
	}
	if (m_strUnits != var.m_strUnits) {
		return false;
	}
	if (m_vecAuxDimNames != var.m_vecAuxDimNames) {
		return false;
	}
	if (m_fOp != var.m_fOp) {
		return false;
	}
	if (m_fReductionOp != var.m_fReductionOp) {
		return false;
	}
	if (m_nSpecifiedDim != var.m_nSpecifiedDim) {
		return false;
	}
	return true;
*/
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::ToString() const {

	// Verify initialization
	if (m_pvarreg == NULL) {
		_EXCEPTIONT("Invalid VariableRegistry");
	}
/*
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
*/
	_EXCEPTION();
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::LoadGridData(
	size_t sTime,
	DataArray1D<float> ** ppdata,
	bool * pfReloaded
) {
	if (m_iTimeDimIx == (-1)) {
		return LoadGridData(m_iterAuxBegin, ppdata, pfReloaded);

	} else if (m_iterAuxBegin.size() != 1) {
		_EXCEPTIONT("Multiple auxiliary indices found");

	} else {
		VariableAuxIndex ixAux = m_iterAuxBegin;
		ixAux[m_iTimeDimIx] = static_cast<long>(sTime);
	
		return LoadGridData(ixAux, ppdata, pfReloaded);
	}
}

///////////////////////////////////////////////////////////////////////////////
/*
std::string Variable::MultiLoadGridData(
	const std::vector<VariableAuxIndex> & vecMultiIndex
) {
	UnloadAllGridData();

	for (int d = 0; d < vecMultiIndex.size(); d++) {
		DataArray1D<float> * pdata = new DataArray1D<float>;
		if (pdata == NULL) {
			_EXCEPTIONT("Allocation failure");
		}

		m_mapData.insert(DataMap::value_type(vecMultiIndex[d], pdata));

		std::string strError = _LoadGridData(vecMultiIndex[d], NULL);
		if (strError != "") {
			return strError;
		}
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::MultiLoadGridData3D(
	const VariableAuxIndex & ixAux
) {
	if (ixAux.size() != m_vecAuxDimNames.size()-1) {
		_EXCEPTIONT("Invalid length of auxiliary index");
	}
	if (m_iVerticalDimIx == (-1)) {
		_EXCEPTIONT("Not a 3D variable -- no vertical dimension found");
	}

	// Number of vertical levels
	long lVeticalDimSize = VerticalDimSize();

	// Build the multi-index
	std::vector<VariableAuxIndex> vecMultiIndex;
	vecMultiIndex.resize(lVeticalDimSize);

	for (long k = 0; k < lVeticalDimSize; k++) {
		vecMultiIndex[k] = ixAux;
		vecMultiIndex[k].insert(
			vecMultiIndex[k].begin() + m_iVerticalDimIx, k);
	}

	return MultiLoadGridData(vecMultiIndex);
}
*/
///////////////////////////////////////////////////////////////////////////////

std::string Variable::LoadGridData(
	const VariableAuxIndex & ixAux,
	DataArray1D<float> ** ppdata,
	bool * pfReloaded
) {
	std::string strError;

	// Assume data needs to be reloaded
	if (pfReloaded != NULL) {
		(*pfReloaded) = true;
	}

	// Check for single time index (variable already loaded)
	if (m_iterAuxBegin.size() != ixAux.size()) {
		_EXCEPTION2("Auxiliary index size mismatch (%lu/%lu)",
			m_iterAuxBegin.size(), ixAux.size());
	}

	// Get a pointer to the RecapConfigObject
	RecapConfigObject * pobjConfig = m_pvarreg->GetRecapConfigObject();

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
	DataMap::iterator iterData = m_mapData.find(ixAux);
	if (iterData != m_mapData.end()) {
		if (ppdata != NULL) {
			(*ppdata) = iterData->second;
		}
		if (pfReloaded != NULL) {
			(*pfReloaded) = false;
		}
		return std::string("");
	}

	// Create a new data array and insert into map
	if (m_fLimitToSingleDataInstance) {
		UnloadAllGridData();
	}

	DataArray1D<float> * pdata = new DataArray1D<float>;
	if (pdata == NULL) {
		_EXCEPTIONT("Out of memory");
	}

	std::pair<DataMap::iterator,bool> prSuccess =
		m_mapData.insert(DataMap::value_type(ixAux, pdata));
	if (!prSuccess.second) {
		_EXCEPTIONT("Unable to insert new datafile into DataMap");
	}
	if (ppdata != NULL) {
		(*ppdata) = pdata;
	}

	//std::cout << "Now contains " << m_mapData.size() << " data instances" << std::endl;

	// Reference to data
	DataArray1D<float> & data = (*pdata);

	// Allocate data
	data.Allocate(mesh.sDOFCount);

	// Check for special variable names
	if (m_pvarinfo == NULL) {
		if ((m_strName == "lon") && (mesh.dCenterLon.GetRows() == mesh.sDOFCount)) {
			Announce("Extracting grid data [%s] [NoTime]", m_strName.c_str());
			for (size_t i = 0; i < mesh.sDOFCount; i++) {
				data[i] = static_cast<float>(mesh.dCenterLon[i]);
			}

			return std::string("");
		}
		if ((m_strName == "lat") && (mesh.dCenterLat.GetRows() == mesh.sDOFCount)) {
			Announce("Extracting grid data [%s] [NoTime]", m_strName.c_str());
			for (size_t i = 0; i < mesh.sDOFCount; i++) {
				data[i] = static_cast<float>(mesh.dCenterLat[i]);
			}

			return std::string("");
		}
	}

	// Get the data directly from a variable
	if (!m_fOp) {
		strError =
			pobjFileList->LoadData_float(
				m_strName,
				ixAux,
				data);

		if (strError != "") return strError;

		return std::string("");
	}

	// Evaluate the climate mean operator
	if (m_strOpName == "_CLIMMEAN") {

		if (m_vecArg.size() != 1) {
			_EXCEPTION1("_CLIMMEAN expects one argument: %i given",
				m_vecArg.size());
		}

		Variable * pvar = NULL;
		strError = pobjConfig->GetVariable(m_vecArg[0], &pvar);
		if (strError != "") return strError;

		if (pvar->m_iterAuxBegin.size()-1 != ixAux.size()) {
			_EXCEPTION2("Dimension size mismatch in reduce operation (%lu / %lu)",
				pvar->m_iterAuxBegin.size()-1, ixAux.size());
		}

		// Calculate climate mean
		size_t sTimeCount = pobjFileList->GetTimeCount();
		for (size_t t = 0; t < sTimeCount; t++) {
			VariableAuxIndex ixAuxArg = ixAux;
			ixAuxArg.insert(ixAuxArg.begin() + pvar->m_iTimeDimIx, t);

			DataArray1D<float> * pdataVar = NULL;
			strError = pvar->LoadGridData(ixAuxArg, &pdataVar);
			if (strError != "") return strError;

			for (int i = 0; i < data.GetRows(); i++) {
				data[i] += (*pdataVar)[i];
			}
		}
		for (int i = 0; i < data.GetRows(); i++) {
			data[i] /= static_cast<double>(sTimeCount);
		}

	// Evaluate the vector magnitude operator
	} else if (m_strOpName == "_VECMAG") {
		if (m_vecArg.size() != 2) {
			_EXCEPTION1("_VECMAG expects two arguments: %i given",
				m_vecArg.size());
		}
		Variable * pvarLeft = NULL;
		Variable * pvarRight = NULL;

		strError = pobjConfig->GetVariable(m_vecArg[0], &pvarLeft);
		if (strError != "") return strError;

		strError = pobjConfig->GetVariable(m_vecArg[1], &pvarRight);
		if (strError != "") return strError;

		DataArray1D<float> * pdataLeft = NULL;
		strError = pvarLeft->LoadGridData(ixAux, &pdataLeft);
		if (strError != "") return strError;

		DataArray1D<float> * pdataRight = NULL;
		strError = pvarRight->LoadGridData(ixAux, &pdataRight);
		if (strError != "") return strError;

		for (int i = 0; i < data.GetRows(); i++) {
			data[i] =
				sqrt((*pdataLeft)[i] * (*pdataLeft)[i]
					+ (*pdataRight)[i] * (*pdataRight)[i]);
		}

	// Evaluate the difference operator
	} else if (m_strOpName == "_DIFF") {
		if (m_vecArg.size() != 2) {
			_EXCEPTION1("_DIFF expects two arguments: %i given",
				m_vecArg.size());
		}
		Variable * pvarLeft = NULL;
		Variable * pvarRight = NULL;

		strError = pobjConfig->GetVariable(m_vecArg[0], &pvarLeft);
		if (strError != "") return strError;

		strError = pobjConfig->GetVariable(m_vecArg[1], &pvarRight);
		if (strError != "") return strError;

		DataArray1D<float> * pdataLeft = NULL;
		strError = pvarLeft->LoadGridData(ixAux, &pdataLeft);
		if (strError != "") return strError;

		DataArray1D<float> * pdataRight = NULL;
		strError = pvarRight->LoadGridData(ixAux, &pdataRight);
		if (strError != "") return strError;

		for (int i = 0; i < data.GetRows(); i++) {
			data[i] = (*pdataLeft)[i] - (*pdataRight)[i];
		}

	// Evaluate the product operator
	} else if (m_strOpName == "_PROD") {
		if (m_vecArg.size() != 2) {
			_EXCEPTION1("_PROD expects two arguments: %i given",
				m_vecArg.size());
		}
		Variable * pvarLeft = NULL;
		Variable * pvarRight = NULL;

		strError = pobjConfig->GetVariable(m_vecArg[0], &pvarLeft);
		if (strError != "") return strError;

		strError = pobjConfig->GetVariable(m_vecArg[1], &pvarRight);
		if (strError != "") return strError;

		DataArray1D<float> * pdataLeft = NULL;
		strError = pvarLeft->LoadGridData(ixAux, &pdataLeft);
		if (strError != "") return strError;

		DataArray1D<float> * pdataRight = NULL;
		strError = pvarRight->LoadGridData(ixAux, &pdataRight);
		if (strError != "") return strError;

		for (int i = 0; i < data.GetRows(); i++) {
			data[i] = (*pdataLeft)[i] * (*pdataRight)[i];
		}

/*
	// Evaluate the absolute value operator
	} else if (m_strOpName == "_ABS") {
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
	} else if (m_strOpName == "_AVG") {
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

	// Evaluate the Coriolis parameter operator
	} else if (m_strOpName == "_F") {
		if (m_varArg.size() != 0) {
			_EXCEPTION1("_F expects zero arguments: %i given",
				m_varArg.size());
		}

		for (int i = 0; i < m_data.GetRows(); i++) {
			m_data[i] = 2.0 * 7.2921e-5 * sin(mesh.dLat[i]);
		}

	// Evaluate the mean operator
	} else if (m_strOpName == "_MEAN") {
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
*/

	} else {
		_EXCEPTION1("Unexpected operator \"%s\"", m_strName.c_str());
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::AllocateGridData(
	const VariableAuxIndex & ixAux,
	DataArray1D<float> ** ppdata
) {
	// Check aux index
	if (ixAux.size() != m_vecAuxDimNames.size()) {
		_EXCEPTIONT("Invalid auxiliary index");
	}

	// Get a pointer to the RecapConfigObject
	RecapConfigObject * pobjConfig = m_pvarreg->GetRecapConfigObject();

	// Get Mesh
	GridObject * pobjGrid = pobjConfig->GetGrid();
	if (pobjGrid == NULL) {
		return std::string("ERROR: Invalid configuration (missing grid)");
	}

	const Mesh & mesh = pobjGrid->GetMesh();

	// Find data instance in DataMap
	DataMap::iterator iterData = m_mapData.find(ixAux);
	if (iterData != m_mapData.end()) {
		(*ppdata) = iterData->second;
	}

	// Remove other data instances
	if (m_fLimitToSingleDataInstance) {
		UnloadAllGridData();
	}

	// Create a new data array and insert into map
	DataArray1D<float> * pdata = new DataArray1D<float>;
	if (pdata == NULL) {
		_EXCEPTIONT("Out of memory");
	}
	if (ppdata != NULL) {
		(*ppdata) = pdata;
	}

	pdata->Allocate(mesh.sDOFCount);

	std::pair<DataMap::iterator,bool> prSuccess =
		m_mapData.insert(DataMap::value_type(ixAux, pdata));
	if (!prSuccess.second) {
		_EXCEPTIONT("Unable to insert new datafile into DataMap");
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::WriteGridData(
	const VariableAuxIndex & ixAux
) {
	// Get a pointer to the RecapConfigObject
	RecapConfigObject * pobjConfig = m_pvarreg->GetRecapConfigObject();

	// Get FileList
	FileListObject * pobjFileList = pobjConfig->GetFileList();
	if (pobjFileList == NULL) {
		_EXCEPTIONT("Invalid configuration: Missing file_list");
	}

	DataMap::iterator iterData = m_mapData.find(ixAux);
	if (iterData == m_mapData.end()) {
		return std::string("ERROR: AuxIndex not found in data");
	}

	std::string strError =
		pobjFileList->WriteData_float(
			m_strName,
			ixAux,
			*(iterData->second));

	return strError;
}

///////////////////////////////////////////////////////////////////////////////

void Variable::UnloadAllGridData() {
	for (DataMap::iterator iter = m_mapData.begin(); iter != m_mapData.end(); iter++) {
/*
		{
			std::string strLoading =
				std::string("UNLOAD [") + m_strName + std::string("]");

			VariableAuxIndex vecAuxIndices = iter->first;
			for (size_t d = 0; d < vecAuxIndices.size(); d++) {
				strLoading +=
					std::string(" [")
					+ std::to_string(vecAuxIndices[d])
					+ std::string("]");
			}
		
			Announce(strLoading.c_str());
		}
*/
		iter->second->Deallocate();
		delete iter->second;
	}
	m_mapData.clear();
}

///////////////////////////////////////////////////////////////////////////////

void Variable::UnloadGridData(
	const VariableAuxIndex & ixAux,
	bool fNoThrow
) {
	DataMap::iterator iter = m_mapData.find(ixAux);
	if (iter == m_mapData.end()) {
		if (fNoThrow) {
			return;
		}
		_EXCEPTIONT("Data instance not found.");
	}
/*
	{
		std::string strLoading =
			std::string("UNLOAD [") + m_strName + std::string("]");

		for (size_t d = 0; d < ixAux.size(); d++) {
			strLoading +=
				std::string(" [")
				+ std::to_string(ixAux[d])
				+ std::string("]");
		}
	
		Announce(strLoading.c_str());
	}
*/
	iter->second->Deallocate();
	m_mapData.erase(iter);
}

///////////////////////////////////////////////////////////////////////////////

