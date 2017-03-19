///////////////////////////////////////////////////////////////////////////////
///
///	\file    PointSearch.cpp
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

#include "PointSearch.h"
#include "PointDataObject.h"

#include "Variable.h"
#include "Exception.h"
#include "Announce.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "TimeObj.h"

#include "RecapConfigObject.h"

#include "kdtree.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <set>
#include <queue>

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

namespace HPointwise {

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a thresholding operator.
///	</summary>
class ThresholdOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		GreaterThan,
		LessThan,
		GreaterThanEqualTo,
		LessThanEqualTo,
		EqualTo,
		NotEqualTo
	};

public:
	///	<summary>
	///		Initialize from an Object.
	///	</summary>
	std::string InitializeFromObject(
		const ObjectRegistry & objreg,
		RecapConfigObject * pobjConfig,
		const Object * pobj
	) {
		// Variable
		StringObject * pobjThresholdCmdVar =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("var")));
		if (pobjThresholdCmdVar == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"var\" property");
		}

		std::string strError =
			pobjConfig->GetVariable(
				pobjThresholdCmdVar->Value(),
				&m_pvar);

		if (strError != "") {
			return strError;
		}

		// Operator
		StringObject * pobjThresholdCmdOp =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("op")));
		if (pobjThresholdCmdOp == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"op\" property");
		}

		bool fSuccessThresholdOp =
			SetOperatorFromString(pobjThresholdCmdOp->Value());
		if (!fSuccessThresholdOp) {
			return std::string("Invalid threshold operator in ")
				+ pobj->ChildName("op");
		}

		// Value 
		StringObject * pobjThresholdValue =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("value")));
		if (pobjThresholdValue == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"value\" property");
		}

		bool fSuccessValue =
			pobjThresholdValue->ToUnit(
				m_pvar->Units(), &m_dValue);
		if (!fSuccessValue) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("value")
				+ std::string(" to ")
				+ m_pvar->Units();
		}

		// Distance
		StringObject * pobjThresholdCmdDist =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("dist")));
		if (pobjThresholdCmdDist == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"dist\" property");
		}

		bool fSuccessDist =
			pobjThresholdCmdDist->ToUnit(
				"deg", &m_dDistance);
		if (!fSuccessDist) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("dist")
				+ std::string(" to great-circle distance");
		}

		return std::string("");
	}

public:
	///	<summary>
	///		Convert to string.
	///	</summary>
	std::string ToString() const {

		if (m_pvar == NULL) {
			return std::string("Invalid ThresholdOp");
		}

		// Output announcement
		std::string strDescription = m_pvar->Name();
		if (m_eOp == GreaterThan) {
			strDescription += " is greater than ";
		} else if (m_eOp == LessThan) {
			strDescription += " is less than ";
		} else if (m_eOp == GreaterThanEqualTo) {
			strDescription += " is greater than or equal to ";
		} else if (m_eOp == LessThanEqualTo) {
			strDescription += " is less than or equal to ";
		} else if (m_eOp == EqualTo) {
			strDescription += " is equal to ";
		} else if (m_eOp == NotEqualTo) {
			strDescription += " is not equal to ";
		}

		strDescription += std::to_string(m_dValue)
			+ std::string(" within ")
			+ std::to_string(m_dDistance)
			+ std::string(" degrees");

		return strDescription;
	}

	///	<summary>
	///		Set the operator from a string.
	///	</summary>
	bool SetOperatorFromString(
		const std::string & strOp
	) {
		if (strOp == ">") {
			m_eOp = GreaterThan;
		} else if (strOp == "<") {
			m_eOp = LessThan;
		} else if (strOp == ">=") {
			m_eOp = GreaterThanEqualTo;
		} else if (strOp == "<=") {
			m_eOp = LessThanEqualTo;
		} else if (strOp == "=") {
			m_eOp = EqualTo;
		} else if (strOp == "!=") {
			m_eOp = NotEqualTo;
		} else {
			return false;
		}
		return true;
	}

public:
	///	<summary>
	///		Pointer to variable to use for threshold op.
	///	</summary>
	Variable * m_pvar;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Distance to search for threshold value (deg)
	///	</summary>
	double m_dDistance;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class describing a general closed contour operation.
///	</summary>
class ClosedContourOp {

public:
	///	<summary>
	///		Initialize from an Object.
	///	</summary>
	std::string InitializeFromObject(
		const ObjectRegistry & objreg,
		RecapConfigObject * pobjConfig,
		const Object * pobj
	) {
		// Variable
		StringObject * pobjClosedContourVar =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("var")));
		if (pobjClosedContourVar == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"var\" property");
		}

		std::string strError =
			pobjConfig->GetVariable(
				pobjClosedContourVar->Value(),
				&m_pvar);

		if (strError != "") {
			return strError;
		}

		// Magnitude
		StringObject * pobjClosedContourMag =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("mag")));
		if (pobjClosedContourMag == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"mag\" property");
		}

		bool fSuccessMag =
			pobjClosedContourMag->ToUnit(
				m_pvar->Units(), &m_dDeltaAmount);
		if (!fSuccessMag) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("mag")
				+ std::string(" to ")
				+ m_pvar->Units();
		}

		// Distance
		StringObject * pobjClosedContourDist =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("dist")));
		if (pobjClosedContourDist == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"dist\" property");
		}

		bool fSuccessDist =
			pobjClosedContourDist->ToUnit(
				"deg", &m_dDistance);
		if (!fSuccessDist) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("dist")
				+ std::string(" to great-circle distance");
		}

		// Search distance
		StringObject * pobjClosedContourMinMaxDist =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("minmaxdist")));
		if (pobjClosedContourMinMaxDist == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"minmaxdist\" property");
		}

		bool fSuccessMinMaxDist =
			pobjClosedContourMinMaxDist->ToUnit(
				"deg", &m_dMinMaxDist);
		if (!fSuccessMinMaxDist) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("minmaxdist")
				+ std::string(" to great-circle distance");
		}

		return std::string("");
	}

	///	<summary>
	///		Convert to string.
	///	</summary>
	std::string ToString() const {

		if (m_pvar == NULL) {
			return std::string("Invalid ClosedContourOp");
		}

		if (m_dDeltaAmount < 0.0) {
			std::string strDescription =
				m_pvar->Name() + std::string(" decreases by ")
				+ std::to_string(-m_dDeltaAmount)
				+ m_pvar->Units()
				+ std::string(" over ")
				+ std::to_string(m_dDistance)
				+ std::string("deg (max search ")
				+ std::to_string(m_dMinMaxDist)
				+ std::string("deg)");

			return strDescription;

		} else {
			std::string strDescription =
				m_pvar->Name() + std::string(" increases by ")
				+ std::to_string(m_dDeltaAmount)
				+ m_pvar->Units()
				+ std::string(" over ")
				+ std::to_string(m_dDistance)
				+ std::string("deg (min search ")
				+ std::to_string(m_dMinMaxDist)
				+ std::string("deg)");

			return strDescription;
		}

		return std::string("");
	}

public:
	///	<summary>
	///		Pointer to variable to use for closed contour op.
	///	</summary>
	Variable * m_pvar;

	///	<summary>
	///		Threshold amount.  If positive this represents a minimum
	///		increase.  If negative this represents a minimum decrease.
	///	</summary>
	double m_dDeltaAmount;

	///	<summary>
	///		Threshold distance.
	///	</summary>
	double m_dDistance;

	///	<summary>
	///		Distance to search for min or max.
	///	</summary>
	double m_dMinMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing an output operator.
///	</summary>
class OutputOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		Max,
		Min,
		Avg,
		MaxDist,
		MinDist
	};

public:
	///	<summary>
	///		Initialize from an Object.
	///	</summary>
	std::string InitializeFromObject(
		const ObjectRegistry & objreg,
		RecapConfigObject * pobjConfig,
		const Object * pobj
	) {
		// Initialize output from a single variable
		const StringObject * pobjString =
			dynamic_cast<const StringObject *>(pobj);
		if (pobjString != NULL) {
			std::string strError =
				pobjConfig->GetVariable(
					pobjString->Value(),
					&m_pvar);

			if (strError != "") {
				return strError;
			}

			m_eOp = Max;
			m_dDistance = 0.0;

			return std::string("");
		}

		// Variable
		StringObject * pobjOutputCmdVar =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("var")));
		if (pobjOutputCmdVar == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"var\" property");
		}

		std::string strError =
			pobjConfig->GetVariable(
				pobjOutputCmdVar->Value(),
				&m_pvar);

		if (strError != "") {
			return strError;
		}

		// Operator
		StringObject * pobjOutputCmdOp =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("op")));
		if (pobjOutputCmdOp == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"op\" property");
		}

		bool fSuccessOutputOp =
			SetOperatorFromString(pobjOutputCmdOp->Value());
		if (!fSuccessOutputOp) {
			return std::string("Invalid output operator in ")
				+ pobj->ChildName("op");
		}

		// Distance
		StringObject * pobjOutputCmdDist =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("dist")));
		if (pobjOutputCmdDist == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"dist\" property");
		}

		bool fSuccessDist =
			pobjOutputCmdDist->ToUnit(
				"deg", &m_dDistance);
		if (!fSuccessDist) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("dist")
				+ std::string(" to great-circle distance");
		}

		return std::string("");
	}

public:
	///	<summary>
	///		Convert to string.
	///	</summary>
	std::string ToString() const {

		if (m_pvar == NULL) {
			return std::string("Invalid OutputOp");
		}

		std::string strDescription;

		if (m_eOp == Max) {
			strDescription += "Maximum of ";
		} else if (m_eOp == Min) {
			strDescription += "Minimum of ";
		} else if (m_eOp == Avg) {
			strDescription += "Average of ";
		} else if (m_eOp == MaxDist) {
			strDescription += "Distance to maximum of ";
		} else if (m_eOp == MinDist) {
			strDescription += "Distance to minimum of ";
		}

		strDescription +=
			m_pvar->Name()
			+ std::string(" within ")
			+ std::to_string(m_dDistance)
			+ " degrees";

		return strDescription;
	}

	///	<summary>
	///		Assign the operator from a string.
	///	</summary>
	bool SetOperatorFromString(
		const std::string & strOp
	) {
		if (strOp == "max") {
			m_eOp = Max;
		} else if (strOp == "min") {
			m_eOp = Min;
		} else if (strOp == "avg") {
			m_eOp = Avg;
		} else if (strOp == "maxdist") {
			m_eOp = MaxDist;
		} else if (strOp == "mindist") {
			m_eOp = MinDist;
		} else {
			return false;
		}
		return true;
	}

public:
	///	<summary>
	///		Pointer to variable to use for output.
	///	</summary>
	Variable * m_pvar;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Distance to use when applying operation.
	///	</summary>
	double m_dDistance;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all minima in the given DataArray1D.
///	</summary>
template <typename real>
void FindAllLocalMinima(
	const Mesh & mesh,
	const DataArray1D<real> & data,
	std::set<int> & setMinima
) {
	const int sFaces = mesh.adjlist.size();
	for (int f = 0; f < sFaces; f++) {
		
		bool fMinimum = true;

		real dValue = data[f];
		int sNeighbors = mesh.adjlist[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (data[mesh.adjlist[f][n]] < dValue) {
				fMinimum = false;
				break;
			}
		}

		if (fMinimum) {
			setMinima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all maxima in the given DataArray2D.
///	</summary>
template <typename real>
void FindAllLocalMaxima(
	const Mesh & mesh,
	const DataArray1D<real> & data,
	std::set<int> & setMaxima
) {
	const int sFaces = mesh.adjlist.size();
	for (int f = 0; f < sFaces; f++) {

		bool fMaximum = true;

		real dValue = data[f];
		int sNeighbors = mesh.adjlist[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (data[mesh.adjlist[f][n]] > dValue) {
				fMaximum = false;
				break;
			}
		}

		if (fMaximum) {
			setMaxima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the minimum/maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
template <typename real>
void FindLocalMinMax(
	const Mesh & mesh,
	bool fMinimum,
	const DataArray1D<real> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	real & dMaxValue,
	float & dRMax
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Check index argument
	if (ix0 >= mesh.dLat.GetRows()) {
		_EXCEPTIONT("ix0 out of range of mesh.dLat");
	}
	if (ix0 >= mesh.dLon.GetRows()) {
		_EXCEPTIONT("ix0 out of range of mesh.dLon");
	}

	// Initialize the maximum to the central location
	ixExtremum = ix0;
	dMaxValue = data[ix0];
	dRMax = 0.0;

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ixExtremum);

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Latitude and longitude at the origin (in radians)
	double dLat0 = mesh.dLat[ix0] * M_PI / 180.0;
	double dLon0 = mesh.dLon[ix0] * M_PI / 180.0;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		double dLatThis = mesh.dLat[ix] * M_PI / 180.0;
		double dLonThis = mesh.dLon[ix] * M_PI / 180.0;

		// Great circle distance to this element
		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

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

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local extremum
		if (fMinimum) {
			if (data[ix] < dMaxValue) {
				ixExtremum = ix;
				dMaxValue = data[ix];
				dRMax = dR;
			}

		} else {
			if (data[ix] > dMaxValue) {
				ixExtremum = ix;
				dMaxValue = data[ix];
				dRMax = dR;
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < mesh.adjlist[ix].size(); n++) {
			queueNodes.push(mesh.adjlist[ix][n]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field has a closed contour about this point.
///	</summary>
template <typename real>
bool HasClosedContour(
	const Mesh & mesh,
	const DataArray1D<real> & dataState,
	const int ix0,
	double dDeltaAmt,
	double dDeltaDist,
	double dMinMaxDist
) {
	// Verify arguments
	if (dDeltaAmt == 0.0) {
		_EXCEPTIONT("Closed contour amount must be non-zero");
	}
	if (dDeltaDist <= 0.0) {
		_EXCEPTIONT("Closed contour distance must be positive");
	}

	// Find min/max near point
	int ixOrigin;

	if (dMinMaxDist == 0.0) {
		ixOrigin = ix0;

	// Find a local minimum / maximum
	} else {
		real dValue;
		float dR;

		FindLocalMinMax<real>(
			mesh,
			(dDeltaAmt > 0.0),
			dataState,
			ix0,
			dMinMaxDist,
			ixOrigin,
			dValue,
			dR);
	}

	//printf("%lu %lu : %lu %lu : %1.5f %1.5f\n", ix0 % grid.m_nGridDim[1], ix0 / grid.m_nGridDim[1], ixOrigin % grid.m_nGridDim[1], ixOrigin / grid.m_nGridDim[1], dataState[ix0], dataState[ixOrigin]);

	// Set of visited nodes
	std::set<int> setNodesVisited;

	// Set of nodes to visit
	std::queue<int> queueToVisit;
	queueToVisit.push(ixOrigin);

	// Reference value
	real dRefValue = dataState[ixOrigin];

	// Latitude and longitude at origin (in radians)
	const double dLat0 = mesh.dLat[ixOrigin] * M_PI / 180.0;
	const double dLon0 = mesh.dLon[ixOrigin] * M_PI / 180.0;

	Announce(2, "Checking (%lu) : (%1.5f %1.5f)",
		ixOrigin, dLat0, dLon0);

	// Build up nodes
	while (queueToVisit.size() != 0) {
		int ix = queueToVisit.front();
		queueToVisit.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element
		double dLatThis = mesh.dLat[ix] * M_PI / 180.0;
		double dLonThis = mesh.dLon[ix] * M_PI / 180.0;

		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

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

		Announce(2, "-- (%lu) : (%1.5f %1.5f) : dx %1.5f",
			ix, dLatThis, dLonThis, dR);

		//printf("-- (%lu %lu) %1.5f %1.5f\n", ix % grid.m_nGridDim[1], ix / grid.m_nGridDim[1], dR, dataState[ix] - dRefValue);

		// Check great circle distance
		if (dR > dDeltaDist) {
			Announce(2, "Failed criteria; returning");
			AnnounceEndBlock(2, NULL);
			return false;
		}

		// Verify sufficient increase in value
		if (dDeltaAmt > 0.0) {
			if (dataState[ix] - dRefValue >= dDeltaAmt) {
				continue;
			}

		// Verify sufficient decrease in value
		} else {
			if (dRefValue - dataState[ix] >= -dDeltaAmt) {
				continue;
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < mesh.adjlist[ix].size(); n++) {
			queueToVisit.push(mesh.adjlist[ix][n]);
		}
	}

	// Report success with criteria
	Announce(2, "Passed criteria; returning");
	AnnounceEndBlock(2, NULL);
	return true;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field satisfies the threshold.
///	</summary>
template <typename real>
bool SatisfiesThreshold(
	const Mesh & mesh,
	const DataArray1D<real> & dataState,
	const int ix0,
	const ThresholdOp::Operation op,
	const double dTargetValue,
	const double dMaxDist
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ix0);

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = mesh.dLat[ix0] * M_PI / 180.0;
	double dLon0 = mesh.dLon[ix0] * M_PI / 180.0;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element (in degrees)
		double dLatThis = mesh.dLat[ix] * M_PI / 180.0;
		double dLonThis = mesh.dLon[ix] * M_PI / 180.0;

		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

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

		if ((ix != ix0) && (dR > dMaxDist)) {
			continue;
		}

		// Value at this location
		double dValue = dataState[ix];

		// Apply operator
		if (op == ThresholdOp::GreaterThan) {
			if (dValue > dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThan) {
			if (dValue < dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::GreaterThanEqualTo) {
			if (dValue >= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThanEqualTo) {
			if (dValue <= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::EqualTo) {
			if (dValue == dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::NotEqualTo) {
			if (dValue != dTargetValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		// Special case: zero distance
		if (dMaxDist == 0.0) {
			return false;
		}

		// Add all neighbors of this point
		for (int n = 0; n < mesh.adjlist[ix].size(); n++) {
			queueNodes.push(mesh.adjlist[ix][n]);
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
template <typename real>
void FindLocalAverage(
	const Mesh & mesh,
	const DataArray1D<real> & data,
	int ix0,
	double dMaxDist,
	float & dAverage
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ix0);

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Latitude and longitude at the origin (in radians)
	double dLat0 = mesh.dLat[ix0] * M_PI / 180.0;
	double dLon0 = mesh.dLon[ix0] * M_PI / 180.0;

	// Number of points
	float dSum = 0.0;
	int nCount = 0;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		double dLatThis = mesh.dLat[ix] * M_PI / 180.0;
		double dLonThis = mesh.dLon[ix] * M_PI / 180.0;

		// Great circle distance to this element (in degrees)
		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

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

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local extremum
		dSum += data[ix];
		nCount++;

		// Add all neighbors of this point
		for (int n = 0; n < mesh.adjlist[ix].size(); n++) {
			queueNodes.push(mesh.adjlist[ix][n]);
		}
	}

	dAverage = dSum / static_cast<float>(nCount);
}

///////////////////////////////////////////////////////////////////////////////

class PointSearchParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PointSearchParam() :
		fpLog(NULL),
		strVariableSearchBy(),
		fSearchByMinima(false),
		dMaxLatitude(0.0),
		dMinLatitude(0.0),
		dMinAbsLatitude(0.0),
		dMaxLongitude(0.0),
		dMinLongitude(0.0),
		dMergeDist(0.0),
		pvecClosedContourOp(NULL),
		pvecNoClosedContourOp(NULL),
		pvecThresholdOp(NULL),
		pvecOutputOp(NULL),
		nTimeStride(1),
		fRegional(false),
		fOutputHeader(false),
		iVerbosityLevel(0)
	{ }

public:
	// Log
	FILE * fpLog;

	// Variable to search on
	std::string strVariableSearchBy;

	// Serach on minima
	bool fSearchByMinima;

	// Maximum latitude for detection
	double dMaxLatitude;

	// Minimum latitude for detection
	double dMinLatitude;

	// Minimum absolute value of latitude for detection
	double dMinAbsLatitude;

	// Maximum longitude for detection
	double dMaxLongitude;

	// Minimum longitude for detection
	double dMinLongitude;

	// Merge distance
	double dMergeDist;

	// Vector of closed contour operators
	std::vector<ClosedContourOp> * pvecClosedContourOp;

	// Vector of no closed contour operators
	std::vector<ClosedContourOp> * pvecNoClosedContourOp;

	// Vector of threshold operators
	std::vector<ThresholdOp> * pvecThresholdOp;

	// Vector of output operators
	std::vector<OutputOp> * pvecOutputOp;

	// Time stride
	int nTimeStride;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Output header
	bool fOutputHeader;

	// Verbosity level
	int iVerbosityLevel;

};

///////////////////////////////////////////////////////////////////////////////

std::string PointSearch(
	RecapConfigObject * pobjConfig,
	size_t sTime,
	const PointSearchParam & param,
	PointDataObject * pobjPointData
) {
	// Set the Announce buffer
	if (param.fpLog == NULL) {
		_EXCEPTIONT("Invalid log buffer");
	}

	AnnounceSetOutputBuffer(param.fpLog);
	//AnnounceOutputOnAllRanks();

	// Check minimum longitude / latitude
	if ((param.dMinLongitude < 0.0) || (param.dMinLongitude >= 360.0)) {
		_EXCEPTIONT("minlon must be in the range [0,360]");
	}
	if ((param.dMaxLongitude < 0.0) || (param.dMaxLongitude >= 360.0)) {
		_EXCEPTIONT("maxlon must be in the range [0,360]");
	}
	if ((param.dMaxLatitude < -90.0) || (param.dMaxLatitude > 90.0)) {
		_EXCEPTIONT("maxlat must be in the range [-90,90]");
	}
	if ((param.dMinLatitude < -90.0) || (param.dMinLatitude > 90.0)) {
		_EXCEPTIONT("minlat must be in the range [-90,90]");
	}
	if ((param.dMinAbsLatitude < 0.0) || (param.dMinAbsLatitude > 90.0)) {
		_EXCEPTIONT("minabslat must be in the range [0,90]");
	}

	// Dereference pointers to operators
	std::vector<ClosedContourOp> & vecClosedContourOp =
		*(param.pvecClosedContourOp);

	std::vector<ClosedContourOp> & vecNoClosedContourOp =
		*(param.pvecNoClosedContourOp);

	std::vector<ThresholdOp> & vecThresholdOp =
		*(param.pvecThresholdOp);

	std::vector<OutputOp> & vecOutputOp =
		*(param.pvecOutputOp);

	// Get a reference to the Mesh
	const Mesh & mesh = pobjConfig->GetGrid()->GetMesh();

	// Load the data for the search variable
	if (param.strVariableSearchBy == "") {
		_EXCEPTION();
	}

	Variable * pvarSearchBy = NULL;
	std::string strError =
		pobjConfig->GetVariable(
			param.strVariableSearchBy,
			&pvarSearchBy);

	if (strError != "") {
		return strError;
	}

	pvarSearchBy->LoadGridData(pobjConfig, sTime);

	const DataArray1D<float> & dataSearch = pvarSearchBy->GetData();

	// Tag all minima
	std::set<int> setCandidates;

	if (param.fSearchByMinima) {
		FindAllLocalMinima<float>(mesh, dataSearch, setCandidates);
	} else {
		FindAllLocalMaxima<float>(mesh, dataSearch, setCandidates);
	}

	// Total number of candidates
	int nTotalCandidates = setCandidates.size();

	int nRejectedLocation = 0;
	int nRejectedTopography = 0;
	int nRejectedMerge = 0;

	std::vector<int> vecRejectedClosedContour;
	vecRejectedClosedContour.resize(vecClosedContourOp.size());

	std::vector<int> vecRejectedNoClosedContour;
	vecRejectedNoClosedContour.resize(vecNoClosedContourOp.size());

	std::vector<int> vecRejectedThreshold;
	vecRejectedThreshold.resize(vecThresholdOp.size());

	// Eliminate based on interval
	if ((param.dMinLatitude != param.dMaxLatitude) ||
	    (param.dMinLongitude != param.dMaxLongitude) ||
		(param.dMinAbsLatitude != 0.0)
	) {
		std::set<int> setNewCandidates;

		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();
		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Get latitude and longitude in degrees
			double dLat = mesh.dLat[*iterCandidate];
			double dLon = mesh.dLon[*iterCandidate];

			if (param.dMinLatitude != param.dMaxLatitude) {
				if (dLat < param.dMinLatitude) {
					nRejectedLocation++;
					continue;
				}
				if (dLat > param.dMaxLatitude) {
					nRejectedLocation++;
					continue;
				}
			}
			if (param.dMinLongitude != param.dMaxLongitude) {
				if (dLon < 0.0) {
					int iLonShift = static_cast<int>(dLon / 360.0);
					dLon += static_cast<double>(iLonShift + 1) * 360.0;
				}
				if (dLon >= 360.0) {
					int iLonShift = static_cast<int>(dLon / 360.0);
					dLon -= static_cast<double>(iLonShift - 1) * 360.0;
				}
				if (param.dMinLongitude < param.dMaxLongitude) {
					if (dLon < param.dMinLongitude) {
						nRejectedLocation++;
						continue;
					}
					if (dLon > param.dMaxLongitude) {
						nRejectedLocation++;
						continue;
					}

				} else {
					if ((dLon > param.dMaxLongitude) &&
					    (dLon < param.dMinLongitude)
					) {
						nRejectedLocation++;
						continue;
					}
				}
			}
			if (param.dMinAbsLatitude != 0.0) {
				if (fabs(dLat) < param.dMinAbsLatitude) {
					nRejectedLocation++;
					continue;
				}
			}
			setNewCandidates.insert(*iterCandidate);
		}

		setCandidates = setNewCandidates;
	}

	// Eliminate based on merge distance
	if (param.dMergeDist != 0.0) {
		std::set<int> setNewCandidates;

		// Calculate chord distance
		double dSphDist =
			2.0 * sin(0.5 * param.dMergeDist / 180.0 * M_PI);

		// Create a new KD Tree containing all nodes
		kdtree * kdMerge = kd_create(3);

		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();
		for (; iterCandidate != setCandidates.end(); iterCandidate++) {
			double dLat = mesh.dLat[*iterCandidate] * M_PI / 180.0;
			double dLon = mesh.dLon[*iterCandidate] * M_PI / 180.0;

			double dX = cos(dLon) * cos(dLat);
			double dY = sin(dLon) * cos(dLat);
			double dZ = sin(dLat);

			kd_insert3(kdMerge, dX, dY, dZ, (void*)(&(*iterCandidate)));
		}

		// Loop through all candidates find set of nearest neighbors
		iterCandidate = setCandidates.begin();
		for (; iterCandidate != setCandidates.end(); iterCandidate++) {
			double dLat = mesh.dLat[*iterCandidate] * M_PI / 180.0;
			double dLon = mesh.dLon[*iterCandidate] * M_PI / 180.0;

			double dX = cos(dLon) * cos(dLat);
			double dY = sin(dLon) * cos(dLat);
			double dZ = sin(dLat);

			// Find all neighbors within dSphDist
			kdres * kdresMerge =
				kd_nearest_range3(kdMerge, dX, dY, dZ, dSphDist);

			// Number of neighbors
			int nNeighbors = kd_res_size(kdresMerge);
			if (nNeighbors == 0) {
				setNewCandidates.insert(*iterCandidate);

			} else {
				double dValue =
					static_cast<double>(dataSearch[*iterCandidate]);

				bool fExtrema = true;
				for (;;) {
					int * ppr = (int *)(kd_res_item_data(kdresMerge));

					if (param.fSearchByMinima) {
						if (static_cast<double>(dataSearch[*ppr]) < dValue) {
							fExtrema = false;
							break;
						}

					} else {
						if (static_cast<double>(dataSearch[*ppr]) > dValue) {
							fExtrema = false;
							break;
						}
					}

					int iHasMore = kd_res_next(kdresMerge);
					if (!iHasMore) {
						break;
					}
				}

				if (fExtrema) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedMerge++;
				}
			}

			kd_res_free(kdresMerge);
		}

		// Destroy the KD Tree
		kd_free(kdMerge);

		// Update set of pressure minima
		setCandidates = setNewCandidates;
	}

	// Eliminate based on thresholds
	for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

		std::set<int> setNewCandidates;

		// Load the search variable data
		vecThresholdOp[tc].m_pvar->LoadGridData(pobjConfig, sTime);

		const DataArray1D<float> & dataState =
			vecThresholdOp[tc].m_pvar->GetData();

		// Loop through all pressure minima
		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();

		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Determine if the threshold is satisfied
			bool fSatisfiesThreshold =
				SatisfiesThreshold<float>(
					mesh,
					dataState,
					*iterCandidate,
					vecThresholdOp[tc].m_eOp,
					vecThresholdOp[tc].m_dValue,
					vecThresholdOp[tc].m_dDistance
				);

			// If not rejected, add to new pressure minima array
			if (fSatisfiesThreshold) {
				setNewCandidates.insert(*iterCandidate);
			} else {
				vecRejectedThreshold[tc]++;
			}
		}

		setCandidates = setNewCandidates;
	}

	// Eliminate based on closed contours
	for (int ccc = 0; ccc < vecClosedContourOp.size(); ccc++) {
		std::set<int> setNewCandidates;

		// Load the search variable data
		vecClosedContourOp[ccc].m_pvar->LoadGridData(pobjConfig, sTime);

		const DataArray1D<float> & dataState =
			vecClosedContourOp[ccc].m_pvar->GetData();

		// Loop through all pressure minima
		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();

		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Determine if a closed contour is present
			bool fHasClosedContour =
				HasClosedContour<float>(
					mesh,
					dataState,
					*iterCandidate,
					vecClosedContourOp[ccc].m_dDeltaAmount,
					vecClosedContourOp[ccc].m_dDistance,
					vecClosedContourOp[ccc].m_dMinMaxDist
				);

			// If not rejected, add to new pressure minima array
			if (fHasClosedContour) {
				setNewCandidates.insert(*iterCandidate);
			} else {
				vecRejectedClosedContour[ccc]++;
			}
		}

		setCandidates = setNewCandidates;
	}

	// Eliminate based on no closed contours
	for (int ccc = 0; ccc < vecNoClosedContourOp.size(); ccc++) {
		std::set<int> setNewCandidates;

		// Load the search variable data
		vecNoClosedContourOp[ccc].m_pvar->LoadGridData(pobjConfig, sTime);

		const DataArray1D<float> & dataState =
			vecNoClosedContourOp[ccc].m_pvar->GetData();

		// Loop through all pressure minima
		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();

		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Determine if a closed contour is present
			bool fHasClosedContour =
				HasClosedContour<float>(
					mesh,
					dataState,
					*iterCandidate,
					vecNoClosedContourOp[ccc].m_dDeltaAmount,
					vecNoClosedContourOp[ccc].m_dDistance,
					vecNoClosedContourOp[ccc].m_dMinMaxDist
				);

			// If a closed contour is present, reject this candidate
			if (fHasClosedContour) {
				vecRejectedNoClosedContour[ccc]++;
			} else {
				setNewCandidates.insert(*iterCandidate);
			}
		}

		setCandidates = setNewCandidates;
	}

	Announce("Total candidates: %i", setCandidates.size());
	Announce("Rejected (  location): %i", nRejectedLocation);
	Announce("Rejected (topography): %i", nRejectedTopography);
	Announce("Rejected (    merged): %i", nRejectedMerge);

	for (int tc = 0; tc < vecRejectedThreshold.size(); tc++) {
		Announce("Rejected (thresh. %s): %i",
				vecThresholdOp[tc].m_pvar->Name().c_str(),
				vecRejectedThreshold[tc]);
	}

	for (int ccc = 0; ccc < vecRejectedClosedContour.size(); ccc++) {
		Announce("Rejected (contour %s): %i",
				vecClosedContourOp[ccc].m_pvar->Name().c_str(),
				vecRejectedClosedContour[ccc]);
	}

	for (int ccc = 0; ccc < vecRejectedNoClosedContour.size(); ccc++) {
		Announce("Rejected (nocontour %s): %i",
				vecNoClosedContourOp[ccc].m_pvar->Name().c_str(),
				vecRejectedNoClosedContour[ccc]);
	}

	// Write results to PointDataObject
	if (pobjPointData != NULL) {

		// Resize the results structure
		pobjPointData->SetFieldCount(2, vecOutputOp.size(), 0);
		pobjPointData->Resize(setCandidates.size());

		// Add time and location for all candidates
		{
			// Add all candidate data to PointDataObject
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			int iCandidateCount = 0;
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				pobjPointData->DataInt(iCandidateCount,0) = sTime;
				pobjPointData->DataInt(iCandidateCount,1) = *iterCandidate;

				iCandidateCount++;
			}
		}

		// Apply output operators
		for (int outc = 0; outc < vecOutputOp.size(); outc++) {

			// Load the search variable data
			vecOutputOp[outc].m_pvar->LoadGridData(pobjConfig, sTime);

			const DataArray1D<float> & dataState =
				vecOutputOp[outc].m_pvar->GetData();

			// Add all candidate data to PointDataObject
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			int iCandidateCount = 0;
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				int ixExtremum;
				float dValue;
				float dRMax;

				if (vecOutputOp[outc].m_eOp == OutputOp::Max) {
					FindLocalMinMax<float>(
						mesh,
						false,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					pobjPointData->DataFloat(iCandidateCount,outc) = dValue;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::MaxDist) {
					FindLocalMinMax<float>(
						mesh,
						false,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					pobjPointData->DataFloat(iCandidateCount,outc) = dRMax;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::Min) {
					FindLocalMinMax<float>(
						mesh,
						true,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					pobjPointData->DataFloat(iCandidateCount,outc) = dValue;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::MinDist) {
					FindLocalMinMax<float>(
						mesh,
						true,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					pobjPointData->DataFloat(iCandidateCount,outc) = dRMax;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::Avg) {
					FindLocalAverage<float>(
						mesh,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						dValue);

					pobjPointData->DataFloat(iCandidateCount,outc) = dValue;

				} else {
					_EXCEPTIONT("Invalid Output operator");
				}

				iCandidateCount++;
			}
		}
/*
		// Output all candidates
		iCandidateCount = 0;

		std::set<int>::const_iterator iterCandidate = setCandidates.begin();
		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			if (grid.m_nGridDim.size() == 1) {
				fprintf(fpOutput, "\t%i", *iterCandidate);

			} else if (grid.m_nGridDim.size() == 2) {
				fprintf(fpOutput, "\t%i\t%i",
					(*iterCandidate) % static_cast<int>(grid.m_nGridDim[1]),
					(*iterCandidate) / static_cast<int>(grid.m_nGridDim[1]));
			}

			fprintf(fpOutput, "\t%3.6f\t%3.6f",
				mesh.dLon[*iterCandidate] * 180.0 / M_PI,
				mesh.dLat[*iterCandidate] * 180.0 / M_PI);

			for (int outc = 0; outc < vecOutputOp.size(); outc++) {
				fprintf(fpOutput, "\t%3.6e",
					dOutput[iCandidateCount][outc]);
			}

			fprintf(fpOutput, "\n");

			iCandidateCount++;
		}
*/
	}

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string PointSearchFunction::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	AnnounceStartBlock("BEGIN point_search");

	if (vecCommandLine.size() != 2) {
		return std::string("point_search(config, param) expects 2 arguments");
	}

	// Check that first argument is of type RecapConfigObject
	RecapConfigObject * pobjConfig =
		dynamic_cast<RecapConfigObject *>(objreg.GetObject(vecCommandLine[0]));
	if (pobjConfig == NULL) {
		return std::string("First argument to point_search must be of type "
			"recap_configuration");
	}

	// Check for valid configuration
	if (!pobjConfig->IsValid()) {
		return std::string("First argument is an invalid recap_configuration "
			"object: Missing \"data\" or \"grid\" properties");
	}

	// Construct the adjacency list (needed for this routine)
	Mesh & mesh = pobjConfig->GetGrid()->GetMesh();

	mesh.ConstructAdjacencyList();

	// Parameters
	Object * pobjParam =
		dynamic_cast<Object *>(objreg.GetObject(vecCommandLine[1]));
	if (pobjParam == NULL) {
		return std::string("Second argument to point_search must be of type parameter_list");
	}

	// Construct PointSearchParam object
	PointSearchParam dcuparam;

	StringObject * pobjSearchByMin =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("searchbymin")));

	StringObject * pobjSearchByMax =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("searchbymax")));

	if (((pobjSearchByMin == NULL) && (pobjSearchByMax == NULL)) ||
	    ((pobjSearchByMin != NULL) && (pobjSearchByMax != NULL))
	) {
		return std::string("Exactly one of searchbymin or searchbymax "
			"must be specified in param");
	}

	if (pobjSearchByMin != NULL) {
		dcuparam.fSearchByMinima = true;
		dcuparam.strVariableSearchBy = pobjSearchByMin->Value();
	}
	if (pobjSearchByMax != NULL) {
		dcuparam.fSearchByMinima = false;
		dcuparam.strVariableSearchBy = pobjSearchByMax->Value();
	}

	// Maximum latitude for detection
	StringObject * pobjMaximumLatitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("maxlat")));
	if (pobjMaximumLatitude != NULL) {
		bool fSuccess =
			pobjMaximumLatitude->ToUnit(
				"deg", &(dcuparam.dMaxLatitude));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("maxlat")
				+ std::string(" to \"deg\"");
		}
	}

	// Minimum latitude for detection
	StringObject * pobjMinimumLatitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("minlat")));
	if (pobjMinimumLatitude != NULL) {
		bool fSuccessMag =
			pobjMinimumLatitude->ToUnit(
				"deg", &(dcuparam.dMinLatitude));
		if (!fSuccessMag) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("minlat")
				+ std::string(" to \"deg\"");
		}
	}

	// Minimum absolute latitude for detection
	StringObject * pobjMinimumAbsLatitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("minabslat")));
	if (pobjMinimumAbsLatitude != NULL) {
		bool fSuccess =
			pobjMinimumAbsLatitude->ToUnit(
				"deg", &(dcuparam.dMinAbsLatitude));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("minabslat")
				+ std::string(" to \"deg\"");
		}
	}

	// Maximum longitude for detection
	StringObject * pobjMaximumLongitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("maxlon")));
	if (pobjMaximumLongitude != NULL) {
		bool fSuccess =
			pobjMaximumLongitude->ToUnit(
				"deg", &(dcuparam.dMaxLongitude));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("maxlon")
				+ std::string(" to \"deg\"");
		}
	}

	// Minimum longitude for detection
	StringObject * pobjMinimumLongitude =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("minlon")));
	if (pobjMinimumLongitude != NULL) {
		bool fSuccess =
			pobjMinimumLongitude->ToUnit(
				"deg", &(dcuparam.dMinLongitude));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("minlon")
				+ std::string(" to \"deg\"");
		}
	}

	// Merge distance for detection
	StringObject * pobjMergeDist =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("mergedist")));
	if (pobjMergeDist != NULL) {
		bool fSuccess =
			pobjMergeDist->ToUnit(
				"deg", &(dcuparam.dMergeDist));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("mergedist")
				+ std::string(" to \"deg\"");
		}
	}

	// Construct all closed contour objects
	std::vector<ClosedContourOp> vecClosedContourOps;

	ListObject * pobjClosedContourList =
		dynamic_cast<ListObject *>(
			objreg.GetObject(pobjParam->ChildName("closedcontourcmd")));

	if (pobjClosedContourList != NULL) {
		for (size_t i = 0; i < pobjClosedContourList->ChildrenCount(); i++) {
			ClosedContourOp opCC;

			Object * pobjChild = pobjClosedContourList->GetChild(i);

			std::string strError =
				opCC.InitializeFromObject(
					objreg, pobjConfig, pobjChild);

			if (strError != "") {
				return strError;
			}

			vecClosedContourOps.push_back(opCC);

			Announce(opCC.ToString().c_str());
		}
	}

	dcuparam.pvecClosedContourOp = &vecClosedContourOps;

	// Construct all no closed contour objects
	std::vector<ClosedContourOp> vecNoClosedContourOps;

	ListObject * pobjNoClosedContourList =
		dynamic_cast<ListObject *>(
			objreg.GetObject(pobjParam->ChildName("noclosedcontourcmd")));

	if (pobjNoClosedContourList != NULL) {
		for (size_t i = 0; i < pobjNoClosedContourList->ChildrenCount(); i++) {
			ClosedContourOp opCC;

			Object * pobjChild = pobjNoClosedContourList->GetChild(i);

			std::string strError =
				opCC.InitializeFromObject(
					objreg, pobjConfig, pobjChild);

			if (strError != "") {
				return strError;
			}

			vecNoClosedContourOps.push_back(opCC);

			Announce(opCC.ToString().c_str());
		}
	}

	dcuparam.pvecNoClosedContourOp = &vecNoClosedContourOps;

	// Threshold operators
	std::vector<ThresholdOp> vecThresholdOps;

	ListObject * pobjThresholdCmdList =
		dynamic_cast<ListObject *>(
			objreg.GetObject(pobjParam->ChildName("thresholdcmd")));

	if (pobjThresholdCmdList != NULL) {
		for (size_t i = 0; i < pobjThresholdCmdList->ChildrenCount(); i++) {
			ThresholdOp opThreshold;

			Object * pobjChild = pobjNoClosedContourList->GetChild(i);

			std::string strError =
				opThreshold.InitializeFromObject(
					objreg, pobjConfig, pobjChild);

			if (strError != "") {
				return strError;
			}

			vecThresholdOps.push_back(opThreshold);

			Announce(opThreshold.ToString().c_str());
		}
	}

	dcuparam.pvecThresholdOp = &vecThresholdOps;

	// Output operators
	std::vector<OutputOp> vecOutputOps;

	ListObject * pobjOutputCmdList =
		dynamic_cast<ListObject *>(
			objreg.GetObject(pobjParam->ChildName("output")));

	if (pobjOutputCmdList != NULL) {
		for (size_t i = 0; i < pobjOutputCmdList->ChildrenCount(); i++) {
			OutputOp opOutput;

			Object * pobjChild = pobjOutputCmdList->GetChild(i);

			std::string strError =
				opOutput.InitializeFromObject(
					objreg, pobjConfig, pobjChild);

			if (strError != "") {
				return strError;
			}

			vecOutputOps.push_back(opOutput);
		}
	}

	dcuparam.pvecOutputOp = &vecOutputOps;

	// Distribute time steps over ranks and run the point search
	dcuparam.fpLog = stdout;

	std::vector<size_t> vecTimeIndices;
	pobjConfig->GetFileList()->GetOnRankTimeIndices(vecTimeIndices);

	std::vector<PointDataObject *> vecpobjPointData;
	vecpobjPointData.resize(vecTimeIndices.size(), NULL);

	for (size_t t = 0; t < vecTimeIndices.size(); t++) {
		if (ppReturn != NULL) {
			vecpobjPointData[t] = new PointDataObject("");
		}

		std::string strError =
			PointSearch(
				pobjConfig,
				vecTimeIndices[t],
				dcuparam,
				vecpobjPointData[t]);

		if (strError != "") {
			return strError;
		}
	}

	// Check if a return value is needed
	if (ppReturn != NULL) {

		// Concatenate all PointDataObjects into a single return value
		PointDataObject * pobjPointDataCombined = new PointDataObject("");

		pobjPointDataCombined->Concatenate(vecpobjPointData);

		if (pobjPointDataCombined->GetIntFieldCount() != 2) {
			_EXCEPTIONT("Logic error");
		}
		if (pobjPointDataCombined->GetFloatFieldCount() != vecOutputOps.size()) {
			_EXCEPTIONT("Logic error");
		}

		// Set field headers
		pobjPointDataCombined->SetIntFieldHeader(0, "time_ix");
		pobjPointDataCombined->SetIntFieldHeader(1, "space_ix");

		for (size_t i = 0; i < vecOutputOps.size(); i++) {
			pobjPointDataCombined->SetFloatFieldHeader(
				i, vecOutputOps[i].m_pvar->Name());

			pobjPointDataCombined->SetFloatFieldUnits(
				i, vecOutputOps[i].m_pvar->Units());
		}

		// Cleanup
		for (size_t s = 0; s < vecpobjPointData.size(); s++) {
			delete vecpobjPointData[s];
		}

		(*ppReturn) = pobjPointDataCombined;

		Announce("Total candidates over all timesteps: %i",
			pobjPointDataCombined->GetRows());
/*
		Announce("%lu int fields, %lu float fields",
			pobjPointDataCombined->GetIntFieldCount(),
			pobjPointDataCombined->GetFloatFieldCount());
*/
	}

	AnnounceEndBlock("Done");
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

}

