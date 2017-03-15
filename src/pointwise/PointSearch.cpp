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

#if defined(TEMPEST_MPIOMP)
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
	///		Convert to string.
	///	</summary>
	std::string ToString() const {
/*
		// Output announcement
		std::string strDescription = var.ToString(varreg);
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

		char szBuffer[128];
		sprintf(szBuffer, "%f within %f degrees",
			m_dValue, m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
*/
		return std::string("");
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
/*
		if (m_dDeltaAmount < 0.0) {
			Announce("%s decreases by %f over %f degrees"
				   " (max search %f deg)",
				var.ToString(varreg).c_str(),
				-m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);

		} else {
			Announce("%s increases by %f over %f degrees"
					" (min search %f deg)",
				var.ToString(varreg).c_str(),
				m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);
		}
*/
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
	///		Convert to string.
	///	</summary>
	std::string ToString() const {
/*
		sprintf(szBuffer, "%s", var.ToString(varreg).c_str());
		strDescription += szBuffer;

		sprintf(szBuffer, " within %f degrees", m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
*/
		return std::string("");
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

	// Latitude and longitude at the origin
	double dLat0 = mesh.dLat[ix0];
	double dLon0 = mesh.dLon[ix0];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		double dLatThis = mesh.dLat[ix];
		double dLonThis = mesh.dLon[ix];

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

	const double dLat0 = mesh.dLat[ixOrigin];
	const double dLon0 = mesh.dLon[ixOrigin];

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
		double dLatThis = mesh.dLat[ix];
		double dLonThis = mesh.dLon[ix];

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
	double dLat0 = mesh.dLat[ix0];
	double dLon0 = mesh.dLon[ix0];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element
		double dLatThis = mesh.dLat[ix];
		double dLonThis = mesh.dLon[ix];

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

	// Latitude and longitude at the origin
	double dLat0 = mesh.dLat[ix0];
	double dLon0 = mesh.dLon[ix0];

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

		double dLatThis = mesh.dLat[ix];
		double dLonThis = mesh.dLon[ix];

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
		pvarSearchBy(NULL),
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
	Variable * pvarSearchBy;

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

void PointSearch(
	size_t iFileIx,
	size_t iTimeIx,
	const Mesh & mesh,
	RecapConfigObject * pobjConfig,
	const PointSearchParam & param
) {
	// Set the Announce buffer
	if (param.fpLog == NULL) {
		_EXCEPTIONT("Invalid log buffer");
	}

	AnnounceSetOutputBuffer(param.fpLog);
	AnnounceOutputOnAllRanks();

	// Check minimum longitude / latitude
	if ((param.dMinLongitude < 0.0) || (param.dMinLongitude >= 360.0)) {
		_EXCEPTIONT("Invalid MinLongitude");
	}
	if ((param.dMaxLongitude < 0.0) || (param.dMaxLongitude >= 360.0)) {
		_EXCEPTIONT("Invalid MaxLongitude");
	}
	if ((param.dMaxLatitude < -90.0) || (param.dMaxLatitude > 90.0)) {
		_EXCEPTIONT("--maxlat must in the range [-90,90]");
	}
	if ((param.dMinLatitude < -90.0) || (param.dMinLatitude > 90.0)) {
		_EXCEPTIONT("--minlat must in the range [-90,90]");
	}
	if ((param.dMinAbsLatitude < 0.0) || (param.dMinAbsLatitude > 90.0)) {
		_EXCEPTIONT("--minabslat must in the range [0,90]");
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

	// Load the data for the search variable
	if (param.pvarSearchBy == NULL) {
		_EXCEPTION();
	}
	const DataArray1D<float> & dataSearch = param.pvarSearchBy->GetData();

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

	DataArray1D<int> vecRejectedClosedContour(vecClosedContourOp.size());
	DataArray1D<int> vecRejectedNoClosedContour(vecNoClosedContourOp.size());
	DataArray1D<int> vecRejectedThreshold(vecThresholdOp.size());

	// Eliminate based on interval
	if ((param.dMinLatitude != param.dMaxLatitude) ||
	    (param.dMinLongitude != param.dMaxLongitude) ||
		(param.dMinAbsLatitude != 0.0)
	) {
		std::set<int> setNewCandidates;

		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();
		for (; iterCandidate != setCandidates.end(); iterCandidate++) {
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
					int iLonShift = static_cast<int>(dLon / (2.0 * M_PI));
					dLon += static_cast<double>(iLonShift + 1) * 2.0 * M_PI;
				}
				if (dLon >= 2.0 * M_PI) {
					int iLonShift = static_cast<int>(dLon / (2.0 * M_PI));
					dLon -= static_cast<double>(iLonShift - 1) * 2.0 * M_PI;
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
/*
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
			double dLat = mesh.dLat[*iterCandidate];
			double dLon = mesh.dLon[*iterCandidate];

			double dX = cos(dLon) * cos(dLat);
			double dY = sin(dLon) * cos(dLat);
			double dZ = sin(dLat);

			kd_insert3(kdMerge, dX, dY, dZ, (void*)(&(*iterCandidate)));
		}

		// Loop through all candidates find set of nearest neighbors
		iterCandidate = setCandidates.begin();
		for (; iterCandidate != setCandidates.end(); iterCandidate++) {
			double dLat = mesh.dLat[*iterCandidate];
			double dLon = mesh.dLon[*iterCandidate];

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
		Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);
		var.LoadGridData(varreg, vecFiles, grid, t);
		const DataArray1D<float> & dataState = var.GetData();

		// Loop through all pressure minima
		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();

		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Determine if the threshold is satisfied
			bool fSatisfiesThreshold =
				SatisfiesThreshold<float>(
					grid,
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
		Variable & var = varreg.Get(vecClosedContourOp[ccc].m_varix);
		var.LoadGridData(varreg, vecFiles, grid, t);
		const DataArray1D<float> & dataState = var.GetData();

		// Loop through all pressure minima
		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();

		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Determine if a closed contour is present
			bool fHasClosedContour =
				HasClosedContour<float>(
					grid,
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
		Variable & var = varreg.Get(vecNoClosedContourOp[ccc].m_varix);
		var.LoadGridData(varreg, vecFiles, grid, t);
		const DataArray1D<float> & dataState = var.GetData();

		// Loop through all pressure minima
		std::set<int>::const_iterator iterCandidate
			= setCandidates.begin();

		for (; iterCandidate != setCandidates.end(); iterCandidate++) {

			// Determine if a closed contour is present
			bool fHasClosedContour =
				HasClosedContour<float>(
					grid,
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

	for (int tc = 0; tc < vecRejectedThreshold.GetRows(); tc++) {
		Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);

		Announce("Rejected (thresh. %s): %i",
				var.m_strName.c_str(),
				vecRejectedThreshold[tc]);
	}

	for (int ccc = 0; ccc < vecRejectedClosedContour.GetRows(); ccc++) {
		Variable & var = varreg.Get(vecClosedContourOp[ccc].m_varix);

		Announce("Rejected (contour %s): %i",
				var.m_strName.c_str(),
				vecRejectedClosedContour[ccc]);
	}

	for (int ccc = 0; ccc < vecRejectedNoClosedContour.GetRows(); ccc++) {
		Variable & var = varreg.Get(vecNoClosedContourOp[ccc].m_varix);

		Announce("Rejected (nocontour %s): %i",
				var.m_strName.c_str(),
				vecRejectedNoClosedContour[ccc]);
	}

	// Write results to file
	{
		// Parse time information
		//NcVar * varDate = ncInput.get_var("date");
		//NcVar * varDateSec = ncInput.get_var("datesec");

		int nDateYear;
		int nDateMonth;
		int nDateDay;
		int nDateHour;

		NcAtt * attTimeUnits = varTime->get_att("units");
		if (attTimeUnits == NULL) {
			_EXCEPTIONT("Variable \"time\" has no \"units\" attribute");
		}

		std::string strTimeUnits = attTimeUnits->as_string(0);

		std::string strTimeCalendar = "standard";
		NcAtt * attTimeCalendar = varTime->get_att("calendar");
		if (attTimeCalendar != NULL) {
			strTimeCalendar = attTimeCalendar->as_string(0);
		}

		ParseTimeDouble(
			strTimeUnits,
			strTimeCalendar,
			dTime[t],
			nDateYear,
			nDateMonth,
			nDateDay,
			nDateHour);

		// Write time information
		fprintf(fpOutput, "%i\t%i\t%i\t%i\t%i\n",
			nDateYear,
			nDateMonth,
			nDateDay,
			static_cast<int>(setCandidates.size()),
			nDateHour);

		// Write candidate information
		int iCandidateCount = 0;

		// Apply output operators
		DataArray2D<float> dOutput(setCandidates.size(), vecOutputOp.size());
		for (int outc = 0; outc < vecOutputOp.size(); outc++) {

			// Load the search variable data
			Variable & var = varreg.Get(vecOutputOp[outc].m_varix);
			var.LoadGridData(varreg, vecFiles, grid, t);
			const DataArray1D<float> & dataState = var.GetData();

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			iCandidateCount = 0;
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				int ixExtremum;
				float dValue;
				float dRMax;

				if (vecOutputOp[outc].m_eOp == OutputOp::Max) {
					FindLocalMinMax<float>(
						grid,
						false,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					dOutput[iCandidateCount][outc] = dValue;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::MaxDist) {
					FindLocalMinMax<float>(
						grid,
						false,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					dOutput[iCandidateCount][outc] = dRMax;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::Min) {
					FindLocalMinMax<float>(
						grid,
						true,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					dOutput[iCandidateCount][outc] = dValue;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::MinDist) {
					FindLocalMinMax<float>(
						grid,
						true,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						ixExtremum,
						dValue,
						dRMax);

					dOutput[iCandidateCount][outc] = dRMax;

				} else if (vecOutputOp[outc].m_eOp == OutputOp::Avg) {
					FindLocalAverage<float>(
						grid,
						dataState,
						*iterCandidate,
						vecOutputOp[outc].m_dDistance,
						dValue);

					dOutput[iCandidateCount][outc] = dValue;

				} else {
					_EXCEPTIONT("Invalid Output operator");
				}

				iCandidateCount++;
			}
		}

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
	}
*/
	AnnounceEndBlock("Done");
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

	// Check validity of configuration
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
/*
	// Grid object
	GridObject * pobjGrid =
		dynamic_cast<GridObject *>(pobjConfig->GetChild("data"));

	// Grid object
	GridObject * pobjGrid =
		dynamic_cast<GridObject *>(objreg.GetObject(vecCommandLine[0]));
	if (pobjGrid == NULL) {
		return std::string("First argument to point_search must be of type grid");
	}

	// FileList object
	FileListObject * pobjFileList =
		dynamic_cast<FileListObject *>(objreg.GetObject(vecCommandLine[1]));
	if (pobjFileList == NULL) {
		return std::string("Second argument to point_search must be of type file_list");
	}

	// Lookup object
	VariableLookupObject * pobjVarLookup =
		dynamic_cast<VariableLookupObject *>(objreg.GetObject(vecCommandLine[2]));
	if (pobjVarLookup == NULL) {
		return std::string("Third argument to point_search must be of type variable_lookup");
	}
*/
	// Parameters
	Object * pobjParam =
		dynamic_cast<Object *>(objreg.GetObject(vecCommandLine[1]));
	if (pobjParam == NULL) {
		return std::string("Fourth argument to point_search must be of type parameter_list");
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
		}
	}

	dcuparam.pvecNoClosedContourOp = &vecNoClosedContourOps;

	// Distribute time steps over ranks

	// Check if a return value is needed
	if (ppReturn != NULL) {
		(*ppReturn) = new StringObject("NULL", "NULL");
	}

	AnnounceEndBlock("Done");
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

/*
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {
	// Parameters for DetectCycloneUnstructured
	PointSearchParam ptsparam;

	// Input dat file
	std::string strInputFile;

	// Input list file
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutput;

	// Output file list
	std::string strOutputFileList;

	// Variable to search for the minimum
	std::string strSearchByMin;

	// Variable to search for the maximum
	std::string strSearchByMax;

	// Closed contour commands
	std::string strClosedContourCmd;

	// Closed contour commands
	std::string strNoClosedContourCmd;

	// Threshold commands
	std::string strThresholdCmd;

	// Output commands
	std::string strOutputCmd;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutput, "out", "");
		CommandLineString(strOutputFileList, "out_file_list", "");
		CommandLineStringD(strSearchByMin, "searchbymin", "", "(default PSL)");
		CommandLineString(strSearchByMax, "searchbymax", "");
		CommandLineDoubleD(ptsparam.dMinLongitude, "minlon", 0.0, "(degrees)");
		CommandLineDoubleD(ptsparam.dMaxLongitude, "maxlon", 0.0, "(degrees)");
		CommandLineDoubleD(ptsparam.dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineDoubleD(ptsparam.dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(ptsparam.dMinAbsLatitude, "minabslat", 0.0, "(degrees)");
		CommandLineDoubleD(ptsparam.dMergeDist, "mergedist", 0.0, "(degrees)");
		CommandLineStringD(strClosedContourCmd, "closedcontourcmd", "", "[var,delta,dist,minmaxdist;...]");
		CommandLineStringD(strNoClosedContourCmd, "noclosedcontourcmd", "", "[var,delta,dist,minmaxdist;...]");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		CommandLineStringD(strOutputCmd, "outputcmd", "", "[var,op,dist;...]");
		CommandLineInt(ptsparam.nTimeStride, "timestride", 1);
		CommandLineBool(ptsparam.fRegional, "regional");
		CommandLineBool(ptsparam.fOutputHeader, "out_header");
		CommandLineInt(ptsparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(ptsparam.iVerbosityLevel);

	// Check input
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}

	// Check output
	if ((strOutput.length() != 0) && (strOutputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out) or (--out_data_list)"
			" may be specified");
	}

	// Load input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		std::ifstream ifInputFileList(strInputFileList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputFiles.push_back(strFileLine);
		}
	}

	// Load output file list
	std::vector<std::string> vecOutputFiles;

	if (strOutputFileList.length() != 0) {

		std::ifstream ifOutputFileList(strOutputFileList.c_str());
		if (!ifOutputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strOutputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecOutputFiles.push_back(strFileLine);
		}

		if (vecOutputFiles.size() != vecInputFiles.size()) {
			_EXCEPTIONT("File --in_file_list must match --out_file_list");
		}
	}

	// Only one of search by min or search by max should be specified
	if ((strSearchByMin == "") && (strSearchByMax == "")) {
		strSearchByMin = "PSL";
	}
	if ((strSearchByMin != "") && (strSearchByMax != "")) {
		_EXCEPTIONT("Only one of --searchbymin or --searchbymax can"
			" be specified");
	}

	ptsparam.fSearchByMinima = false;
	{
		Variable varSearchByArg;
		if (strSearchByMin != "") {
			varSearchByArg.ParseFromString(varreg, strSearchByMin);
			ptsparam.fSearchByMinima = true;
		}
		if (strSearchByMax != "") {
			varSearchByArg.ParseFromString(varreg, strSearchByMax);
			ptsparam.fSearchByMinima = false;
		}

		ptsparam.ixSearchBy = varreg.FindOrRegister(varSearchByArg);
	}

	// Parse the closed contour command string
	std::vector<ClosedContourOp> vecClosedContourOp;
	ptsparam.pvecClosedContourOp = &vecClosedContourOp;

	if (strClosedContourCmd != "") {
		AnnounceStartBlock("Parsing closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strClosedContourCmd.length(); i++) {

			if ((i == strClosedContourCmd.length()) ||
				(strClosedContourCmd[i] == ';') ||
				(strClosedContourCmd[i] == ':')
			) {
				std::string strSubStr =
					strClosedContourCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecClosedContourOp.size());
				vecClosedContourOp.resize(iNextOp + 1);
				vecClosedContourOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the no closed contour command string
	std::vector<ClosedContourOp> vecNoClosedContourOp;
	ptsparam.pvecNoClosedContourOp = &vecNoClosedContourOp;

	if (strNoClosedContourCmd != "") {
		AnnounceStartBlock("Parsing no closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strNoClosedContourCmd.length(); i++) {

			if ((i == strNoClosedContourCmd.length()) ||
				(strNoClosedContourCmd[i] == ';') ||
				(strNoClosedContourCmd[i] == ':')
			) {
				std::string strSubStr =
					strNoClosedContourCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecNoClosedContourOp.size());
				vecNoClosedContourOp.resize(iNextOp + 1);
				vecNoClosedContourOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the threshold operator command string
	std::vector<ThresholdOp> vecThresholdOp;
	ptsparam.pvecThresholdOp = &vecThresholdOp;

	if (strThresholdCmd != "") {
		AnnounceStartBlock("Parsing threshold operations");

		int iLast = 0;
		for (int i = 0; i <= strThresholdCmd.length(); i++) {

			if ((i == strThresholdCmd.length()) ||
				(strThresholdCmd[i] == ';') ||
				(strThresholdCmd[i] == ':')
			) {
				std::string strSubStr =
					strThresholdCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecThresholdOp.size());
				vecThresholdOp.resize(iNextOp + 1);
				vecThresholdOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the output operator command string
	std::vector<OutputOp> vecOutputOp;
	ptsparam.pvecOutputOp = &vecOutputOp;

	if (strOutputCmd != "") {
		AnnounceStartBlock("Parsing output operations");

		int iLast = 0;
		for (int i = 0; i <= strOutputCmd.length(); i++) {

			if ((i == strOutputCmd.length()) ||
				(strOutputCmd[i] == ';') ||
				(strOutputCmd[i] == ':')
			) {
				std::string strSubStr =
					strOutputCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecOutputOp.size());
				vecOutputOp.resize(iNextOp + 1);
				vecOutputOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Check minimum/maximum latitude/longitude
	if ((ptsparam.dMaxLatitude < -90.0) || (ptsparam.dMaxLatitude > 90.0)) {
		_EXCEPTIONT("--maxlat must in the range [-90,90]");
	}
	if ((ptsparam.dMinLatitude < -90.0) || (ptsparam.dMinLatitude > 90.0)) {
		_EXCEPTIONT("--minlat must in the range [-90,90]");
	}
	if (ptsparam.dMinLatitude > ptsparam.dMaxLatitude) {
		_EXCEPTIONT("--minlat must be less than --maxlat");
	}

	ptsparam.dMaxLatitude *= M_PI / 180.0;
	ptsparam.dMinLatitude *= M_PI / 180.0;

	if (ptsparam.dMinLongitude < 0.0) {
		int iMinLongitude =
			static_cast<int>(-ptsparam.dMinLongitude / 360.0);
		ptsparam.dMinLongitude +=
			static_cast<double>(iMinLongitude + 1) * 360.0;
	}
	if (ptsparam.dMinLongitude >= 360.0) {
		int iMinLongitude =
			static_cast<int>(ptsparam.dMinLongitude / 360.0);
		ptsparam.dMinLongitude -=
			static_cast<double>(iMinLongitude - 1) * 360.0;
	}
	if (ptsparam.dMaxLongitude < 0.0) {
		int iMaxLongitude =
			static_cast<int>(-ptsparam.dMaxLatitude / 360.0);
		ptsparam.dMaxLongitude +=
			static_cast<double>(iMaxLongitude + 1) * 360.0;
	}
	if (ptsparam.dMaxLongitude >= 360.0) {
		int iMaxLongitude =
			static_cast<int>(ptsparam.dMaxLongitude / 360.0);
		ptsparam.dMaxLongitude -=
			static_cast<double>(iMaxLongitude - 1) * 360.0;
	}

	ptsparam.dMaxLongitude *= M_PI / 180.0;
	ptsparam.dMinLongitude *= M_PI / 180.0;

#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
#endif

	AnnounceStartBlock("Begin search operation");
	if (vecInputFiles.size() != 1) {
		if (vecOutputFiles.size() != 0) {
			Announce("Output will be written following --out_file_list");
		} else if (strOutput == "") {
			Announce("Output will be written to outXXXXXX.dat");
		} else {
			Announce("Output will be written to %sXXXXXX.dat",
				strOutput.c_str());
		}
		Announce("Logs will be written to logXXXXXX.txt");
	}

	// Loop over all files to be processed
	for (int f = 0; f < vecInputFiles.size(); f++) {
#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}
#endif
		// Generate output file name
		std::string strOutputFile;
		if (vecInputFiles.size() == 1) {
			ptsparam.fpLog = stdout;

			if (strOutput == "") {
				strOutputFile = "out.dat";
			} else {
				strOutputFile = strOutput;
			}

		} else {
			char szFileIndex[32];
			sprintf(szFileIndex, "%06i", f);

			if (vecOutputFiles.size() != 0) {
				strOutputFile = vecOutputFiles[f];
			} else {
				if (strOutput == "") {
					strOutputFile =
						"out" + std::string(szFileIndex) + ".dat";
				} else {
					strOutputFile =
						strOutput + std::string(szFileIndex) + ".dat";
				}
			}

			std::string strLogFile = "log" + std::string(szFileIndex) + ".txt";
			ptsparam.fpLog = fopen(strLogFile.c_str(), "w");
		}

		// Perform PointSearch
		PointSearch(
			f,
			vecInputFiles[f],
			strOutputFile,
			strConnectivity,
			varreg,
			ptsparam);

		// Close the log file
		if (vecInputFiles.size() != 1) {
			fclose(ptsparam.fpLog);
		}
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}
*/

}

