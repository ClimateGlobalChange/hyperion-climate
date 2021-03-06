///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchNodes.cpp
///	\author  Paul Ullrich
///	\version March 16, 2017
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

#include "StitchNodes.h"
#include "RecapConfigObject.h"
#include "PointDataObject.h"
#include "PathDataObject.h"

#include "kdtree.h"

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

namespace HPointwise {

///////////////////////////////////////////////////////////////////////////////

typedef std::vector< std::vector<std::string> > TimesVector;

typedef std::vector< std::vector< std::vector<std::string> > > CandidateVector;

///////////////////////////////////////////////////////////////////////////////

class PathSegment {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PathSegment(
		int iTime0,
		int iCandidate0,
		int iTime1,
		int iCandidate1
	) {
		m_iTime[0] = iTime0;
		m_iTime[1] = iTime1;
		
		m_iCandidate[0] = iCandidate0;
		m_iCandidate[1] = iCandidate1;
	}

public:
	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator< (const PathSegment & seg) const {
		std::pair<int,int> pr0(m_iTime[0], m_iCandidate[0]);
		std::pair<int,int> pr1(seg.m_iTime[0], seg.m_iCandidate[0]);
		return (pr0 < pr1);
	}

public:
	///	<summary>
	///		Begin and end time.
	///	</summary>
	int m_iTime[2];

	///	<summary>
	///		Begin and end candidate.
	///	</summary>
	int m_iCandidate[2];
};

typedef std::set<PathSegment> PathSegmentSet;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing the trajectory of a pointwise feature.
///	</summary>
class Path {

public:
	///	<summary>
	///		Array of times.
	///	</summary>
	std::vector<int> m_iTimes;

	///	<summary>
	///		Array of candidates.
	///	</summary>
	std::vector<int> m_iCandidates;
};

///	<summary>
///		A vector of Path objects.
///	</summary>
typedef std::vector<Path> PathVector;

///////////////////////////////////////////////////////////////////////////////

class PathThresholdOp {

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
		NotEqualTo,
		AbsGreaterThanEqualTo,
		AbsLessThanEqualTo
	};

public:
	///	<summary>
	///		Initialize from an Object.
	///	</summary>
	std::string InitializeFromObject(
		const ObjectRegistry & objreg,
		RecapConfigObject * pobjConfig,
		const PointDataObject * pobjPointData,
		const Object * pobj
	) {
		// Column
		StringObject * pobjColumnName =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("var")));
		if (pobjColumnName == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"var\" property");
		}

		bool fColFound = false;
		for (size_t i = 0; i < pobjPointData->GetFloatFieldCount(); i++) {
			if (pobjPointData->GetFloatFieldHeader(i) == pobjColumnName->Value()) {
				m_iColumn = i;
				fColFound = true;
				break;
			}
		}
		if (!fColFound) {
			return std::string("ERROR: No column \"")
				+ pobjColumnName->Value()
				+ std::string("\" in point_data ")
				+ pobjPointData->Name();
		}

		std::string strUnits = pobjPointData->GetFloatFieldUnits(m_iColumn);

		// Operator
		StringObject * pobjPathThresholdOp =
			dynamic_cast<StringObject *>(
				objreg.GetObject(pobj->ChildName("op")));
		if (pobjPathThresholdOp == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"op\" property");
		}

		bool fSuccessThresholdOp =
			SetOperatorFromString(pobjPathThresholdOp->Value());
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
				strUnits, &m_dValue);
		if (!fSuccessValue) {
			return std::string("Cannot convert ")
				+ pobj->ChildName("value")
				+ std::string(" to ")
				+ strUnits;
		}

		// Count
		IntegerObject * pobjThresholdCount =
			dynamic_cast<IntegerObject *>(
				objreg.GetObject(pobj->ChildName("count")));
		if (pobjThresholdCount == NULL) {
			return pobj->Name()
				+ std::string(" has invalid \"count\" property");
		}

		m_nMinimumCount = pobjThresholdCount->Value();

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
		} else if (strOp == "|>=") {
			m_eOp = AbsGreaterThanEqualTo;
		} else if (strOp == "|<=") {
			m_eOp = AbsLessThanEqualTo;
		} else {
			return false;
		}
		return true;
	}

public:
	///	<summary>
	///		Verify that the specified path satisfies the threshold op.
	///	</summary>
	bool Apply(
		const Path & path,
		const PointDataObject * pobjPointData,
		const std::vector<int> & vecFirstCandidate
	) {
		int nCount = 0;
		for (int p = 0; p < path.m_iTimes.size(); p++) {
			int t = path.m_iTimes[p];
			int i = path.m_iCandidates[p];
			int s = vecFirstCandidate[t]+i;

			if (s >= pobjPointData->GetRows()) {
				_EXCEPTIONT("Logic error");
			}
			if (m_iColumn >= pobjPointData->GetFloatFieldCount()) {
				_EXCEPTIONT("Logic error");
			}

			float dCandidateValue =
				pobjPointData->DataFloat(s,m_iColumn);

			if ((m_eOp == GreaterThan) &&
				(dCandidateValue > m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == LessThan) &&
				(dCandidateValue < m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == GreaterThanEqualTo) &&
				(dCandidateValue >= m_dValue)
			) {
				nCount++;
			
			} else if (
				(m_eOp == LessThanEqualTo) &&
				(dCandidateValue <= m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == EqualTo) &&
				(dCandidateValue == m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == NotEqualTo) &&
				(dCandidateValue != m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == AbsGreaterThanEqualTo) &&
				(fabs(dCandidateValue) >= m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == AbsLessThanEqualTo) &&
				(fabs(dCandidateValue) <= m_dValue)
			) {
				nCount++;
			}
		}

		// Check that the criteria is satisfied for all segments
		if (m_nMinimumCount == (-1)) {
			if (nCount == (int)(path.m_iTimes.size())) {
				return true;
			} else {
				return false;
			}
		}

		// Check total count against min count
		if (nCount >= m_nMinimumCount) {
			return true;
		} else {
			return false;
		}
	}

protected:
	///	<summary>
	///		Active column in PointDataObject.
	///	</summary>
	size_t m_iColumn;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Minimum number of segments that need to satisfy the op.
	///	</summary>
	int m_nMinimumCount;
};

///////////////////////////////////////////////////////////////////////////////

class StitchNodesParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	StitchNodesParam() :
		dRange(0.0),
		nMinDuration(0),
		dMinEndpointDistance(0.0),
		dMinPathDistance(0.0),
		nMaxGapSize(0),
		nTimeStride(0),
		pvecPathThresholdOp(NULL)
	{ }

public:
	// Range (in degrees)
	double dRange;

	// Minimum path length
	int nMinDuration;

	// Minimum distance between endpoints of path (in degrees)
	double dMinEndpointDistance;

	// Minimum path length (in degrees)
	double dMinPathDistance;

	// Maximum time gap (in time steps)
	int nMaxGapSize;

	// Time step stride
	int nTimeStride;

	// Vector of path headers
	std::string strPathHeaderFormat;

	// Vector of path threshold operators
	std::vector<PathThresholdOp> * pvecPathThresholdOp;
};

///////////////////////////////////////////////////////////////////////////////

std::string StitchNodes(
	RecapConfigObject * pobjConfig,
	const StitchNodesParam & param,
	const PointDataObject * pobjPointData,
	PathDataObject * pobjPathData
) {

	// Get the FileList object
	const FileListObject * pobjFileList = pobjConfig->GetFileList();

	// Number of times
	const size_t sTimes = pobjFileList->GetTimeCount();

	// Time stride from the point data
	const int nTimeStride = pobjPointData->GetTimeStride();

	// Null pointer
	int * noptr = NULL;

	// Path threshold operators
	if (param.pvecPathThresholdOp == NULL) {
		_EXCEPTION();
	}
	std::vector<PathThresholdOp> & vecPathThresholdOp =
		*(param.pvecPathThresholdOp);

	// Get the latitude and longitude indices in PointDataObject
	size_t iLatIndex = (-1);
	size_t iLonIndex = (-1);
	for (size_t i = 0; i < pobjPointData->GetFloatFieldCount(); i++) {
		if (pobjPointData->GetFloatFieldHeader(i) == "lat") {
			iLatIndex = i;
		}
		if (pobjPointData->GetFloatFieldHeader(i) == "lon") {
			iLonIndex = i;
		}
	}
	if (iLatIndex == (-1)) {
		return std::string("ERROR: point_data ")
			+ pobjPointData->Name()
			+ std::string(" missing \"lat\" column, "
				"which is required for stitch_nodes()");
	}
	if (iLonIndex == (-1)) {
		return std::string("ERROR: point_data ")
			+ pobjPointData->Name()
			+ std::string(" missing \"lon\" column, "
				"which is required for stitch_nodes()");
	}

	// Number of candidates at each time
	std::vector<int> vecCandidateCount;
	vecCandidateCount.resize(sTimes);

	// Vector of nodal points for each candidate
	std::vector< std::vector<Node> > vecNodes;
	vecNodes.resize(sTimes);

	// Create kdtree at each time
	AnnounceStartBlock("Creating KD trees at each time level");

	// Vector of first candidates for each time
	std::vector<int> vecFirstCandidate;
	vecFirstCandidate.resize(sTimes);

	// Vector of KD trees
	std::vector<kdtree *> vecKDTrees;
	vecKDTrees.resize(sTimes, NULL);

	for (size_t s = 0; s < pobjPointData->GetRows(); s++) {
		int t = pobjPointData->DataInt(s, 0);

		if ((t < 0) || (t >= sTimes)) {
			_EXCEPTIONT("Logic error");
		}

		vecCandidateCount[t]++;
	}
	for (size_t s = 0; s < pobjPointData->GetRows(); s++) {
		int t = pobjPointData->DataInt(s, 0);

		if (vecKDTrees[t] == NULL) {
			vecKDTrees[t] = kd_create(3);
			vecNodes[t].reserve(vecCandidateCount[t]);
			vecFirstCandidate[t] = s;
		}

		double dLat = pobjPointData->DataFloat(s,iLatIndex);
		double dLon = pobjPointData->DataFloat(s,iLonIndex);

		dLat *= M_PI / 180.0;
		dLon *= M_PI / 180.0;

		double dX = sin(dLon) * cos(dLat);
		double dY = cos(dLon) * cos(dLat);
		double dZ = sin(dLat);

		int iCandidate = vecNodes[t].size();

		vecNodes[t].push_back(Node(dX, dY, dZ));

		kd_insert3(vecKDTrees[t], dX, dY, dZ, reinterpret_cast<void*>(noptr+iCandidate));
	}

	AnnounceEndBlock("Done");

	// Create set of path segments
	AnnounceStartBlock("Populating set of path segments");

	std::vector<PathSegmentSet> vecPathSegmentsSet;
	vecPathSegmentsSet.resize(sTimes-1);

	// Loop through all time levels
	for (int t = 0; t < sTimes-1; t++) {

		// Loop through all points at the current time level
		for (int i = 0; i < vecNodes[t].size(); i++) {

			//if (vecNodes[t+1].size() == 0) {
			//	break;
			//}

			double dX = vecNodes[t][i].x;
			double dY = vecNodes[t][i].y;
			double dZ = vecNodes[t][i].z;

			const int nMaxFutureTime = nTimeStride * (param.nMaxGapSize + 1);
			for (int g = nTimeStride; g <= nMaxFutureTime; g += nTimeStride) {

				if (vecNodes[t+g].size() == 0) {
					continue;
				}
				if (t+g >= sTimes) {
					break;
				}

				if (vecKDTrees[t+g] == NULL) {
					continue;
				}

				kdres * set = kd_nearest3(vecKDTrees[t+g], dX, dY, dZ);

				if (kd_res_size(set) == 0) {
					kd_res_free(set);
					break;
				}

				int iRes =
					  reinterpret_cast<int*>(kd_res_item_data(set))
					- reinterpret_cast<int*>(noptr);

				kd_res_free(set);

				// Great circle distance between points
				double dDX = vecNodes[t+g][iRes].x - dX;
				double dDY = vecNodes[t+g][iRes].y - dY;
				double dDZ = vecNodes[t+g][iRes].z - dZ;

				double dCoord = 0.5 * sqrt(dDX * dDX + dDY * dDY + dDZ * dDZ);
				double dR;
				if (dCoord >= 1.0) {
					dR = 180.0;
				} else {
					dR = 2.0 * asin(dCoord) * 180.0 / M_PI;
				}
				if (dR != dR) {
					_EXCEPTIONT("NaN value detected");
				}

				// Verify great circle distance satisfies range requirement
				if (dR <= param.dRange) {

					// Insert new path segment into set of path segments
					vecPathSegmentsSet[t].insert(
						PathSegment(t, i, t+g, iRes));

					break;
				}
			}
		}
	}

	AnnounceEndBlock("Done");

	// Work forwards to find all paths
	AnnounceStartBlock("Constructing paths");

	std::vector<Path> vecPaths;

	int nRejectedMinLengthPaths = 0;
	int nRejectedMinEndpointDistPaths = 0;
	int nRejectedMinPathDistPaths = 0;
	int nRejectedThresholdPaths = 0;

	// Loop through all times
	for (int t = 0; t < sTimes-1; t++) {

		// Loop through all remaining segments
		while (vecPathSegmentsSet[t].size() > 0) {

			// Create a new path
			Path path;

			PathSegmentSet::iterator iterSeg
				= vecPathSegmentsSet[t].begin();

			path.m_iTimes.push_back(iterSeg->m_iTime[0]);
			path.m_iCandidates.push_back(iterSeg->m_iCandidate[0]);

			int tx = t;

			for (;;) {
				path.m_iTimes.push_back(iterSeg->m_iTime[1]);
				path.m_iCandidates.push_back(iterSeg->m_iCandidate[1]);

				int txnext = iterSeg->m_iTime[1];

				if (txnext >= sTimes-1) {
					vecPathSegmentsSet[tx].erase(iterSeg);
					break;
				}

				PathSegment segFind(
					iterSeg->m_iTime[1], iterSeg->m_iCandidate[1], 0, 0);

				vecPathSegmentsSet[tx].erase(iterSeg);

				iterSeg = vecPathSegmentsSet[txnext].find(segFind);

				if (iterSeg == vecPathSegmentsSet[txnext].end()) {
					break;
				}

				tx = txnext;
			}

			// Reject path due to minimum length
			if (path.m_iTimes.size() < param.nMinDuration) {
				nRejectedMinLengthPaths++;
				continue;
			}

			// Reject path due to minimum endpoint distance
			if (param.dMinEndpointDistance > 0.0) {
				int nT = path.m_iTimes.size();

				int iTime0 = path.m_iTimes[0];
				int iRes0  = path.m_iCandidates[0];

				int iTime1 = path.m_iTimes[nT-1];
				int iRes1  = path.m_iCandidates[nT-1];

				double dDX = vecNodes[iTime0][iRes0].x - vecNodes[iTime1][iRes1].x;
				double dDY = vecNodes[iTime0][iRes0].y - vecNodes[iTime1][iRes1].y;
				double dDZ = vecNodes[iTime0][iRes0].z - vecNodes[iTime1][iRes1].z;

				double dCoord = 0.5 * sqrt(dDX * dDX + dDY * dDY + dDZ * dDZ);
				double dR;
				if (dCoord >= 1.0) {
					dR = 180.0;
				} else {
					dR = 2.0 * asin(dCoord) * 180.0 / M_PI;
				}
				if (dR != dR) {
					_EXCEPTIONT("NaN value detected");
				}

				if (dR < param.dMinEndpointDistance) {
					nRejectedMinEndpointDistPaths++;
					continue;
				}
			}

			// Reject path due to minimum total path distance
			if (param.dMinPathDistance > 0.0) {
				double dTotalPathDistance = 0.0;
				for (int i = 0; i < path.m_iTimes.size() - 1; i++) {
					int iTime0 = path.m_iTimes[i];
					int iRes0 = path.m_iCandidates[i];

					int iTime1 = path.m_iTimes[i+1];
					int iRes1 = path.m_iCandidates[i+1];

					double dDX = vecNodes[iTime0][iRes0].x - vecNodes[iTime1][iRes1].x;
					double dDY = vecNodes[iTime0][iRes0].y - vecNodes[iTime1][iRes1].y;
					double dDZ = vecNodes[iTime0][iRes0].z - vecNodes[iTime1][iRes1].z;

					double dCoord = 0.5 * sqrt(dDX * dDX + dDY * dDY + dDZ * dDZ);
					double dR;
					if (dCoord >= 1.0) {
						dR = 180.0;
					} else {
						dR = 2.0 * asin(dCoord) * 180.0 / M_PI;
					}
					if (dR != dR) {
						_EXCEPTIONT("NaN value detected");
					}

					dTotalPathDistance += dR;
				}

				if (dTotalPathDistance < param.dMinPathDistance) {
					nRejectedMinPathDistPaths++;
					continue;
				}
			}

			// Reject path due to threshold
			bool fOpResult = true;
			for (int x = 0; x < vecPathThresholdOp.size(); x++) {
				fOpResult =
					vecPathThresholdOp[x].Apply(
						path,
						pobjPointData,
						vecFirstCandidate);

				if (!fOpResult) {
					break;
				}
			}
			if (!fOpResult) {
				nRejectedThresholdPaths++;
				continue;
			}

			// Add path to array of paths
			vecPaths.push_back(path);
		}
	}

	Announce("Paths rejected (minlength): %i", nRejectedMinLengthPaths);
	Announce("Paths rejected (minendpointdist): %i", nRejectedMinEndpointDistPaths);
	Announce("Paths rejected (minpathdist): %i", nRejectedMinPathDistPaths);
	Announce("Paths rejected (threshold): %i", nRejectedThresholdPaths);
	Announce("Total paths found: %i", vecPaths.size());
	AnnounceEndBlock("Done");

	// Store results
	if (pobjPathData != NULL) {
		pobjPathData->SetHeaders(*pobjPointData);

		std::vector< MultiTypeDataArray2D > & vecPathData =
			pobjPathData->GetData();

		std::vector<std::string> & vecPathHeaders =
			pobjPathData->GetHeaders();

		vecPathHeaders.resize(vecPaths.size());

		vecPathData.resize(vecPaths.size());

		for (int i = 0; i < vecPaths.size(); i++) {
			vecPathData[i].SetFieldCount(
				pobjPointData->GetIntFieldCount(),
				pobjPointData->GetFloatFieldCount(),
				pobjPointData->GetDoubleFieldCount());

			vecPathData[i].Resize(vecPaths[i].m_iTimes.size());

			if (vecPaths[i].m_iTimes.size() == 0) {
				_EXCEPTIONT("Invalid path");
			}

			const Time & time =
				pobjFileList->GetTime(vecPaths[i].m_iTimes[0]);

			// Construct path headers
			std::string strPathHeader;
			size_t sFormatLen = param.strPathHeaderFormat.length();
			for (size_t j = 0; j < sFormatLen; j++) {
				if (param.strPathHeaderFormat[j] == '%') {

					// Number of times in path
					if (param.strPathHeaderFormat.substr(j,6) == "%count") {
						strPathHeader +=
							std::to_string(
								static_cast<unsigned long>(vecPaths[i].m_iTimes.size()));
						j += 5;

					// Path index
					} else if (param.strPathHeaderFormat.substr(j,6) == "%path") {
						strPathHeader +=
							std::to_string(static_cast<unsigned long>(i));
						j += 4;

					// Year from first path
					} else if (param.strPathHeaderFormat.substr(j,5) == "%year") {
						strPathHeader +=
							std::to_string(static_cast<long long>(time.GetYear()));
						j += 4;

					// Month from first path
					} else if (param.strPathHeaderFormat.substr(j,6) == "%month") {
						strPathHeader +=
							std::to_string(static_cast<long long>(time.GetMonth()));
						j += 5;

					// Day from first path
					} else if (param.strPathHeaderFormat.substr(j,4) == "%day") {
						strPathHeader +=
							std::to_string(static_cast<long long>(time.GetDay()));
						j += 3;

					// Hour from first path
					} else if (param.strPathHeaderFormat.substr(j,5) == "%hour") {
						strPathHeader +=
							std::to_string(
								static_cast<long long>(time.GetSecond() / 3600));
						j += 4;

					// Second from first path
					} else if (param.strPathHeaderFormat.substr(j,7) == "%second") {
						strPathHeader +=
							std::to_string(
								static_cast<long long>(time.GetSecond() % 3600));
						j += 6;

					} else {
						return std::string("ERROR: Invalid \% sequence in header");
					}

				} else {
					strPathHeader += param.strPathHeaderFormat[j];
				}
			}
			vecPathHeaders[i] = strPathHeader;
/*
			for (int j = 0; j < param.pvecPathHeaders->size(); j++) {
				const std::string & strHeader = (*param.pvecPathHeaders)[j];

				const Time & time =
					pobjFileList->GetTime(vecPaths[i].m_iTimes[0]);

				if (strHeader == "%count") {
					vecPathHeaders[i].push_back(vecPaths[i].m_iTimes.size());
				} else if (strHeader == "%year") {
					vecPathHeaders[i].push_back(time.GetYear());
				} else if (strHeader == "%month") {
					vecPathHeaders[i].push_back(time.GetMonth());
				} else if (strHeader == "%day") {
					vecPathHeaders[i].push_back(time.GetDay());
				} else if (strHeader == "%hour") {
					vecPathHeaders[i].push_back(time.GetSecond() / 3600);
				} else if (strHeader == "%second") {
					vecPathHeaders[i].push_back(time.GetSecond() % 3600);
				} else {
					_EXCEPTIONT("Invalid header");
				}
			}
*/

			for (int t = 0; t < vecPaths[i].m_iTimes.size(); t++) {
				int iTime = vecPaths[i].m_iTimes[t];
				int iCandidate = vecPaths[i].m_iCandidates[t];

				int tlocal = vecPaths[i].m_iTimes[t];
				int iGlobalIx = vecFirstCandidate[tlocal] + iCandidate;

				for (int j = 0; j < pobjPointData->GetIntFieldCount(); j++) {
					vecPathData[i].DataInt(t,j) =
						pobjPointData->DataInt(iGlobalIx,j);
				}
				for (int j = 0; j < pobjPointData->GetFloatFieldCount(); j++) {
					vecPathData[i].DataFloat(t,j) =
						pobjPointData->DataFloat(iGlobalIx,j);
				}
				for (int j = 0; j < pobjPointData->GetDoubleFieldCount(); j++) {
					vecPathData[i].DataDouble(t,j) =
						pobjPointData->DataDouble(iGlobalIx,j);
				}
			}
		}
	}

	// Cleanup
	AnnounceStartBlock("Cleanup");

	for (int t = 0; t < vecKDTrees.size(); t++) {
		kd_free(vecKDTrees[t]);
	}

	AnnounceEndBlock("Done");

	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

std::string StitchNodesFunction::Call(
	const ObjectRegistry & objreg,
	const std::vector<std::string> & vecCommandLine,
	const std::vector<ObjectType> & vecCommandLineType,
	Object ** ppReturn
) {
	AnnounceStartBlock("BEGIN stitch_nodes");

	if (vecCommandLine.size() != 3) {
		return std::string("stitch_nodes(config, candidates, param) expects 3 arguments");
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

	// Candidates
	PointDataObject * pobjPointData =
		dynamic_cast<PointDataObject *>(objreg.GetObject(vecCommandLine[1]));
	if (pobjPointData == NULL) {
		return std::string("Second argument to point_search must be of type point_data");
	}

	// Parameters
	Object * pobjParam =
		dynamic_cast<Object *>(objreg.GetObject(vecCommandLine[2]));
	if (pobjParam == NULL) {
		return std::string("Third argument to point_search must be of type parameter_list");
	}

	StitchNodesParam snparam;

	// Range (in degrees)
	StringObject * pobjRange =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("range")));
	if (pobjRange != NULL) {
		bool fSuccess =
			pobjRange->ToUnit(
				"deg", &(snparam.dRange));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("range")
				+ std::string(" to \"deg\"");
		}
	}

	// Minimum duration (in time units)
	IntegerObject * pobjMinDuration =
		dynamic_cast<IntegerObject *>(
			objreg.GetObject(pobjParam->ChildName("minduration")));
	if (pobjMinDuration != NULL) {
		snparam.nMinDuration = pobjMinDuration->Value();
	}

	// Minimum distance between endpoints (in degrees)
	StringObject * pobjMinEndpointDist =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("minendpointdist")));
	if (pobjMinEndpointDist != NULL) {
		bool fSuccess =
			pobjMinEndpointDist->ToUnit(
				"deg", &(snparam.dMinEndpointDistance));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("minendpointdist")
				+ std::string(" to \"deg\"");
		}
	}

	// Minimum distance between endpoints (in degrees)
	StringObject * pobjMinPathDist =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("minpathdist")));
	if (pobjMinPathDist != NULL) {
		bool fSuccess =
			pobjMinPathDist->ToUnit(
				"deg", &(snparam.dMinPathDistance));
		if (!fSuccess) {
			return std::string("Cannot convert ")
				+ pobjParam->ChildName("minpathdist")
				+ std::string(" to \"deg\"");
		}
	}

	// Maximum gap size (in time units)
	IntegerObject * pobjMaxGapSize =
		dynamic_cast<IntegerObject *>(
			objreg.GetObject(pobjParam->ChildName("maxgap")));
	if (pobjMaxGapSize != NULL) {
		snparam.nMaxGapSize = pobjMaxGapSize->Value();
	}

	// Time stride (in time units)
	IntegerObject * pobjTimeStride =
		dynamic_cast<IntegerObject *>(
			objreg.GetObject(pobjParam->ChildName("timestride")));
	if (pobjTimeStride != NULL) {
		snparam.nTimeStride = pobjTimeStride->Value();
	}

	// Threshold operators
	std::vector<PathThresholdOp> vecPathThresholdOps;

	ListObject * pobjPathThresholdCmdList =
		dynamic_cast<ListObject *>(
			objreg.GetObject(pobjParam->ChildName("thresholdcmd")));

	if (pobjPathThresholdCmdList != NULL) {
		for (size_t i = 0; i < pobjPathThresholdCmdList->ChildrenCount(); i++) {
			PathThresholdOp opPathThreshold;

			Object * pobjChild = pobjPathThresholdCmdList->GetChild(i);

			std::string strError =
				opPathThreshold.InitializeFromObject(
					objreg, pobjConfig, pobjPointData, pobjChild);

			if (strError != "") {
				return strError;
			}

			vecPathThresholdOps.push_back(opPathThreshold);
		}
	}

	snparam.pvecPathThresholdOp = &vecPathThresholdOps;

	// Path header format
	StringObject * pobjPathHeaderFormat =
		dynamic_cast<StringObject *>(
			objreg.GetObject(pobjParam->ChildName("pathhead")));

	if (pobjPathHeaderFormat != NULL) {
		snparam.strPathHeaderFormat = pobjPathHeaderFormat->Value();
	} else {
		snparam.strPathHeaderFormat = "%count";
	}

	// Gather data to root node
	pobjPointData->Gather();

	// Execute a serial process
	int nRank = 0;

#if defined(HYPERION_MPIOMP)
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
#endif

	if (nRank == 0) {
		PathDataObject * pobjPathData = NULL;
		if (ppReturn != NULL) {
			pobjPathData = new PathDataObject("");
		}

		std::string strError =
			StitchNodes(
				pobjConfig,
				snparam,
				pobjPointData,
				pobjPathData
			);

		if (ppReturn != NULL) {
			(*ppReturn) = pobjPathData;
		}
		if (strError != "") {
			return strError;
		}

	} else {
		if (ppReturn != NULL) {
			(*ppReturn) = new PathDataObject("");
		}
	}

#if defined(HYPERION_MPIOMP)
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	AnnounceEndBlock("Done");
	return std::string("");
}

///////////////////////////////////////////////////////////////////////////////

};


