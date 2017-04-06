///////////////////////////////////////////////////////////////////////////////
///
///	\file    PointDataObject.h
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

#ifndef _POINTDATAOBJECT_H_
#define _POINTDATAOBJECT_H_

#include "Announce.h"
#include "Object.h"
#include "MultiTypeDataArray2D.h"
#include "MultiTypeDataHeaders.h"

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class PointDataObject :
	public DistributedObject,
	public MultiTypeDataArray2D,
	public MultiTypeDataHeaders
{

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PointDataObject(const std::string & strName) :
		m_nTimeStride(1),
		DistributedObject(strName),
		MultiTypeDataArray2D(),
		MultiTypeDataHeaders()
	{ }

	///	<summary>
	///		Self-duplicator.
	///	</summary>
	virtual Object * Duplicate(
		const std::string & strDuplicateName,
		ObjectRegistry & objreg
	) const {
		return _Duplicate(
			new PointDataObject(strDuplicateName),
			objreg);
	}

	///	<summary>
	///		Call a member function of this Object.
	///	</summary>
	virtual std::string Call(
		const ObjectRegistry & objreg,
		const std::string & strFunctionName,
		const std::vector<std::string> & vecCommandLine,
		const std::vector<ObjectType> & vecCommandLineType,
		Object ** ppReturn
	) {
		if (strFunctionName == "output_csv") {
			if ((vecCommandLineType.size() != 1) ||
			    (vecCommandLineType[0] != ObjectType_String)
			) {
				return std::string("ERROR: Invalid parameters to function \"output_csv\"");
			}
			return OutputCSV(vecCommandLine[0]);
		}
		return
			Object::Call(
				objreg,
				strFunctionName,
				vecCommandLine,
				vecCommandLineType,
				ppReturn);
	}


	///	<summary>
	///		Output this object as a CSV file.
	///	</summary>
	std::string OutputCSV(
		const std::string & strCSVOutputFilename
	);

public:
	///	<summary>
	///		Set the number of fields used by the object.
	///	</summary>
	virtual void SetFieldCount(
		size_t nIntFields = 0,
		size_t nFloatFields = 0,
		size_t nDoubleFields = 0
	) {
		MultiTypeDataArray2D::SetFieldCount(
			nIntFields, nFloatFields, nDoubleFields);

		MultiTypeDataHeaders::SetFieldCount(
			nIntFields, nFloatFields, nDoubleFields);
	}

	///	<summary>
	///		Concatenate multiple PointDataObjects.
	///	</summary>
	virtual void Concatenate(
		const std::vector<PointDataObject *> & vecpobjPointData
	) {
		if (vecpobjPointData.size() == 0) {
			return;
		}

		m_nTimeStride = vecpobjPointData[0]->m_nTimeStride;

		std::vector<MultiTypeDataArray2D *> vecpMultiTypeData;
		vecpMultiTypeData.resize(vecpobjPointData.size());
		for (size_t i = 0; i < vecpobjPointData.size(); i++) {
			vecpMultiTypeData[i] = vecpobjPointData[i];
			if (m_nTimeStride != vecpobjPointData[i]->GetTimeStride()) {
				_EXCEPTIONT("TimeStride mismatch");
			}
		}
		MultiTypeDataArray2D::Concatenate(vecpMultiTypeData);
	}

	///	<summary>
	///		Gather data to the head rank.
	///	</summary>
	virtual void Gather() {
		if (IsDistributed()) {
			MultiTypeDataArray2D::Gather();
		}

		DistributedObject::UnsetDistributed();
	}

public:
	///	<summary>
	///		Get the time stride.
	///	</summary>
	int GetTimeStride() const {
		return m_nTimeStride;
	}

	///	<summary>
	///		Set the time stride.
	///	</summary>
	void SetTimeStride(int nTimeStride) {
		m_nTimeStride = nTimeStride;
	}

protected:
	///	<summary>
	///		Time stride used for this point data.
	///	</summary>
	int m_nTimeStride;
};

///////////////////////////////////////////////////////////////////////////////

#endif

