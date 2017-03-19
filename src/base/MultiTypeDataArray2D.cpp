///////////////////////////////////////////////////////////////////////////////
///
///	\file    MultiTypeDataArray2D.cpp
///	\author  Paul Ullrich
///	\version March 17, 2017
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

#include "MultiTypeDataArray2D.h"

#include "DataArray1D.h"

#if defined(HYPERION_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

void MultiTypeDataArray2D::SetFieldCount(
	size_t nIntFields,
	size_t nFloatFields,
	size_t nDoubleFields
) {
	m_nIntFields = nIntFields;
	m_nFloatFields = nFloatFields;
	m_nDoubleFields = nDoubleFields;

	m_dDataInt.Deallocate();
	m_dDataFloat.Deallocate();
	m_dDataDouble.Deallocate();

	if (m_nRows != 0) {
		if (m_nIntFields != 0) {
			m_dDataInt.SetSize(m_nRows, m_nIntFields);
		}
		if (m_nFloatFields != 0) {
			m_dDataFloat.SetSize(m_nRows, m_nFloatFields);
		}
		if (m_nDoubleFields != 0) {
			m_dDataDouble.SetSize(m_nRows, m_nDoubleFields);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void MultiTypeDataArray2D::Resize(
	size_t nRows
) {
	if (nRows == m_nRows) {
		return;
	}

	if (nRows == 0) {
		if (m_nIntFields != 0) {
			m_dDataInt.Deallocate();
		}
		if (m_nFloatFields != 0) {
			m_dDataFloat.Deallocate();
		}
		if (m_nDoubleFields != 0) {
			m_dDataDouble.Deallocate();
		}
		return;
	}

	if (m_nRows == 0) {
		m_nRows = nRows;
		if (m_nIntFields != 0) {
			m_dDataInt.Allocate(m_nRows, m_nIntFields);
		}
		if (m_nFloatFields != 0) {
			m_dDataFloat.Allocate(m_nRows, m_nFloatFields);
		}
		if (m_nDoubleFields != 0) {
			m_dDataDouble.Allocate(m_nRows, m_nDoubleFields);
		}
		return;
	}

	int nCopy = m_nRows;
	if (nCopy > nRows) {
		nCopy = nRows;
	}

	if (m_nIntFields != 0) {
		DataArray2D<int> dTempInt = m_dDataInt;
		m_dDataInt.Deallocate();
		m_dDataInt.Allocate(nRows, m_nIntFields);
		memcpy(&(m_dDataInt(0,0)), &(dTempInt(0,0)),
			sizeof(int) * nCopy * m_nIntFields);
	}

	if (m_nFloatFields != 0) {
		DataArray2D<float> dTempFloat = m_dDataFloat;
		m_dDataFloat.Deallocate();
		m_dDataFloat.Allocate(nRows, m_nFloatFields);
		memcpy(&(m_dDataFloat(0,0)), &(dTempFloat(0,0)),
			sizeof(float) * nCopy * m_nFloatFields);
	}
	if (m_nDoubleFields != 0) {
		DataArray2D<double> dTempDouble = m_dDataDouble;
		m_dDataDouble.Deallocate();
		m_dDataDouble.Allocate(nRows, m_nDoubleFields);
		memcpy(&(m_dDataDouble(0,0)), &(dTempDouble(0,0)),
			sizeof(double) * nCopy * m_nDoubleFields);
	}

	m_nRows = nRows;
}

///////////////////////////////////////////////////////////////////////////////

void MultiTypeDataArray2D::Concatenate(
	const std::vector<MultiTypeDataArray2D *> & vecpobjPointData
) {
	if (vecpobjPointData.size() == 0) {
		return;
	}

	size_t sRows = 0;
	for (size_t s = 0; s < vecpobjPointData.size(); s++) {
		sRows += vecpobjPointData[s]->GetRows();
	}

	SetFieldCount(
		vecpobjPointData[0]->m_nIntFields,
		vecpobjPointData[0]->m_nFloatFields,
		vecpobjPointData[0]->m_nDoubleFields);

	Resize(sRows);

	sRows = 0;
	for (size_t s = 0; s < vecpobjPointData.size(); s++) {
		if (vecpobjPointData[s]->m_nIntFields != m_nIntFields) {
			_EXCEPTIONT("int field count mismatch");
		}
		if (vecpobjPointData[s]->m_nFloatFields != m_nFloatFields) {
			_EXCEPTIONT("float field count mismatch");
		}
		if (vecpobjPointData[s]->m_nDoubleFields != m_nDoubleFields) {
			_EXCEPTIONT("double field count mismatch");
		}

		memcpy(
			&(m_dDataInt(sRows,0)),
			&(vecpobjPointData[s]->m_dDataInt(0,0)),
			m_nIntFields * vecpobjPointData[s]->m_nRows * sizeof(int));

		memcpy(
			&(m_dDataFloat(sRows,0)),
			&(vecpobjPointData[s]->m_dDataFloat(0,0)),
			m_nFloatFields * vecpobjPointData[s]->m_nRows * sizeof(float));

		memcpy(
			&(m_dDataDouble(sRows,0)),
			&(vecpobjPointData[s]->m_dDataDouble(0,0)),
			m_nDoubleFields * vecpobjPointData[s]->m_nRows * sizeof(double));

		sRows += vecpobjPointData[s]->GetRows();
	}
}

///////////////////////////////////////////////////////////////////////////////

void MultiTypeDataArray2D::Gather() {

#if defined(HYPERION_MPIOMP)
	// Record the total number of rows across all threads
	int nLocalRows = GetRows();

	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	int nSize;
	MPI_Comm_size(MPI_COMM_WORLD, &nSize);

	// Gather total row count on each process
	DataArray1D<int> vecAllRows;
	DataArray1D<int> vecAllRowsInt;
	DataArray1D<int> vecAllRowsIntDispl;
	DataArray1D<int> vecAllRowsFloat;
	DataArray1D<int> vecAllRowsFloatDispl;
	DataArray1D<int> vecAllRowsDouble;
	DataArray1D<int> vecAllRowsDoubleDispl;

	if (nRank == 0) {
		vecAllRows.Allocate(nSize);
		vecAllRowsInt.Allocate(nSize);
		vecAllRowsIntDispl.Allocate(nSize);
		vecAllRowsFloat.Allocate(nSize);
		vecAllRowsFloatDispl.Allocate(nSize);
		vecAllRowsDouble.Allocate(nSize);
		vecAllRowsDoubleDispl.Allocate(nSize);
	}

	int iError =
		MPI_Gather(
			&nLocalRows,
			1,
			MPI_INT,
			&(vecAllRows[0]),
			1,
			MPI_INT,
			0,
			MPI_COMM_WORLD);

	if (iError != MPI_SUCCESS) {
		_EXCEPTION1("MPI_Gather error: %i", iError);
	}

	// Determine recv buffer sizes and displacements
	int nTotalRows = 0;
	if (nRank == 0) {
		for (int i = 0; i < nSize; i++) {
			nTotalRows += vecAllRows[i];
			vecAllRowsInt[i] = vecAllRows[i] * m_nIntFields;
			vecAllRowsFloat[i] = vecAllRows[i] * m_nFloatFields;
			vecAllRowsDouble[i] = vecAllRows[i] * m_nDoubleFields;
		}

		vecAllRowsIntDispl[0] = 0;
		vecAllRowsFloatDispl[0] = 0;
		vecAllRowsDoubleDispl[0] = 0;
		for (int i = 1; i < nSize; i++) {
			vecAllRowsIntDispl[i] =
				vecAllRowsIntDispl[i-1] + vecAllRowsInt[i-1];
			vecAllRowsFloatDispl[i] =
				vecAllRowsFloatDispl[i-1] + vecAllRowsFloat[i-1];
			vecAllRowsDoubleDispl[i] =
				vecAllRowsDoubleDispl[i-1] + vecAllRowsDouble[i-1];
		}

		Resize(nTotalRows);
	}

	// Send total number of rows to all processes
	iError =
		MPI_Scatter(
			&nTotalRows,
			1,
			MPI_INT,
			&nTotalRows,
			1,
			MPI_INT,
			0,
			MPI_COMM_WORLD);

	if (iError != MPI_SUCCESS) {
		_EXCEPTION1("MPI_Scatter error: %i", iError);
	}

	if (nTotalRows == 0) {
		return;
	}

	// Gather int fields
	if (m_nIntFields != 0) {
		iError =
			MPI_Gatherv(
				&(m_dDataInt(0,0)),
				nLocalRows * m_nIntFields,
				MPI_INT,
				&(m_dDataInt(0,0)),
				&(vecAllRowsInt[0]),
				&(vecAllRowsIntDispl[0]),
				MPI_INT,
				0,
				MPI_COMM_WORLD);

		if (iError != MPI_SUCCESS) {
			_EXCEPTION1("MPI_Gatherv error: %i", iError);
		}
	}

	if (m_nFloatFields != 0) {
		iError =
			MPI_Gatherv(
				&(m_dDataFloat(0,0)),
				nLocalRows * m_nFloatFields,
				MPI_FLOAT,
				&(m_dDataFloat(0,0)),
				&(vecAllRowsFloat[0]),
				&(vecAllRowsFloatDispl[0]),
				MPI_FLOAT,
				0,
				MPI_COMM_WORLD);

		if (iError != MPI_SUCCESS) {
			_EXCEPTION1("MPI_Gatherv error: %i", iError);
		}
	}

	if (m_nDoubleFields != 0) {
		iError =
			MPI_Gatherv(
				&(m_dDataDouble(0,0)),
				nLocalRows * m_nDoubleFields,
				MPI_DOUBLE,
				&(m_dDataDouble(0,0)),
				&(vecAllRowsDouble[0]),
				&(vecAllRowsDoubleDispl[0]),
				MPI_DOUBLE,
				0,
				MPI_COMM_WORLD);

		if (iError != MPI_SUCCESS) {
			_EXCEPTION1("MPI_Gatherv error: %i", iError);
		}
	}

#endif
}
	
///////////////////////////////////////////////////////////////////////////////

void MultiTypeDataArray2D::OutputCSV(
	std::ostream & osOutput
) const {
	bool fInitialComma;

	// Loop through all rows
	for (int i = 0; i < GetRows(); i++) {
		fInitialComma = false;

		// Output integer data
		for (int j = 0; j < GetIntFieldCount(); j++) {
			if (!fInitialComma) {
				fInitialComma = true;
			} else {
				osOutput << ",";
			}
			osOutput << DataInt(i,j);
		}
		for (int j = 0; j < GetFloatFieldCount(); j++) {
			if (!fInitialComma) {
				fInitialComma = true;
			} else {
				osOutput << ",";
			}
			osOutput << DataFloat(i,j);
		}
		for (int j = 0; j < GetDoubleFieldCount(); j++) {
			if (!fInitialComma) {
				fInitialComma = true;
			} else {
				osOutput << ",";
			}
			osOutput << DataDouble(i,j);
		}
		osOutput << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

