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
	int nTotalRows = 0;

	int iError =
		MPI_Reduce(
			&nLocalRows,
			&nTotalRows,
			1,
			MPI_INT,
			MPI_SUM,
			0,
			MPI_COMM_WORLD);

	if (iError != MPI_SUCCESS) {
		_EXCEPTION1("MPI_Reduce error: %i", iError);
	}

	Resize(nTotalRows);

	if (nTotalRows == 0) {
		return;
	}

	// Gather int fields
	if (m_nIntFields != 0) {
		int iError =
			MPI_Gather(
				&(m_dDataInt(0,0)),
				m_nRows * m_nIntFields,
				MPI_INT,
				&(m_dDataInt(0,0)),
				nTotalRows * m_nIntFields,
				MPI_INT,
				0,
				MPI_COMM_WORLD);

		if (iError != MPI_SUCCESS) {
			_EXCEPTION1("MPI_Gather error: %i", iError);
		}
	}

	if (m_nFloatFields != 0) {
		int iError =
			MPI_Gather(
				&(m_dDataFloat(0,0)),
				m_nRows * m_nFloatFields,
				MPI_FLOAT,
				&(m_dDataFloat(0,0)),
				nTotalRows * m_nFloatFields,
				MPI_FLOAT,
				0,
				MPI_COMM_WORLD);

		if (iError != MPI_SUCCESS) {
			_EXCEPTION1("MPI_Gather error: %i", iError);
		}
	}

	if (m_nDoubleFields != 0) {
		int iError =
			MPI_Gather(
				&(m_dDataDouble(0,0)),
				m_nRows * m_nDoubleFields,
				MPI_DOUBLE,
				&(m_dDataDouble(0,0)),
				nTotalRows * m_nDoubleFields,
				MPI_DOUBLE,
				0,
				MPI_COMM_WORLD);

		if (iError != MPI_SUCCESS) {
			_EXCEPTION1("MPI_Gather error: %i", iError);
		}
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////

