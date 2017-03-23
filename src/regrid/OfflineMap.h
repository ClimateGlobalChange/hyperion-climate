///////////////////////////////////////////////////////////////////////////////
///
///	\file    OfflineMap.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#ifndef _OFFLINEMAP_H_
#define _OFFLINEMAP_H_

#include "SparseMatrix.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include <string>
#include <vector>

class Mesh;

namespace HRegrid {

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An offline map between two Meshes.
///	</summary>
class OfflineMap {

public:
	///	<summary>
	///		Initialize the array of input dimensions from arrays.
	///	</summary>
	void InitializeSourceDimensions(
		const std::vector<int> & vecSourceDimSizes,
		const std::vector<std::string> & vecSourceDimNames
	) {
		m_vecSourceDimSizes = vecSourceDimSizes;
		m_vecSourceDimNames = vecSourceDimNames;
	}

	///	<summary>
	///		Initialize the array of output dimensions from arrays.
	///	</summary>
	void InitializeTargetDimensions(
		const std::vector<int> & vecTargetDimSizes,
		const std::vector<std::string> & vecTargetDimNames
	) {
		m_vecTargetDimSizes = vecTargetDimSizes;
		m_vecTargetDimNames = vecTargetDimNames;
	}

	///	<summary>
	///		Initialize the array of input dimensions from a file.
	///	</summary>
	void InitializeSourceDimensionsFromFile(
		const std::string & strSourceMesh
	);

	///	<summary>
	///		Initialize the array of output dimensions from a file.
	///	</summary>
	void InitializeTargetDimensionsFromFile(
		const std::string & strTargetMesh
	);

private:
	///	<summary>
	///		Initialize the coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeCoordinatesFromMeshFV(
		const Mesh & mesh,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray2D<double> & dVertexLon,
		DataArray2D<double> & dVertexLat,
		bool fLatLon
	);

	///	<summary>
	///		Initialize the coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeCoordinatesFromMeshFE(
		const Mesh & mesh,
		int nP,
		const DataArray3D<int> & dataGLLnodes,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray2D<double> & dVertexLon,
		DataArray2D<double> & dVertexLat
	);

	///	<summary>
	///		Initialize the rectilinear coordinate vectors.
	///	</summary>
	void InitializeRectilinearCoordinateVector(
		int nLon,
		int nLat,
		const DataArray2D<double> & dVertexLon,
		const DataArray2D<double> & dVertexLat,
		bool fLonFirst,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray1D<double> & dVectorCenterLon,
		DataArray1D<double> & dVectorCenterLat,
		DataArray2D<double> & dVectorBoundsLon,
		DataArray2D<double> & dVectorBoundsLat
	);

public:
	///	<summary>
	///		Initialize the source coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeSourceCoordinatesFromMeshFV(
		const Mesh & meshSource
	);

	///	<summary>
	///		Initialize the target coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeTargetCoordinatesFromMeshFV(
		const Mesh & meshTarget
	);

	///	<summary>
	///		Initialize the source coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeSourceCoordinatesFromMeshFE(
		const Mesh & meshSource,
		int nP,
		const DataArray3D<int> & dataGLLnodesSource
	);

	///	<summary>
	///		Initialize the target coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeTargetCoordinatesFromMeshFE(
		const Mesh & meshTarget,
		int nP,
		const DataArray3D<int> & dataGLLnodesTarget
	);

public:
	///	<summary>
	///		Copy a list of variables from a source file to target file.
	///	</summary>
	void PreserveVariables(
		const std::string & strSourceDataFile,
		const std::string & strTargetDataFile,
		const std::vector<std::string> & vecPreserveVariables
	);

	///	<summary>
	///		Copy all non-remapped, non-dimensional variables from source file
	///		to target file.
	///	</summary>
	void PreserveAllVariables(
		const std::string & strSourceDataFile,
		const std::string & strTargetDataFile
	);

public:
	///	<summary>
	///		Apply the offline map to a data file.
	///	</summary>
	void Apply(
		const std::string & strSourceDataFile,
		const std::string & strTargetDataFile,
		const std::vector<std::string> & vecVariables,
		const std::string & strNColName,
		bool fTargetDouble = false,
		bool fAppend = false
	);

	///	<summary>
	///		Read the OfflineMap from a NetCDF file.
	///	</summary>
	void Read(
		const std::string & strSource
	);

	///	<summary>
	///		Write the OfflineMap to a NetCDF file.
	///	</summary>
	void Write(
		const std::string & strTarget
	);

public:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	bool IsConsistent(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	bool IsConservative(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	bool IsMonotone(
		double dTolerance
	);

public:
	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	DataArray1D<double> & GetSourceAreas() {
		return m_dSourceAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	const DataArray1D<double> & GetSourceAreas() const {
		return m_dSourceAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	DataArray1D<double> & GetTargetAreas() {
		return m_dTargetAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	const DataArray1D<double> & GetTargetAreas() const {
		return m_dTargetAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	void SetSourceAreas(const DataArray1D<double> & dSourceAreas) {
		m_dSourceAreas = dSourceAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	void SetTargetAreas(const DataArray1D<double> & dTargetAreas) {
		m_dTargetAreas = dTargetAreas;
	}

	///	<summary>
	///		Get the SparseMatrix representation of the OfflineMap.
	///	</summary>
	const SparseMatrix<double> & GetSparseMatrix() const {
		return m_mapRemap;
	}

	///	<summary>
	///		Get the SparseMatrix representation of the OfflineMap.
	///	</summary>
	SparseMatrix<double> & GetSparseMatrix() {
		return m_mapRemap;
	}

public:
	///	<summary>
	///		Set the fill value override.
	///	</summary>
	void SetFillValueOverride(float flFillValueOverride) {
		m_flFillValueOverride = flFillValueOverride;
	}

protected:
	///	<summary>
	///		The SparseMatrix representing this operator.
	///	</summary>
	SparseMatrix<double> m_mapRemap;

	///	<summary>
	///		Vector of areas associated with input degrees of freedom.
	///	</summary>
	DataArray1D<double> m_dSourceAreas;

	///	<summary>
	///		Vector of areas associated with output degrees of freedom.
	///	</summary>
	DataArray1D<double> m_dTargetAreas;

protected:
	///	<summary>
	///		Vector of cell center longitudes on source grid.
	///	</summary>
	DataArray1D<double> m_dSourceCenterLon;

	///	<summary>
	///		Vector of cell center latitudes on source grid.
	///	</summary>
	DataArray1D<double> m_dSourceCenterLat;

	///	<summary>
	///		Vector of cell center longitudes on target grid.
	///	</summary>
	DataArray1D<double> m_dTargetCenterLon;

	///	<summary>
	///		Vector of cell center latitudes on target grid.
	///	</summary>
	DataArray1D<double> m_dTargetCenterLat;

	///	<summary>
	///		Matrix of vertex longitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dSourceVertexLon;

	///	<summary>
	///		Matrix of vertex latitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dSourceVertexLat;

	///	<summary>
	///		Vector containing cell center longitude along "lon" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorSourceCenterLon;

	///	<summary>
	///		Vector containing cell center latitude along "lat" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorSourceCenterLat;

	///	<summary>
	///		Vector containing bounds for longitude along "lon" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorSourceBoundsLon;

	///	<summary>
	///		Vector containing bounds for latitude along "lat" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorSourceBoundsLat;

	///	<summary>
	///		Matrix of vertex longitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dTargetVertexLon;

	///	<summary>
	///		Matrix of vertex latitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dTargetVertexLat;

	///	<summary>
	///		Vector containing cell center longitude along "lon" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorTargetCenterLon;

	///	<summary>
	///		Vector containing cell center latitude along "lat" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorTargetCenterLat;

	///	<summary>
	///		Vector containing bounds for longitude along "lon" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorTargetBoundsLon;

	///	<summary>
	///		Vector containing bounds for latitude along "lat" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorTargetBoundsLat;

protected:
	///	<summary>
	///		Vector of dimension sizes for source.
	///	</summary>
	std::vector<int> m_vecSourceDimSizes;

	///	<summary>
	///		Vector of dimension names for source.
	///	</summary>
	std::vector<std::string> m_vecSourceDimNames;

	///	<summary>
	///		Vector of dimension sizes for target.
	///	</summary>
	std::vector<int> m_vecTargetDimSizes;

	///	<summary>
	///		Vector of dimension names for target.
	///	</summary>
	std::vector<std::string> m_vecTargetDimNames;

	///	<summary>
	///		The fill value override.
	///	</summary>
	float m_flFillValueOverride;
};

///////////////////////////////////////////////////////////////////////////////

};

#endif

