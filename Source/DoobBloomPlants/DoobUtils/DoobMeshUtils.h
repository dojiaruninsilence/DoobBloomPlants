#pragma once

#include "CoreMinimal.h"

/**
 * @namespace DoobMeshUtils
 * Utility namespace for mesh operations.
 */
namespace DoobMeshUtils {
	/**
	 * Smooths the vertices of a mesh by averaging neighboring vertices.
	 *
	 * @param Vertices The array of vertices to be smoothed.
	 * @param Triangles The array of triangle indices defining the mesh.
	 * @param Iterations The number of smoothing iterations to perform.
	 */
	void SmoothMeshVertices(TArray<FVector>& Vertices, const TArray<int32>& Triangles, int32 Iterations);

	/**
	 * Recalculates vertex normals for a mesh.
	 *
	 * @param Vertices The array of vertices defining the mesh.
	 * @param Triangles The array of triangle indices defining the mesh.
	 * @param Normals The output array of recalculated normals.
	 */
	void RecalculateNormals(const TArray<FVector>& Vertices, const TArray<int32>& Triangles, TArray<FVector>& Normals);

	/**
	* Removes degenerate triangles from a mesh.
	*
	* @param Vertices The array of vertices defining the mesh.
	* @param Triangles The array of triangle indices defining the mesh. This array will be modified to exclude degenerate triangles.
	*/
	void RemoveDegenerateTriangles(TArray<FVector>& Vertices, TArray<int32>& Triangles);

	/**
	 * Generates the triangle indices for a cylindrical tube mesh.
	 *
	 * @param TubeVertices The array of vertices defining the tube.
	 * @param ProfileSegments The number of segments in each profile.
	 * @param OutTriangles The output array of triangle indices for the tube mesh.
	 */
	void GenerateTubeMesh(
		const TArray<FVector>& TubeVertices,
		int32 ProfileSegments,
		TArray<int32>& OutTriangles
	);

	FVector CalculateMeshCentroid(const TArray<FVector>& Vertices);

	void RemoveInwardFacingTriangles(
		TArray<FVector>& Vertices,
		TArray<int32>& Triangles
	);
}