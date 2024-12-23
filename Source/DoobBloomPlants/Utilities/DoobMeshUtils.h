#pragma once

#include "CoreMinimal.h"

namespace DoobMeshUtils {
	void SmoothMeshVertices(TArray<FVector>& Vertices, const TArray<int32>& Triangles, int32 Iterations);
	void RecalculateNormals(const TArray<FVector>& Vertices, const TArray<int32>& Triangles, TArray<FVector>& Normals);
	void RemoveDegenerateTriangles(TArray<FVector>& Vertices, TArray<int32>& Triangles);
}