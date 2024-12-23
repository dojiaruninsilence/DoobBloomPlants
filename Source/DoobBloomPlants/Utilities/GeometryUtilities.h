#pragma once

#include "CoreMinimal.h"

namespace GeometryUtilities {
	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, TArray<FVector>& RingVertices);
	void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);
	void ConnectRingArray(const TArray<TArray<FVector>>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);
}