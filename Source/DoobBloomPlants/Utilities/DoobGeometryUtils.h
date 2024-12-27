#pragma once

#include "CoreMinimal.h"

#include "DoobProfileUtils.h"

//namespace DoobProfileUtils {
//    struct F2DProfile;
//}
namespace DoobGeometryUtils {

	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, TArray<FVector>& RingVertices);
	void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);
	void ConnectRingArray(const TArray<TArray<FVector>>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * Generates the intersection curve of two cylindrical shapes based on their 2D profiles and transformations.
     *
     * @param MainProfile The 2D profile of the main cylinder.
     * @param MainTransform The transformation applied to the main cylinder.
     * @param LateralProfile The 2D profile of the lateral cylinder.
     * @param LateralTransform The transformation applied to the lateral cylinder.
     * @param MainSegments Number of parametric segments to sample on the main cylinder.
     * @param LateralSegments Number of parametric segments to sample on the lateral cylinder.
     * @param Threshold Distance threshold for considering points as intersecting.
     * @return An array of 3D points representing the intersection curve.
     */
    TArray<FVector> GenerateIntersectionCurve(
        const DoobProfileUtils::F2DProfile& MainProfile,
        const FTransform& MainTransform,
        const DoobProfileUtils::F2DProfile& LateralProfile,
        const FTransform& LateralTransform,
        int32 MainSegments,
        int32 LateralSegments,
        float Threshold
    );

    // construct the vertices for a cylindrical tube using a 2d profile
    void ConstructTubeFromProfile(
        const DoobProfileUtils::F2DProfile& Profile,
        const FVector& StartPosition,
        const FVector& EndPosition,
        const FVector& Direction,
        const FVector& UpVector,
        int32 NumSegments,
        int32 NumSides,
        float Length,
        TArray<TArray<FVector>>& OutRingVerticesArrays,
        bool ApplyBothAxis = true
    );
}