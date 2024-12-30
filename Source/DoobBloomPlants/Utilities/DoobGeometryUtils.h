#pragma once

#include "CoreMinimal.h"

#include "DoobProfileUtils.h"

namespace DoobGeometryUtils {

    // ------------------------------------------------------ Ring and Tube ------------------------------------------------------------ // 

    struct FRingData {
        TArray<FVector> Vertices;
        FVector Normal;
        FVector Center;
        FVector Direction;
        FVector UpVector;
        float Radius;
        bool bIsClosed = true;
        int32 NumVertices = Vertices.Num();
        FName RingID;
    };

    struct FTubeData {
        DoobProfileUtils::F2DProfile Profile;
        TArray<FRingData> Rings;
        FVector StartPosition;
        FVector EndPosition;
        FVector Direction;
        FVector UpVector;
        float Length;
        int32 NumSegments;
        int32 NumSides;
        int32 NumRings = Rings.Num();
        FName TubeID;
    };

	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, FRingData& RingData);
	void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);
	void ConnectRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);


    // construct the vertices for a cylindrical tube using a 2d profile 
    void ConstructTubeFromProfile(
        const DoobProfileUtils::F2DProfile& Profile,
        const FVector& StartPosition,
        const FVector& Direction,
        const FVector& UpVector,
        int32 NumSegments,
        int32 NumSides,
        float Length,
        FTubeData& OutTube,
        bool ApplyBothAxis
    );

    // ------------------------------------------------------ Planes and Lines Intersections ------------------------------------------------------------ //

    // contain plane equations
    struct FPlaneEquation {
        FVector Normal;
        float D;

        FPlaneEquation();
        FPlaneEquation(const FVector& InNormal, float InD);

        // possibly add methods for point containment or intersection checks here
    };

    // function to calculate plane equation from two rings
    FPlaneEquation CalculatePlaneFromRings(const FRingData& RingA, const FRingData& RingB);

    // calculate a ray intersection with a plane
    bool CalculateIntersectionWithPlane(
        const FPlaneEquation& Plane,
        const FVector& RayOrigin,
        const FVector& RayDirection,
        FVector& OutIntersectionPoint
    );

    bool PointInsideCircle(const FVector& Point, const FVector& Center, float Radius);

    void FilterPlanesAndLines(
        const FTubeData& TubeA,
        const FTubeData& TubeB,
        TArray<int32>& OutPlaneIndicesA,
        TArray<int32>& OutPlaneIndicesB,
        TArray<int32>& OutLineIndicesA,
        TArray<int32>& OutLineIndicesB
    );

    // ------------------------------------------------------ Shape Intersections ------------------------------------------------------------ //

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
}