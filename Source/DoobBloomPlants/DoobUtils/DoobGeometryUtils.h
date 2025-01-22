#pragma once

#include "CoreMinimal.h"

#include "DoobProfileUtils.h"

/**
 * @namespace DoobGeometryUtils
 * Contains utility functions and data structures for geometric operations.
 */
namespace DoobGeometryUtils {

    // ------------------------------------------------------ Ring and Tube ------------------------------------------------------------ // 

    /**
     * @struct FRingData
     * Stores data related to a ring in 3D space, including vertices, normal, and other properties.
     */
    struct FRingData {
        TArray<FVector> Vertices; ///< The vertices forming the ring.
        FVector Normal; ///< Normal vector of the ring.
        FVector Center; ///< Center position of the ring.
        FVector Direction; ///< Direction the ring is facing.
        FVector UpVector; ///< Up vector for the ring.
        float Radius; ///< Radius of the ring.
        bool bIsClosed = true; ///< Whether the ring is closed.
        int32 NumVertices = Vertices.Num(); ///< Number of vertices in the ring.
        FName RingID; ///< Unique identifier for the ring.
        bool bIsComplete = true; ///< Whether the ring is completed or missing vertices
    };

    struct FIntersectionRingData {
        TArray<FVector> CombinedVertices;
        TArray<FVector> MainTubeVertices;
        TArray<FVector> LateralTubeVertices;
        FVector Centroid;
        TArray<FVector> CardinalVertices; // nesw - 0123
        TArray<int32> CardinalIndices;
        /*TArray<FVector> BottomLeftVertices;
        TArray<FVector> BottomRightVertices;
        TArray<FVector> TopRightVertices;
        TArray<FVector> TopLeftVertices;*/
    };

    struct FIntersectionSquareData {
        TArray<FVector> Corners;
        // TArray<int32> CornerRingConnectionIndices;
        TArray<FVector> LeftSideRingConnections;
        TArray<FVector> RightSideRingConnections;
        FVector Center;
        TArray<float> Angles;
        /*TArray<FVector> BottomLeftVertices;
        TArray<FVector> BottomRightVertices;
        TArray<FVector> TopRightVertices;
        TArray<FVector> TopLeftVertices;*/
        TArray<TArray<FVector>> BottomLeftPartialRings;
        TArray<TArray<FVector>> BottomRightPartialRings;
        TArray<TArray<FVector>> TopRightPartialRings;
        TArray<TArray<FVector>> TopLeftPartialRings;
    };

    /**
     * @struct FTubeData
     * Stores data related to a tube constructed from rings and a 2D profile.
     */
    struct FTubeData {
        DoobProfileUtils::F2DProfile Profile;  ///< 2D profile defining the tube's cross-section.
        TArray<FRingData> Rings; ///< Array of rings forming the tube.
        FVector StartPosition; ///< Starting position of the tube.
        FVector EndPosition; ///< Ending position of the tube.
        FVector Direction; ///< Direction the tube is facing.
        FVector UpVector; ///< Up vector for the tube.
        float Length; ///< Length of the tube.
        int32 NumSegments; ///< Number of segments in the tube.
        int32 NumSides; ///< Number of sides per segment.
        int32 NumRings = Rings.Num(); ///< Number of rings forming the tube.
        FName TubeID; ///< Unique identifier for the tube.
    };

    struct FTwoTubeIntersectionData {
        FTubeData MainTube;
        FTubeData LateralTube;
        FIntersectionRingData IntersectionRing;
        FIntersectionSquareData IntersectionSquare;
        FRingData MTAboveIntersectionRing;
        FRingData MTBelowIntersectionRing;
        FRingData MTAboveIntersectionRingPartial;
        FRingData MTBelowIntersectionRingPartial;
        FTubeData LateralTubeIntersectionRings;
        FRingData LateralTubeFirstFullRing;
        FTubeData LateralTubeRemovedVertices;
        TArray<FVector> AllVertices;
        TArray<int32> Triangles;
        TArray<FRingData> MainTubePartialRings;
    };

    /**
     * Generates a ring of vertices around a center point.
     * @param Center Center of the ring.
     * @param Direction Direction the ring is facing.
     * @param UpVector Up vector for orientation.
     * @param Radius Radius of the ring.
     * @param NumSides Number of sides (vertices) in the ring.
     * @param RingData Output data structure to store the ring information.
     */
    void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, FRingData& RingData);

    void GenerateIntersectionRing(
        const FTubeData& MainTube,
        const FTubeData& LateralTube,
        FIntersectionRingData& OutRing,
        float Precision = KINDA_SMALL_NUMBER
    );

    void GenerateHalfIntersectionRing(
        const FTubeData& MainTube,
        const FTubeData& LateralTube,
        TArray<FVector>& OutRingData,
        float Precision = KINDA_SMALL_NUMBER
    );

    void GenerateSquareAroundIntersection(
        const FRingData& LowestRing,
        const FRingData& HighestRing,
        const FIntersectionRingData& IntersectionRing,
        int32 LeftIndex,
        int32 RightIndex,
        TArray<FVector>& OutSquareVertices
    );

    void GenerateAboveBelowIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData);

    void GenerateLateralTubeIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData);

    bool LineSegmentIntersectsTriangle(
        const FVector& LineStart,
        const FVector& LineEnd,
        const FVector& V0,
        const FVector& V1,
        const FVector& V2,
        FVector& OutIntersectionPoint
    );

    FVector ComputeCentroid(const TArray<FVector>& Vertices);
    TArray<FVector> OrderRingVertices(const TArray<FVector>& InputVertices);

    /**
     * Connects two rings of vertices to form triangles between them.
     * @param RingA Vertices of the first ring.
     * @param RingB Vertices of the second ring.
     * @param Vertices Output array of all vertices.
     * @param Triangles Output array of triangle indices.
     * @param BaseIndex Starting index for the vertices.
     */
    void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    void ConnectPartialRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * Connects an array of rings into a continuous tube.
     * @param Rings Array of ring data.
     * @param Vertices Output array of all vertices.
     * @param Triangles Output array of triangle indices.
     * @param BaseIndex Starting index for the vertices.
     */
    void ConnectRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    void ConnectPartialRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    void ConnectPartialRingArrayPaired(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    void ConnectIntersectionRingToSquare(FTwoTubeIntersectionData& TubeIntersectionData, int32& BaseIndex);

    void ConnectIntersectionCornerArrays(const TArray<TArray<FVector>>& CornerVertexArrays, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    void ConnectTwoTubeIntersection(FTwoTubeIntersectionData& TubeIntersectionData);

    /**
     * Constructs a cylindrical tube from a 2D profile and input parameters.
     * @param Profile 2D profile defining the cross-section.
     * @param StartPosition Starting position of the tube.
     * @param Direction Direction the tube is facing.
     * @param UpVector Up vector for orientation.
     * @param NumSegments Number of segments along the tube's length.
     * @param NumSides Number of sides per segment.
     * @param Length Length of the tube.
     * @param OutTube Output data structure to store the tube information.
     * @param ApplyBothAxis Whether to apply the profile curve to both axes.
     */
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

    bool IsPointInsideFrustum(const FRingData& StartRing, const FRingData& EndRing, const FVector& Point);

    bool IsPointInsidePolygon(const TArray<FVector>& RingVertices, const FVector& Point, const FVector& RingCenter);

    void RemoveInternalVertices(const FTubeData& TubeA, FTubeData& TubeB);

    void FindSegmentForPoint(
        const FVector& Point,
        const FTubeData& Tube,
        FRingData& StartRing,
        FRingData& EndRing
    );

    FVector CalculateCenterLinePoint(
        const FVector& Point,
        const FRingData& StartRing,
        const FRingData& EndRing
    );

    float InterpolatedRingRadius(
        const FVector& CenterlinePoint,
        const FRingData& StartRing,
        const FRingData& EndRing
    );

    void FindClosestVertex(const FVector& InputVertex, const TArray<FVector>& Vertices, FVector& OutVertex, int32& OutIndex);

    FVector FindIntersectionOnRing(const TArray<FVector>& RingVertices, const FVector& Direction, const FVector& Center);

    FVector FindIntersectionOnNewRing(const TArray<FVector>& RingVertices, const FVector& Direction, const FVector& Center, const FVector& TargetVertex);

    void FindIntersectionRingCardinalPoints(FIntersectionRingData& IntersectionRing, const FVector& StartCenter, const FVector& EndCenter);

    int32 FindVertexIndex(const TArray<FVector>& Vertices, const FVector& TargetVertex);

    int32 FindRingIndexByCenter(const TArray<FRingData>& Rings, const FVector& Center, float Tolerance = KINDA_SMALL_NUMBER);

    void OrderSquareIntersectionConnections(FTwoTubeIntersectionData& TubeIntersectionData);

    void OrderSquareIntersectionConnectionsOneCorner(
        const FTubeData& MainTubeIntersectionRings,
        const FTubeData& LateralTube,
        const TArray<FVector> RingVertices,
        const TArray<FVector> SquareVertices,
        const int32 StartIndex,
        TArray<TArray<FVector>>& OutVertexArrays,
        const bool Reversed
    );

    void RemoveVerticesByInterpolatedDirections(
        FTwoTubeIntersectionData& TubeIntersectionData,
        FIntersectionSquareData& IntersectionSquare
    );

    void RemoveVerticesInsideSquareByAngle(
        FTwoTubeIntersectionData& TubeIntersectionData,
        FIntersectionSquareData& IntersectionSquare
    );

    void RemoveDuplicateVertices(TArray<FVector>& Vertices, float Tolerance = 0.01f);

    void CalculateSquareCenter(FIntersectionSquareData& IntersectionSquare);

    void GetSquareAngles(FIntersectionSquareData& IntersectionSquare);

    bool IsVertexInsideSquareAngle(const FVector& Vertex, const FIntersectionSquareData& IntersectionSquare);

    bool IsPointOnLine(const FVector& LineStart, const FVector& LineEnd, const FVector& Point, float Tolerance = KINDA_SMALL_NUMBER);

    /**
     * Generates the intersection curve of two cylindrical shapes based on their 2D profiles and transformations.
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