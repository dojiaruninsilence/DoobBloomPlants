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

    /**
     * Connects an array of rings into a continuous tube.
     * @param Rings Array of ring data.
     * @param Vertices Output array of all vertices.
     * @param Triangles Output array of triangle indices.
     * @param BaseIndex Starting index for the vertices.
     */
    void ConnectRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);


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

    // ------------------------------------------------------ Planes and Lines Intersections ------------------------------------------------------------ //

    /**
     * Represents a plane equation defined by a normal vector and a scalar D.
     */
    struct FPlaneEquation {
        FVector Normal; ///< The normal vector of the plane.
        float D; ///< The scalar D in the plane equation.

        /**
         * Default constructor. Initializes a zero plane.
         */
        FPlaneEquation();

        /**
         * Constructs a plane from a normal vector and a scalar D.
         * @param InNormal The normal vector of the plane.
         * @param InD The scalar D in the plane equation.
         */
        FPlaneEquation(const FVector& InNormal, float InD);

        // possibly add methods for point containment or intersection checks here
    };

    /**
     * Calculates a plane equation using two rings.
     * @param RingA The first ring.
     * @param RingB The second ring.
     * @return The calculated plane equation.
     */
    FPlaneEquation CalculatePlaneFromRings(const FRingData& RingA, const FRingData& RingB);

    /**
     * Calculates the intersection of a ray with a plane.
     * @param Plane The plane equation.
     * @param RayOrigin The origin of the ray.
     * @param RayDirection The direction of the ray.
     * @param OutIntersectionPoint The calculated intersection point, if any.
     * @return True if the intersection is found; false otherwise.
     */
    bool CalculateIntersectionWithPlane(
        const FPlaneEquation& Plane,
        const FVector& RayOrigin,
        const FVector& RayDirection,
        FVector& OutIntersectionPoint
    );

    /**
     * Checks if a point lies within or on the boundary of a circle.
     * @param Point The point to check.
     * @param Center The center of the circle.
     * @param Radius The radius of the circle.
     * @return True if the point is inside or on the circle; false otherwise.
     */
    bool PointInsideCircle(const FVector& Point, const FVector& Center, float Radius);

    /**
     * Filters and identifies planes and lines that potentially intersect between two tubes.
     *
     * This function determines the relevant planes and lines in each tube that are within
     * the range of potential intersections. It outputs indices of the intersecting planes
     * and lines for both TubeA and TubeB. The filtering accounts for the geometry of the
     * tubes, including their positions, radii, and alignment.
     *
     * @param TubeA The first tube, containing its rings, center positions, and radii.
     * @param TubeB The second tube, containing its rings, center positions, and radii.
     * @param OutPlaneIndicesA A list of index pairs representing planes in TubeA that
     *                         are near potential intersection points. Each pair indicates
     *                         a plane defined by two rings in TubeA.
     * @param OutPlaneIndicesB A list of index pairs representing planes in TubeB that
     *                         are near potential intersection points. Each pair indicates
     *                         a plane defined by two rings in TubeB.
     * @param OutLineIndicesA A list of indices representing lines in TubeA that are
     *                        near potential intersection planes of TubeB.
     * @param OutLineIndicesB A list of indices representing lines in TubeB that are
     *                        near potential intersection planes of TubeA.
 */
    void FilterPlanesAndLines(
        const FTubeData& TubeA,
        const FTubeData& TubeB,
        TArray<TPair<int32, int32>>& OutPlaneIndicesA,
        TArray<TPair<int32, int32>>& OutPlaneIndicesB,
        TArray<int32>& OutLineIndicesA,
        TArray<int32>& OutLineIndicesB
    );

    // ------------------------------------------------------ Shape Intersections ------------------------------------------------------------ //

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