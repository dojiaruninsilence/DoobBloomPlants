/**
 * @file DoobGeometryUtils.h
 * @brief Utility functions and data structures for geometric operations in the Doob procedural generation system.
 *
 * This file provides tools for creating, manipulating, and analyzing geometric entities such as rings, tubes,
 * intersections, and polygons. It includes:
 * - **Data Structures**: Representing geometric entities (e.g., rings, tubes, intersections).
 * - **Utility Functions**: Operations like generating procedural shapes, connecting geometry, and calculating intersections.
 * - **Transformations and Interpolations**: Tools for manipulating and aligning geometric data.
 *
 * ## Purpose
 * The `DoobGeometryUtils` module is part of the Doob plant-growing game system and serves as the foundation for
 * procedural geometry generation. It provides reusable components for 3D geometry construction and optimization.
 *
 * ## Table of Contents
 * ### 1. Data Structures
 * - FRingData
 * - FIntersectionRingData
 * - FIntersectionSquareData
 * - FTubeData
 * - FTwoTubeIntersectionData
 *
 * ### 2. Ring Operations
 * - GenerateRing
 * - OrderRingVertices
 * - FindIntersectionOnRing
 * - FindIntersectionRingCardinalPoints
 *
 * ### 3. Tube Operations
 * - ConstructTubeFromProfile
 * - RemoveInternalVertices
 * - RemoveVerticesByInterpolatedDirections
 * - RemoveVerticesInsideSquareByAngle
 *
 * ### 4. Intersection Calculations
 * - GenerateIntersectionRing
 * - GenerateHalfIntersectionRing
 * - GenerateSquareAroundIntersection
 * - GenerateAboveBelowIntersectionRings
 * - GenerateLateralTubeIntersectionRings
 * - ConnectTwoTubeIntersection
 *
 * ### 5. Polygon and Shape Utilities
 * - IsPointInsidePolygon
 * - IsPointInsideFrustum
 * - RemoveDuplicateVertices
 *
 * ### 6. Connection Utilities
 * - ConnectRings
 * - ConnectPartialRings
 * - ConnectRingArray
 * - ConnectPartialRingArray
 * - ConnectPartialRingArrayPaired
 * - ConnectIntersectionRingToSquare
 * - ConnectIntersectionCornerArrays
 *
 * ### 7. Vertex and Index Utilities
 * - FindClosestVertex
 * - FindVertexIndex
 * - FindRingIndexByCenter
 *
 * ### 8. Transformations and Interpolations
 * - CalculateCenterLinePoint
 * - InterpolatedRingRadius
 * - CalculateSquareCenter
 * - GetSquareAngles
 *
 * ### 9. Utility Functions
 * - RemoveVerticesByInterpolatedDirections
 * - RemoveVerticesInsideSquareByAngle
 * - OrderSquareIntersectionConnections
 * - OrderSquareIntersectionConnectionsOneCorner
 *
 * ### 10. Geometry Query Functions
 * - IsPointOnLine
 * - LineSegmentIntersectsTriangle
 *
 * ### 11. Profiling and Curve Utilities
 * - GenerateIntersectionCurve
 *
 * ## Usage
 * Include this header file to use the `DoobGeometryUtils` utilities:
 * @code
 * #include "DoobGeometryUtils.h"
 * @endcode
 *
 * @note Ensure you link any dependencies required by this module when compiling.
 */


#pragma once

#include "CoreMinimal.h"

#include "DoobProfileUtils.h"

/**
 * @namespace DoobGeometryUtils
 * Contains utility functions and data structures for geometric operations.
 */
namespace DoobGeometryUtils {

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                                 1. Data Structures                                                                   //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- // 
    /**
     * @struct FRingData
     * Represents a ring in 3D space, including its vertices, orientation, and other geometric properties.
     *
     * This structure is used in the DoobGeometryUtils module for procedural geometry generation and manipulation.
     */
    struct FRingData {
        TArray<FVector> Vertices; ///< The vertices that define the ring in 3D space.
        FVector Normal; ///< The normal vector of the ring, indicating its orientation.
        FVector Center; ///< The center point of the ring in 3D space.
        FVector Direction; ///< The primary direction the ring is facing.
        FVector UpVector; ///< The up vector for the ring, used for alignment.
        float Radius; ///< The radius of the ring, calculated from the center to the vertices.
        bool bIsClosed = true; ///< Whether the ring forms a closed loop.
        int32 NumVertices = Vertices.Num(); ///< The number of vertices forming the ring.
        FName RingID; ///< A unique identifier for the ring, useful for tracking and debugging.
        bool bIsComplete = true; ///< Indicates whether the ring is complete (contains all expected vertices).
    };

    /**
     * @struct FIntersectionRingData
     * Stores data related to an intersection ring formed by the connection of two tubes in 3D space.
     *
     * This structure is used to represent and manage the vertices and properties of the intersection region
     * where a lateral tube meets a main tube.
     */
    struct FIntersectionRingData {
        TArray<FVector> CombinedVertices; ///< All vertices forming the intersection ring.
        TArray<FVector> MainTubeVertices; ///< The vertices contributed by the main tube.
        TArray<FVector> LateralTubeVertices; ///< The vertices contributed by the lateral tube.
        FVector Centroid; ///< The geometric center of the intersection ring.
        TArray<FVector> CardinalVertices; ///< The cardinal vertices (e.g., North, East, South, West) in the intersection ring.
        TArray<int32> CardinalIndices; ///< Indices of the cardinal vertices in the vertex array.
    };

    /**
     * @struct FIntersectionSquareData
     * Represents data for an intersection square, including corner points, side connections,
     * and partial rings associated with each side of the square.
     *
     * This structure is used to manage the geometry and properties of a square-shaped intersection
     * formed by connecting multiple tubes or ring structures.
     */
    struct FIntersectionSquareData {
        TArray<FVector> Corners; ///< The four corners of the intersection square.
        TArray<FVector> LeftSideRingConnections; ///< Connection vertices along the left side of the square.
        TArray<FVector> RightSideRingConnections; ///< Connection vertices along the right side of the square.
        FVector Center; ///< The geometric center of the intersection square.
        TArray<float> Angles; ///< Angles between adjacent corners, in degrees or radians.
        TArray<TArray<FVector>> BottomLeftPartialRings; ///< Partial ring vertices associated with the bottom-left corner.
        TArray<TArray<FVector>> BottomRightPartialRings; ///< Partial ring vertices associated with the bottom-right corner.
        TArray<TArray<FVector>> TopRightPartialRings; ///< Partial ring vertices associated with the top-right corner.
        TArray<TArray<FVector>> TopLeftPartialRings; ///< Partial ring vertices associated with the top-left corner.
    };

    /**
     * @struct FTubeData
     * Represents data for a 3D tube constructed from rings and a 2D profile.
     *
     * This structure defines the geometry, orientation, and properties of a tube, including its
     * cross-sectional profile, rings, and spatial configuration. Tubes are built from a sequence
     * of rings, defined by a 2D profile, and are used as a core element in constructing 3D models.
     */
    struct FTubeData {
        DoobProfileUtils::F2DProfile Profile;  ///< The 2D profile defining the tube's cross-section shape.
        TArray<FRingData> Rings; ///< Array of rings forming the tube, representing its geometry.
        FVector StartPosition; ///< The starting position of the tube in 3D space.
        FVector EndPosition; ///< The ending position of the tube in 3D space.
        FVector Direction; ///< The directional vector indicating the tube's orientation.
        FVector UpVector; ///< The up vector for the tube, used for orientation calculations.
        float Length; ///< The total length of the tube.
        int32 NumSegments; ///< The number of segments that make up the tube.
        int32 NumSides; ///< The number of sides per segment, determining the tube's smoothness.
        int32 NumRings = Rings.Num(); ///< The total number of rings forming the tube (derived from Rings).
        FName TubeID; ///< A unique identifier for the tube, used for referencing and tracking.
    };

    /**
     * @struct FTwoTubeIntersectionData
     * Represents data for the intersection of two 3D tubes.
     *
     * This structure encapsulates all relevant information about the intersection between a main
     * tube and a lateral tube, including geometry, partial rings, and intersection details. It is
     * used to compute and store the intersection's properties, as well as to manage the affected
     * geometry of both tubes.
     */
    struct FTwoTubeIntersectionData {
        FTubeData MainTube; ///< Data representing the main tube in the intersection.
        FTubeData LateralTube; ///< Data representing the lateral tube in the intersection.
        FIntersectionRingData IntersectionRing; ///< The ring formed by the intersection of the two tubes.
        FIntersectionSquareData IntersectionSquare; ///< Data representing a square approximation of the intersection.
        FRingData MTAboveIntersectionRing; ///< Ring of the main tube above the intersection.
        FRingData MTBelowIntersectionRing; ///< Ring of the main tube below the intersection.
        FRingData MTAboveIntersectionRingPartial; ///< Partial ring of the main tube above the intersection.
        FRingData MTBelowIntersectionRingPartial; ///< Partial ring of the main tube below the intersection.
        FTubeData LateralTubeIntersectionRings; ///< Tube data for the rings in the lateral tube around the intersection.
        FRingData LateralTubeFirstFullRing; ///< The first complete ring in the lateral tube after the intersection.
        FTubeData LateralTubeRemovedVertices; ///< Data representing vertices removed from the lateral tube during intersection processing.
        TArray<FVector> AllVertices; ///< Combined array of vertices from all relevant intersection components.
        TArray<int32> Triangles; ///< Triangle indices defining the mesh for the intersection geometry.
        TArray<FRingData> MainTubePartialRings; ///< Partial rings of the main tube near the intersection.
    };


    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                                 2. Ring Operations                                                                   //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Generates a ring of vertices in 3D space based on the provided parameters.
     *
     * This function creates a ring of vertices around a specified center point, aligned along a given
     * direction, with a defined radius and number of sides. It also calculates and stores additional
     * ring properties, such as the normal vector, center, and direction.
     *
     * @param Center The center point of the ring.
     * @param Direction The direction vector the ring is facing.
     * @param UpVector The up vector used for orientation and cross-product calculations.
     * @param Radius The radius of the ring.
     * @param NumSides The number of vertices forming the ring.
     * @param RingData Reference to an FRingData structure where the generated ring data will be stored.
     */
    void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, FRingData& RingData);

    /**
     * @brief Orders the vertices of a ring in counter-clockwise order based on their angular position.
     *
     * This function computes the centroid of the input vertices and calculates the average normal of the plane
     * formed by the vertices. It then projects the vertices onto the plane and orders them based on the angle of each
     * vertex relative to the centroid. The result is a list of vertices sorted by their angular position around the centroid.
     *
     * @param InputVertices A list of FVector points representing the vertices of the ring to be ordered.
     *
     * @return TArray<FVector> A sorted array of FVector points representing the ordered vertices of the ring.
     *
     * @note The sorting is based on the angle of the projected vertices in the plane perpendicular to the average normal.
     */
    TArray<FVector> OrderRingVertices(const TArray<FVector>& InputVertices);

    /**
     * @brief Finds and assigns the cardinal points of the intersection ring.
     *
     * This function calculates the cardinal points (North, East, South, West) of a ring formed by
     * intersection vertices based on their projection in relation to the ring's center and axis.
     * It calculates the highest, lowest, left, and right points based on the positions of the vertices
     * relative to the given start and end centers, and then stores these points in the provided
     * FIntersectionRingData structure.
     *
     * @param IntersectionRing The intersection ring data structure that contains the combined vertices.
     *                         The cardinal points will be assigned to this structure.
     * @param StartCenter The center position of the starting point for the intersection ring.
     * @param EndCenter The center position of the ending point for the intersection ring.
     *
     * @return void
     */
    void FindIntersectionRingCardinalPoints(FIntersectionRingData& IntersectionRing, const FVector& StartCenter, const FVector& EndCenter);

    /**
     * @brief Finds the index of a ring in the array that has a center point within a given tolerance.
     *
     * This function iterates through an array of ring data structures, compares the distance squared
     * between the provided center and each ring's center, and returns the index of the ring that
     * matches the provided center point within the specified tolerance.
     * If no match is found, it returns `INDEX_NONE` (-1).
     *
     * @param Rings The array of FRingData structures that contain ring information.
     * @param Center The center point to compare against the rings' centers.
     * @param Tolerance The tolerance value to allow for small differences in the center position (default is `KINDA_SMALL_NUMBER`).
     *
     * @return The index of the ring in the array that matches the center within the given tolerance,
     *         or `INDEX_NONE` (-1) if no matching ring is found.
     */
    int32 FindRingIndexByCenter(const TArray<FRingData>& Rings, const FVector& Center, float Tolerance = KINDA_SMALL_NUMBER);

    /**
     * @brief Removes vertices from the top and bottom intersection rings of two tubes based on interpolated directions.
     *
     * This function processes the top and bottom intersection rings of two tubes, removes vertices from them based on their
     * direction relative to an interpolated direction from two square corner directions, and reorders the remaining vertices.
     * The function calculates and checks intersections between the vertices and interpolated direction vectors,
     * adjusting the rings accordingly.
     *
     * @param TubeIntersectionData Data structure containing the intersection information of the two tubes, including the
     *        top and bottom intersection rings.
     * @param IntersectionSquare Data structure containing the corners of the intersection square used to compute the direction
     *        vectors for vertex removal.
     *
     * @warning The function assumes the IntersectionSquare has exactly 4 corners. If the number of corners is not 4, a warning
     *          is logged and the function returns early without performing any modifications.
     */
    void RemoveVerticesByInterpolatedDirections(
        FTwoTubeIntersectionData& TubeIntersectionData,
        FIntersectionSquareData& IntersectionSquare
    );

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                                 3. Tube Operations                                                                   //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Constructs a tube mesh from a 2D profile.
     *
     * This function generates a tube mesh in 3D space based on a 2D profile. The tube is constructed starting from a specified position and direction, with each segment of the tube derived from the profile's points. The tube's direction and up vector determine its orientation, and the number of segments and sides control the tube's resolution. The final tube is stored in the `OutTube` parameter, which will include all necessary data such as the rings, direction, position, and length.
     *
     * @param Profile A 2D profile that defines the shape of the tube. It provides a set of points that will be used to generate the cross-sectional rings of the tube.
     * @param StartPosition The starting position for the tube, where the first ring will be placed.
     * @param Direction The direction vector for the tube, which determines the orientation of the tube in 3D space.
     * @param UpVector The up vector used to orient the rings in space, ensuring the correct orientation for each ring.
     * @param NumSegments The number of segments into which the tube should be divided along its length.
     * @param NumSides The number of sides to use for each ring of the tube, determining the resolution of the tube's cross-section.
     * @param Length The total length of the tube.
     * @param OutTube The output tube data that will store the generated rings, direction, up vector, and other related information.
     * @param ApplyBothAxis If true, applies the curve to both the X and Y axes based on the profile, allowing for more complex tube shapes. If false, the curve is applied only to one axis.
     *
     * @note The tube's length and segment placement are adjusted based on the profile's shape, with special handling for the Y values in the profile when `ApplyBothAxis` is enabled.
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

    /**
     * @brief Removes vertices from TubeB that are inside the frustums defined by consecutive rings in TubeA.
     *
     * This function processes the vertices in TubeB and eliminates those that lie within the frustums formed
     * by consecutive rings in TubeA. The resulting TubeB contains only vertices that are outside all frustums
     * in TubeA.
     *
     * @param TubeA The reference tube whose frustums are used to test for vertex inclusion.
     * @param TubeB The tube from which internal vertices are removed. This tube is modified in place.
     *
     * @details The function performs the following steps:
     * - Iterates through each ring in TubeB.
     * - For each vertex in a ring of TubeB, checks if it is inside any frustum defined by consecutive rings in TubeA.
     * - If a vertex is outside all frustums, it is retained; otherwise, it is discarded.
     * - Constructs a temporary tube with the retained vertices and updates TubeB with this tube.
     *
     * @note If all vertices of a ring in TubeB are removed, that ring is excluded from the resulting TubeB.
     * The function ensures that ring properties such as center, direction, and radius are preserved for retained rings.
     */
    void RemoveInternalVertices(const FTubeData& TubeA, FTubeData& TubeB);

    /**
     * @brief Finds the segment in a tube where a given point lies and assigns the corresponding start and end rings.
     *
     * This function iterates through the rings of a tube to determine which segment the given point lies within.
     * It identifies the pair of consecutive rings (start and end) that define the segment containing the point.
     *
     * @param Point The 3D position of the point to check.
     * @param Tube The tube structure (FTubeData) containing rings that define the tube geometry.
     * @param StartRing Output parameter representing the ring at the start of the segment containing the point.
     * @param EndRing Output parameter representing the ring at the end of the segment containing the point.
     *
     * @details The function calculates the projection of the point onto the vector between each pair of consecutive rings
     * and checks whether the projection lies within the bounds of the segment. If the point lies within a segment, the
     * corresponding start and end rings are assigned to the output parameters.
     *
     * @note The Tube must contain at least two rings for the function to perform correctly.
     */
    void FindSegmentForPoint(
        const FVector& Point,
        const FTubeData& Tube,
        FRingData& StartRing,
        FRingData& EndRing
    );

    /**
     * @brief Removes vertices inside a defined square's angular area from the rings of a tube.
     *
     * This function processes the rings of a main tube and removes any vertices that lie inside the
     * angular area defined by the square. It also modifies the rings of the tube and updates the
     * intersection square data, adding connection points at the sides of the square.
     *
     * @param TubeIntersectionData The data structure containing the two tubes' intersection
     *                             information, including the rings of the main tube.
     * @param IntersectionSquare   The data structure representing the intersection square, which
     *                             contains the square's corners and angular data.
     *
     * @warning The function assumes the intersection square contains exactly 4 vertices. If it does
     *          not, a warning message will be logged, and the function will return without making
     *          changes.
     *
     * @note This function modifies the provided TubeIntersectionData by removing vertices from the
     *       rings and adding new connection points to the sides of the square.
     */
    void RemoveVerticesInsideSquareByAngle(
        FTwoTubeIntersectionData& TubeIntersectionData,
        FIntersectionSquareData& IntersectionSquare
    );

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                            4. Intersection Calculations                                                              //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Generates an intersection ring between two tubes (main and lateral).
     *
     * This function calculates the vertices that form the intersection ring between two tubes: the main tube and the lateral tube.
     * The vertices from both tubes are combined, duplicate vertices are removed, and the resulting ring is ordered.
     * Additionally, the centroid of the intersection ring is computed.
     *
     * @param MainTube The main tube in the intersection. This tube serves as the reference for the first half of the ring.
     * @param LateralTube The lateral tube in the intersection. This tube serves as the reference for the second half of the ring.
     * @param OutRingData The resulting intersection ring data, including the vertices and the centroid of the ring.
     * @param Precision A threshold value for determining the precision of the intersection (default value: KINDA_SMALL_NUMBER).
     *
     * @note The function relies on the `GenerateHalfIntersectionRing` function to compute the half-rings for both the main and lateral tubes.
     *       The final intersection ring is built by combining the results from both tubes and eliminating duplicates.
     */
    void GenerateIntersectionRing(
        const FTubeData& MainTube,
        const FTubeData& LateralTube,
        FIntersectionRingData& OutRing,
        float Precision = KINDA_SMALL_NUMBER
    );

    /**
     * @function GenerateHalfIntersectionRing
     * @brief Generates half of the intersection ring between two tubes by calculating intersection points
     *        between the line segments of TubeA and the rectangles of TubeB.
     *
     * This function iterates through the rings of TubeA and checks for intersection points with the
     * rectangles of TubeB. For each line segment of TubeA, it finds the intersection points with the
     * triangles formed by consecutive vertices of TubeB. The resulting intersection points are stored
     * in the output array, OutRingVertices.
     *
     * @param TubeA The first tube (MainTube) involved in the intersection calculation.
     * @param TubeB The second tube (LateralTube) involved in the intersection calculation.
     * @param OutRingVertices The output array where the intersection points will be stored.
     * @param Precision The precision value used for the intersection calculation (default: KINDA_SMALL_NUMBER).
     */
    void GenerateHalfIntersectionRing(
        const FTubeData& MainTube,
        const FTubeData& LateralTube,
        TArray<FVector>& OutRingData,
        float Precision = KINDA_SMALL_NUMBER
    );

    /**
     * @brief Generates a square around the intersection between two rings based on the intersection ring.
     *
     * This function calculates the intersection points and creates a square that lies between the lowest and highest rings.
     * The square is created using the intersection ring and the specified left and right indices to determine the
     * positions of the square's vertices. The output square vertices are stored in `OutSquareVertices`.
     *
     * @param LowestRing The lowest ring in the tube to be used for intersection.
     * @param HighestRing The highest ring in the tube to be used for intersection.
     * @param IntersectionRing The ring containing the intersection vertices.
     * @param LeftIndex The index of the left intersection point in the intersection ring.
     * @param RightIndex The index of the right intersection point in the intersection ring.
     * @param OutSquareVertices The output array where the vertices of the generated square will be stored.
     *
     * @warning If any of the input arrays (`CombinedVertices`, `LowestRing.Vertices`, or `HighestRing.Vertices`) are empty,
     *          a warning is logged, and the function exits early without generating the square.
     */
    void GenerateSquareAroundIntersection(
        const FRingData& LowestRing,
        const FRingData& HighestRing,
        const FIntersectionRingData& IntersectionRing,
        int32 LeftIndex,
        int32 RightIndex,
        TArray<FVector>& OutSquareVertices
    );

    /**
     * @brief Generates the rings above and below the intersection for a two-tube intersection.
     *
     * This function calculates the rings positioned above and below the intersection points of the main tube
     * and lateral tube at the cardinal vertices. The rings are generated by finding the segments in the main tube
     * that correspond to the lowest and highest intersection points, calculating their respective radii,
     * and then generating the above and below intersection rings.
     *
     * @param TubeIntersectionData [in, out] The data structure containing the main tube, lateral tube, and intersection ring data.
     * The generated above and below intersection rings will be stored in the `MTAboveIntersectionRing` and `MTBelowIntersectionRing` members of the structure.
     */
    void GenerateAboveBelowIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData);

    /**
     * @brief Generates lateral tube intersection rings for a given tube intersection data.
     *
     * This function generates the intersection rings for a lateral tube that intersects with a main tube.
     * It calculates the intersection rings at specific points along the lateral tube and stores them in the
     * TubeIntersectionData. The process includes finding the nearest vertices, generating rings, and adjusting
     * the vertices of the intersection rings based on the calculated centers and radii.
     *
     * @param TubeIntersectionData The data structure containing the details of the intersection between the
     *        main tube and the lateral tube. This includes information about the intersection ring, the main tube,
     *        and the lateral tube.
     *
     * @note This function assumes that the TubeIntersectionData contains valid and properly initialized
     *       data for the intersection and tubes.
     */
    void GenerateLateralTubeIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData);

    /**
     * @brief Finds the intersection point of a line (defined by a direction vector) with a ring.
     *
     * This function calculates where a line, originating from a given point and extending in a specific direction,
     * intersects with a polygonal ring (defined by a list of vertices). The ring is assumed to be closed, meaning
     * that the last vertex connects back to the first vertex. The intersection point is computed using the parametric
     * equation of a line and the properties of the ring's segments.
     *
     * @param RingVertices The vertices that define the ring.
     * @param Direction The direction vector of the line in which to search for the intersection.
     * @param Center The center point of the ring.
     *
     * @return The intersection point, or FVector::ZeroVector if no intersection is found.
     *
     * @note This function assumes that the ring is closed and that the input direction is not parallel to the
     *       edges of the ring.
     */
    FVector FindIntersectionOnRing(const TArray<FVector>& RingVertices, const FVector& Direction, const FVector& Center);

    /**
     * @brief Finds the intersection points between a direction and a ring, and returns the closest intersection to a target vertex.
     *
     * This function iterates over the vertices of the provided ring and checks for intersections between the given
     * direction and the line segments forming the ring. If any intersections are found, it returns the intersection
     * point that is closest to a specified target vertex.
     *
     * @param RingVertices An array of FVector representing the vertices of the ring.
     * @param Direction The direction vector to check for intersections with the ring.
     * @param Center The center of the ring used in calculations.
     * @param TargetVertex The vertex used as a reference to find the closest intersection point.
     *
     * @return The FVector representing the closest intersection point to the target vertex, or FVector::ZeroVector if no intersections are found.
     *
     * @note The function handles wrapping of the ring vertices and checks for parallel lines to avoid division by zero.
     */
    FVector FindIntersectionOnNewRing(const TArray<FVector>& RingVertices, const FVector& Direction, const FVector& Center, const FVector& TargetVertex);

    /**
     * @brief Generates the intersection curve between two profiles.
     *
     * This function calculates the intersection curve between two 2D profiles (main and lateral) by evaluating them along
     * their respective parametric segments and checking the distance between points. If the distance between points from the
     * main and lateral profiles is within a given threshold, the midpoint is added to the intersection curve.
     *
     * @param MainProfile The 2D profile for the main cylinder.
     * @param MainTransform The transformation applied to the main profile.
     * @param LateralProfile The 2D profile for the lateral cylinder.
     * @param LateralTransform The transformation applied to the lateral profile.
     * @param MainSegments The number of segments to divide the main profile into.
     * @param LateralSegments The number of segments to divide the lateral profile into.
     * @param Threshold The maximum distance between points from the main and lateral profiles for them to be considered part of the intersection curve.
     *
     * @return A TArray of FVector points representing the intersection curve between the main and lateral profiles.
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

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                           5. Polygon and Shape Utilities                                                             //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Checks if a line segment intersects a triangle.
     *
     * This function performs the Möller–Trumbore intersection algorithm to determine if a
     * line segment intersects a triangle in 3D space. If an intersection occurs, the point
     * of intersection is calculated and returned in the OutIntersectionPoint.
     *
     * @param LineStart The starting point of the line segment.
     * @param LineEnd The ending point of the line segment.
     * @param V0 The first vertex of the triangle.
     * @param V1 The second vertex of the triangle.
     * @param V2 The third vertex of the triangle.
     * @param OutIntersectionPoint The point of intersection (if one exists) between the line segment and the triangle.
     *
     * @return true if the line segment intersects the triangle, false otherwise.
     *
     * @note The function uses floating-point arithmetic, so it may not be perfectly precise
     * due to potential rounding errors in some edge cases.
     */
    bool LineSegmentIntersectsTriangle(
        const FVector& LineStart,
        const FVector& LineEnd,
        const FVector& V0,
        const FVector& V1,
        const FVector& V2,
        FVector& OutIntersectionPoint
    );

    /**
     * @brief Computes the centroid (geometric center) of a set of vertices.
     *
     * This function calculates the centroid of a polygon or shape by averaging the positions
     * of all the vertices provided. The centroid is the point at which the shape would balance
     * if it were a physical object.
     *
     * @param Vertices The array of vertices (FVector) whose centroid is to be computed.
     * @return FVector The computed centroid of the vertices.
     */
    FVector ComputeCentroid(const TArray<FVector>& Vertices);

    /**
     * @brief Calculates the center point of a square defined by its corner vertices.
     *
     * This function computes the center of a square by averaging the positions
     * of its four corner vertices. The resulting center point is stored in the
     * `Center` attribute of the provided `FIntersectionSquareData` object.
     *
     * @param[in,out] IntersectionSquare The square whose center is being calculated.
     *                                   The calculated center is stored in the `Center` field.
     *
     * @note The square must have exactly four vertices in its `Corners` array.
     */
    void CalculateSquareCenter(FIntersectionSquareData& IntersectionSquare);

    /**
     * @brief Calculates the angles of the corners of a square relative to its center.
     *
     * This function calculates the angles of the square's corners in the 2D plane after transforming
     * the square into local space where the normal vector of the square defines the Z-axis.
     * The angles are calculated in relation to the center of the square, and are sorted in ascending order.
     *
     * @param IntersectionSquare Reference to the FIntersectionSquareData object containing the square's corners
     *        and center. The calculated angles will be stored in the `Angles` property of this object.
     */
    void GetSquareAngles(FIntersectionSquareData& IntersectionSquare);

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                               6. Connection Utilities                                                                //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Connects two rings of vertices by creating triangles between them.
     *
     * This function connects two rings (RingA and RingB) by creating a series of triangles
     * between the corresponding vertices in the two rings. It adds the vertices and triangle
     * indices to the provided arrays `Vertices` and `Triangles`, respectively. The index
     * for the next set of vertices is updated through `BaseIndex`.
     *
     * @param RingA The first ring of vertices.
     * @param RingB The second ring of vertices.
     * @param Vertices An array that stores the combined vertices of RingA and RingB.
     * @param Triangles An array that stores the triangle indices formed by connecting RingA and RingB.
     * @param BaseIndex The starting index for adding vertices from RingA and RingB into the `Vertices` array.
     *
     * @note Assumes that RingA and RingB have the same number of vertices.
     */
    void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * @brief Connects two partially defined rings (RingA and RingB) by generating triangles between their vertices.
     *
     * This function connects two partial rings of vertices, `RingA` and `RingB`, and creates the corresponding triangles
     * to form a smooth transition between them. It appends the vertices of both rings to the `Vertices` array and updates
     * the `Triangles` array with indices that define the triangles formed by connecting the rings.
     *
     * @param RingA The first ring of vertices.
     * @param RingB The second ring of vertices.
     * @param Vertices The array to append the vertices of both rings to.
     * @param Triangles The array to append the triangle indices.
     * @param BaseIndex The base index used for computing the vertex indices for both rings.
     */
    void ConnectPartialRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * @brief Connects an array of rings by iterating through consecutive rings and generating triangles to connect them.
     *
     * This function iterates through an array of rings (`Rings`), connecting each consecutive pair of rings by generating
     * triangles between their vertices. It calls the `ConnectRings` function for each consecutive ring pair to append the
     * vertices and triangle indices to the `Vertices` and `Triangles` arrays, respectively. The `BaseIndex` is updated
     * during the connection of each pair of rings.
     *
     * @param Rings The array of ring data, where each element contains a ring's vertices to be connected.
     * @param Vertices The array to append the vertices of the connected rings to.
     * @param Triangles The array to append the triangle indices generated between consecutive rings.
     * @param BaseIndex The base index used for computing the vertex indices for each ring.
     */
    void ConnectRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * @brief Connects an array of partial rings by iterating through consecutive rings and generating triangles to connect them.
     *
     * This function iterates through an array of partial rings (`Rings`), connecting each consecutive pair of rings by generating
     * triangles between their vertices. It calls the `ConnectPartialRings` function for each consecutive pair of partial rings to append
     * the vertices and triangle indices to the `Vertices` and `Triangles` arrays, respectively. The `BaseIndex` is updated
     * during the connection of each pair of partial rings.
     *
     * @param Rings The array of partial ring data, where each element contains a ring's vertices to be connected.
     * @param Vertices The array to append the vertices of the connected partial rings to.
     * @param Triangles The array to append the triangle indices generated between consecutive partial rings.
     * @param BaseIndex The base index used for computing the vertex indices for each ring.
     */
    void ConnectPartialRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * @brief Connects an array of partial rings in pairs, generating triangles between consecutive rings in each pair.
     *
     * This function iterates through an array of partial rings (`Rings`), connecting each pair of consecutive rings (i.e., Ring[i] and Ring[i+1]) by generating
     * triangles between their vertices. It calls the `ConnectPartialRings` function for each pair of consecutive partial rings and appends the
     * vertices and triangle indices to the `Vertices` and `Triangles` arrays, respectively. The `BaseIndex` is updated during the connection of each pair of rings.
     * The function operates on consecutive rings in pairs, thus incrementing by 2 during the iteration.
     *
     * @param Rings The array of partial ring data, where each element contains a ring's vertices to be connected.
     * @param Vertices The array to append the vertices of the connected partial rings to.
     * @param Triangles The array to append the triangle indices generated between consecutive partial rings.
     * @param BaseIndex The base index used for computing the vertex indices for each ring.
     */
    void ConnectPartialRingArrayPaired(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * @brief Connects arrays of corner vertices from multiple intersection arrays, generating triangles between corresponding corners.
     *
     * This function takes an array of arrays of corner vertices (`CornerVertexArrays`) and connects each consecutive pair of vertex arrays
     * by generating triangles between their corresponding corners. The function iterates through the `CornerVertexArrays` and processes
     * each pair of consecutive corner arrays, creating two triangles for each pair of vertices from `CurrentVertexArray` and `NextVertexArray`.
     * The `Vertices` and `Triangles` arrays are populated with the resulting geometry, and `BaseIndex` is updated during the process.
     *
     * @param CornerVertexArrays The array of arrays containing corner vertex data for intersections, where each inner array represents a set of corner vertices.
     * @param Vertices The array to append the vertices from each intersection corner array to.
     * @param Triangles The array to append the triangle indices generated between corresponding corners of the consecutive arrays.
     * @param BaseIndex The base index used to compute the vertex indices for each set of corner arrays.
     */
    void ConnectIntersectionCornerArrays(const TArray<TArray<FVector>>& CornerVertexArrays, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);

    /**
     * @brief Connects intersection ring data to a square by generating triangles for each corner array of the intersection square.
     *
     * This function processes the four corner arrays of the intersection square, which are stored in the `TubeIntersectionData.IntersectionSquare` object.
     * For each corner (top right, bottom right, bottom left, and top left), the function calls `ConnectIntersectionCornerArrays` to connect the corner vertices
     * and generate triangles. The resulting vertices and triangles are appended to the `TubeIntersectionData.AllVertices` and `TubeIntersectionData.Triangles` arrays, respectively.
     * The `BaseIndex` is updated during the process to maintain the correct indexing for vertices and triangles.
     *
     * @param TubeIntersectionData The data structure containing intersection details, including corner ring data and the square’s geometry.
     * @param BaseIndex The base index used to compute the vertex indices when connecting corner arrays.
     */
    void ConnectIntersectionRingToSquare(FTwoTubeIntersectionData& TubeIntersectionData, int32& BaseIndex);

    /**
     * @brief Connects the two tubes at the intersection, generating the necessary geometry for both tubes' rings.
     *
     * This function handles the connection of two intersecting tubes by processing and connecting their rings. It classifies and organizes the rings into three categories:
     * below the intersection, partial intersection, and above the intersection. It also processes the lateral tube intersection rings and generates the corresponding
     * geometry. The function utilizes multiple helper functions like `ConnectRingArray`, `ConnectPartialRingArray`, and `ConnectIntersectionRingToSquare`
     * to generate vertices and triangles, appending them to the `TubeIntersectionData.AllVertices` and `TubeIntersectionData.Triangles` arrays.
     * Additionally, the function manages the addition of full and reversed lateral tube rings for further geometry generation.
     *
     * @param TubeIntersectionData The data structure that contains information about the main tube, intersection, lateral tube, and the vertices and triangles for the tubes' rings.
     */
    void ConnectTwoTubeIntersection(FTwoTubeIntersectionData& TubeIntersectionData);

    /**
     * Orders the intersection connections for a square shape at a specific corner, involving the intersection of a main tube and a lateral tube.
     * This function calculates the proper order of vertices along the intersection path, considering the relationship between the rings of the main tube
     * and the lateral tube, and adds the square's vertices as the final set.
     *
     * @param MainTubeIntersectionRings The intersection data for the main tube, containing the rings involved in the connection process.
     * @param LateralTube The data of the lateral tube intersecting with the main tube.
     * @param RingVertices The vertices of the current ring to start processing.
     * @param SquareVertices The vertices defining the square at the intersection.
     * @param StartIndex The index to start the connection process within the main tube's rings.
     * @param OutVertexArrays The resulting array of vertex arrays, where each array represents a set of vertices at a particular step of the process.
     * @param Reversed A boolean flag to determine the direction of traversal; if true, the process moves backward through the rings, otherwise it moves forward.
     *
     * @warning The function assumes that the `MainTubeIntersectionRings` contains at least one ring; otherwise, a warning will be logged.
     */
    void OrderSquareIntersectionConnectionsOneCorner(
        const FTubeData& MainTubeIntersectionRings,
        const FTubeData& LateralTube,
        const TArray<FVector> RingVertices,
        const TArray<FVector> SquareVertices,
        const int32 StartIndex,
        TArray<TArray<FVector>>& OutVertexArrays,
        const bool Reversed
    );

    /**
     * @brief Orders the intersection connections for the square intersection at the corner.
     *
     * This function calculates the intersection points for a square formed by two tubes.
     * It iterates through the vertices of the intersection ring and reorders them based on cardinal points (right, bottom, left, top).
     * For each segment of the square, it calculates the intersection on new rings, stores the results in temporary arrays,
     * and calls `OrderSquareIntersectionConnectionsOneCorner` for each corner of the intersection.
     *
     * @param TubeIntersectionData The data structure containing the intersection information of the tubes,
     *                              including the intersection rings, cardinal vertices, and related data for the lateral tube.
     *
     * @note The function uses helper functions such as `ReorderedArray`, `FindVertexIndex`, `FindIntersectionOnNewRing`, and others
     *       to handle calculations of direction, centerline points, and the final intersection points.
     *       It also ensures the intersection connections are ordered correctly for each corner of the square.
     */
    void OrderSquareIntersectionConnections(FTwoTubeIntersectionData& TubeIntersectionData);

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                            7. Vertex and Index Utilities                                                             //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Finds the closest vertex to a given input vertex from an array of vertices.
     *
     * This function iterates through an array of vertices and computes the squared Euclidean distance
     * between each vertex and the input vertex. The closest vertex is identified based on the smallest
     * distance, and its value and index are returned through the output parameters.
     *
     * @param InputVertex The vertex to compare against the array of vertices.
     * @param Vertices The array of vertices from which the closest one will be found.
     * @param OutVertex The closest vertex found. This will be set to the closest vertex on successful completion.
     * @param OutIndex The index of the closest vertex in the input array. This will be set to the index of the closest vertex.
     *
     * @note If the input vertices array is empty, the output vertex will be set to `FVector::ZeroVector`
     *       and the index will be set to `-1`.
     */
    void FindClosestVertex(const FVector& InputVertex, const TArray<FVector>& Vertices, FVector& OutVertex, int32& OutIndex);

    /**
     * @brief Finds the index of a specified vertex in the given array of vertices.
     *
     * This function searches through an array of vertices and returns the index of the vertex
     * that matches the provided target vertex. If the target vertex is not found, it returns -1.
     * The function compares each vertex in the array with the target vertex using the equality
     * operator, which compares their position in 3D space.
     *
     * @param Vertices The array of vertices to search through.
     * @param TargetVertex The vertex to search for in the array.
     *
     * @return The index of the target vertex in the `Vertices` array. If the vertex is not found,
     *         it returns -1.
     */
    int32 FindVertexIndex(const TArray<FVector>& Vertices, const FVector& TargetVertex);

    /**
     * @brief Removes duplicate vertices from an array of vertices within a given tolerance.
     *
     * This function iterates through the input array of vertices and checks each vertex against the others.
     * If two vertices are found to be within a specified tolerance, the second one is considered a duplicate
     * and will be excluded from the result. The input array is updated to contain only unique vertices.
     *
     * @param Vertices The array of vertices to be processed. This array will be modified to contain only unique vertices.
     * @param Tolerance The tolerance distance used to identify duplicate vertices. Vertices within this distance are considered duplicates. Default is 0.01f.
     */
    void RemoveDuplicateVertices(TArray<FVector>& Vertices, float Tolerance = 0.01f);

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                        8. Transformations and Interpolations                                                         //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Calculates the closest point on the centerline of a frustum segment defined by two rings.
     *
     * This function determines the projection of a given point onto the centerline
     * connecting the centers of the start and end rings of a frustum. The result
     * is clamped to lie within the segment between the two rings.
     *
     * @param Point The 3D point for which the closest centerline point is to be calculated.
     * @param StartRing The starting ring of the frustum segment, defining one end of the centerline.
     * @param EndRing The ending ring of the frustum segment, defining the other end of the centerline.
     * @return FVector The closest point on the centerline of the frustum segment.
     *
     * @details
     * - The function projects the vector from the start ring's center to the input point onto
     *   the centerline defined by the centers of the start and end rings.
     * - The projection is clamped to ensure the resulting point lies within the segment of the centerline
     *   bounded by the start and end rings.
     * - The function is useful in operations involving point alignment with a frustum's axis.
     */
    FVector CalculateCenterLinePoint(
        const FVector& Point,
        const FRingData& StartRing,
        const FRingData& EndRing
    );

    /**
     * @brief Computes the interpolated radius of a ring at a given centerline point between two rings.
     *
     * This function calculates the radius of a ring at the position of a given centerline point,
     * using linear interpolation between the radii of the start and end rings.
     *
     * @param CenterlinePoint The point along the centerline where the radius is to be interpolated.
     * @param StartRing The starting ring of the segment, providing its center and radius.
     * @param EndRing The ending ring of the segment, providing its center and radius.
     * @return The interpolated radius of the ring at the specified centerline point.
     *
     * @details
     * - If the distance between the start and end ring centers is zero, the function returns the average radius of the two rings.
     * - The interpolation factor is calculated as the relative distance of the centerline point from the start ring center,
     *   clamped between 0.0 and 1.0.
     */
    float InterpolatedRingRadius(
        const FVector& CenterlinePoint,
        const FRingData& StartRing,
        const FRingData& EndRing
    );

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                                9. Utility Functions                                                                  //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                            10. Geometry Query Functions                                                              //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * @brief Checks whether a point lies inside a frustum defined by two rings.
     *
     * This function determines if a given point is located within a frustum formed by two cross-sectional rings (StartRing and EndRing).
     * The frustum is defined by interpolating the radius and center position between the two rings along their connecting axis.
     *
     * @param StartRing The starting ring of the frustum. It defines one end of the frustum's axis and the radius at that end.
     * @param EndRing The ending ring of the frustum. It defines the other end of the frustum's axis and the radius at that end.
     * @param Point The point to check for containment within the frustum.
     * @return true If the point lies inside the frustum.
     * @return false If the point lies outside the frustum.
     *
     * @details
     * - The function first computes the height vector between the two rings and projects the point onto this vector.
     * - It then checks whether the projected point lies within the bounds of the frustum's height.
     * - The radius at the projected point is interpolated linearly between the radii of the two rings.
     * - If the distance of the point from the frustum's axis exceeds the interpolated radius, the point is outside the frustum.
     * - A final polygonal cross-section check is performed using the vertices of the appropriate ring (based on the height projection).
     *
     * @note The function assumes that the two rings are aligned and form a valid frustum with a smooth transition of radii between them.
     */
    bool IsPointInsideFrustum(const FRingData& StartRing, const FRingData& EndRing, const FVector& Point);

    /**
     * @brief Determines whether a given point lies inside a 3D polygon defined by its vertices.
     *
     * This function checks if a 3D point is within the bounds of a polygon, represented by a set of vertices.
     * It projects the vertices and the point into a local 2D plane, computes a polygon normal, and performs
     * a 2D point-in-polygon test.
     *
     * @param RingVertices The vertices of the polygon, defined in 3D space.
     * @param Point The 3D point to test for inclusion within the polygon.
     * @param RingCenter The center of the polygon, used to project the vertices and point into the local plane.
     * @return true If the point lies inside the polygon.
     * @return false If the point lies outside the polygon or if the polygon is invalid (fewer than 3 vertices).
     *
     * @note This function assumes the polygon is planar and closed. If the input vertices do not define a valid polygon,
     * the function returns false.
     */
    bool IsPointInsidePolygon(const TArray<FVector>& RingVertices, const FVector& Point, const FVector& RingCenter);

    /**
     * @brief Checks if a vertex lies within the angular range of a square.
     *
     * This function converts the vertex position to polar coordinates relative to the center of the square
     * and checks whether its angle falls within the sorted angular range of the square.
     *
     * @param Vertex The vertex to be checked.
     * @param IntersectionSquare The square data containing the center and angle information of the square.
     *
     * @return True if the vertex is within the angular range of the square, false otherwise.
     */
    bool IsVertexInsideSquareAngle(const FVector& Vertex, const FIntersectionSquareData& IntersectionSquare);

    /**
     * @brief Checks if a given point lies on a line segment defined by two points.
     *
     * This function determines whether a point lies on the line segment between two points, with a specified tolerance for precision.
     * It checks for collinearity and whether the point is within the bounds of the line segment.
     *
     * @param LineStart The starting point of the line segment.
     * @param LineEnd The ending point of the line segment.
     * @param Point The point to check.
     * @param Tolerance The tolerance for precision in the calculation. Default is KINDA_SMALL_NUMBER.
     * @return true if the point lies on the line segment, false otherwise.
     */
    bool IsPointOnLine(const FVector& LineStart, const FVector& LineEnd, const FVector& Point, float Tolerance = KINDA_SMALL_NUMBER);

    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
    //                                                          11. Profiling and Curve Utilities                                                           //
    // ---------------------------------------------------------------------------------------------------------------------------------------------------- //
}