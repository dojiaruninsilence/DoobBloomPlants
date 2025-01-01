#include "CoreMinimal.h"
#include "DoobGeometryUtils.h"
#include "Misc/AutomationTest.h"

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFilterPlanesAndLinesTest, "Doob.Geometry.FilterPlanesAndLines", EAutomationTestFlags::EditorContext | EAutomationTestFlags::EngineFilter)
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FCalculateIntersectionWithPlaneTest, "Doob.Geometry.CalculateIntersectionWithPlane", EAutomationTestFlags::EditorContext | EAutomationTestFlags::EngineFilter)

bool FFilterPlanesAndLinesTest::RunTest(const FString& Parameters) {
    // Define parameters for TubeA
    FVector StartPositionA(0, 0, 0);
    FVector DirectionA = FVector::ForwardVector;
    FVector UpVectorA = FVector::UpVector;
    float LengthA = 10.0f;
    int32 NumSegmentsA = 2;
    int32 NumSidesA = 16;
    float StartRadiusA = 50.0f;
    float EndRadiusA = 50.0f;

    // Generate profile and construct TubeA
    DoobProfileUtils::F2DProfile ProfileA = DoobProfileUtils::GenerateLinearProfile(NumSegmentsA, StartRadiusA, EndRadiusA);
    DoobGeometryUtils::FTubeData TubeA;
    DoobGeometryUtils::ConstructTubeFromProfile(ProfileA, StartPositionA, DirectionA, UpVectorA, NumSegmentsA, NumSidesA, LengthA, TubeA, false);

    // Define parameters for TubeB
    FVector StartPositionB(0, 25, 5);
    FVector DirectionB = FVector::ForwardVector;
    FVector UpVectorB = FVector::UpVector;
    float LengthB = 10.0f;
    int32 NumSegmentsB = 2;
    int32 NumSidesB = 16;
    float StartRadiusB = 25.0f;
    float EndRadiusB = 25.0f;

    // Generate profile and construct TubeB
    DoobProfileUtils::F2DProfile ProfileB = DoobProfileUtils::GenerateLinearProfile(NumSegmentsB, StartRadiusB, EndRadiusB);
    DoobGeometryUtils::FTubeData TubeB;
    DoobGeometryUtils::ConstructTubeFromProfile(ProfileB, StartPositionB, DirectionB, UpVectorB, NumSegmentsB, NumSidesB, LengthB, TubeB, false);

    // Log tube details for debugging
    UE_LOG(LogTemp, Log, TEXT("TubeA has %d rings"), TubeA.Rings.Num());
    for (int32 i = 0; i < TubeA.Rings.Num(); ++i) {
        const auto& Ring = TubeA.Rings[i];
        UE_LOG(LogTemp, Log, TEXT("TubeA Ring %d: Center: %s, Radius: %f"), i, *Ring.Center.ToString(), Ring.Radius);
    }

    UE_LOG(LogTemp, Log, TEXT("TubeB has %d rings"), TubeB.Rings.Num());
    for (int32 i = 0; i < TubeB.Rings.Num(); ++i) {
        const auto& Ring = TubeB.Rings[i];
        UE_LOG(LogTemp, Log, TEXT("TubeB Ring %d: Center: %s, Radius: %f"), i, *Ring.Center.ToString(), Ring.Radius);
    }

    // Call the function to test
    TArray<int32> PlaneIndicesA, PlaneIndicesB, LineIndicesA, LineIndicesB;
    DoobGeometryUtils::FilterPlanesAndLines(TubeA, TubeB, PlaneIndicesA, PlaneIndicesB, LineIndicesA, LineIndicesB);

    // Verify results
    TestEqual("PlaneIndicesA size matches expected", PlaneIndicesA.Num(), 2);
    TestEqual("PlaneIndicesB size matches expected", PlaneIndicesB.Num(), 2);
    TestEqual("LineIndicesA size matches expected", LineIndicesA.Num(), 2);
    TestEqual("LineIndicesB size matches expected", LineIndicesB.Num(), 2);

    return true;
}

bool FCalculateIntersectionWithPlaneTest::RunTest(const FString& Parameters) {
    // Test Case 1: Perpendicular Ray and Plane
    {
        DoobGeometryUtils::FPlaneEquation Plane(FVector(0, 0, 1), -10); // Plane Z = 10
        FVector RayOrigin(0, 0, 0); // Ray starts at origin
        FVector RayDirection(0, 0, 1); // Ray points upward along Z-axis
        FVector IntersectionPoint;

        bool bIntersected = DoobGeometryUtils::CalculateIntersectionWithPlane(Plane, RayOrigin, RayDirection, IntersectionPoint);

        TestEqual("Perpendicular ray intersects plane", bIntersected, true);
        TestEqual("Intersection point is correct", IntersectionPoint, FVector(0, 0, 10));
    }

    // Test Case 2: Parallel Ray and Plane
    {
        DoobGeometryUtils::FPlaneEquation Plane(FVector(0, 0, 1), -10); // Plane Z = 10
        FVector RayOrigin(0, 0, 5); // Ray starts below the plane
        FVector RayDirection(1, 0, 0); // Ray points along X-axis (parallel to the plane)
        FVector IntersectionPoint;

        bool bIntersected = DoobGeometryUtils::CalculateIntersectionWithPlane(Plane, RayOrigin, RayDirection, IntersectionPoint);

        TestEqual("Parallel ray does not intersect plane", bIntersected, false);
    }

    // Test Case 3: Ray Originating on the Plane
    {
        DoobGeometryUtils::FPlaneEquation Plane(FVector(0, 0, 1), -10); // Plane Z = 10
        FVector RayOrigin(0, 0, 10); // Ray starts exactly on the plane
        FVector RayDirection(0, 1, 0); // Ray points along Y-axis
        FVector IntersectionPoint;

        // Log initial parameters
        UE_LOG(LogTemp, Log, TEXT("Testing Ray Originating on Plane"));
        UE_LOG(LogTemp, Log, TEXT("Plane: Normal=(%s), D=%f"), *Plane.Normal.ToString(), Plane.D);
        UE_LOG(LogTemp, Log, TEXT("Ray Origin: %s, Direction: %s"), *RayOrigin.ToString(), *RayDirection.ToString());


        bool bIntersected = DoobGeometryUtils::CalculateIntersectionWithPlane(Plane, RayOrigin, RayDirection, IntersectionPoint);

        // Log results
        UE_LOG(LogTemp, Log, TEXT("Intersection Result: %s"), bIntersected ? TEXT("True") : TEXT("False"));
        if (bIntersected) {
            UE_LOG(LogTemp, Log, TEXT("Intersection Point: %s"), *IntersectionPoint.ToString());
        }

        TestEqual("Ray originating on plane intersects plane", bIntersected, true);
        TestEqual("Intersection point is at origin", IntersectionPoint, RayOrigin);
    }

    // Test Case 4: Skew Ray
    {
        DoobGeometryUtils::FPlaneEquation Plane(FVector(0, 0, 1), -10); // Plane Z = 10
        FVector RayOrigin(0, 0, 0); // Ray starts at origin
        FVector RayDirection(1, 1, 1); // Ray points diagonally
        FVector IntersectionPoint;

        bool bIntersected = DoobGeometryUtils::CalculateIntersectionWithPlane(Plane, RayOrigin, RayDirection, IntersectionPoint);

        TestEqual("Skew ray intersects plane", bIntersected, true);
        TestEqual("Intersection point is correct", IntersectionPoint, FVector(10, 10, 10));
    }

    return true;
}