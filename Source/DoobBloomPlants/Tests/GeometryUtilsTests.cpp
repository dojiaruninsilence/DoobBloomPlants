#include "CoreMinimal.h"
#include "DoobGeometryUtils.h"
#include "Misc/AutomationTest.h"

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFilterPlanesAndLinesTest, "Doob.Geometry.FilterPlanesAndLines", EAutomationTestFlags::EditorContext | EAutomationTestFlags::EngineFilter)

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
