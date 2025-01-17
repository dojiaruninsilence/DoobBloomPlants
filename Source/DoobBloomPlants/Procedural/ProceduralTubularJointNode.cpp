// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralTubularJointNode.h"
#include "DrawDebugHelpers.h"

// Sets default values
AProceduralTubularJointNode::AProceduralTubularJointNode()
{
	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	// create the procedural mesh component
	TubularJointNodeMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("ProceduralMesh"));
	RootComponent = TubularJointNodeMesh;

	// intialize properties
	MainTubeSegments = 2;
	LateralTubeSegments = 2;
	MainTubeRadius = 50.0f;
	LateralTubeRadius = 25.0f;
	bIsMainTubeClosed = true;
	bIsLateralTubeClosed = true;

	//MainTubeProfile = DoobProfileUtils::GenerateEggShapedCylinderProfile(MainTubeSegments, 0.0f, 300.0f, 0.0f, 0.35f);

	// Generate a default circular profile
	MainTubeProfile = DoobProfileUtils::GenerateLinearProfile(MainTubeSegments, 200.0f, 200.0f);
	MainTubeTransform = FTransform::Identity;

	LateralTubeProfile = DoobProfileUtils::GenerateLinearProfile(LateralTubeSegments, 100.0f, 100.0f);
	LateralTubeTransform = FTransform::Identity;
}

// Called when the game starts or when spawned
void AProceduralTubularJointNode::BeginPlay()
{
	Super::BeginPlay();
	BuildMainTube();
}

// Called every frame
void AProceduralTubularJointNode::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

void AProceduralTubularJointNode::BuildMainTube() {
	DoobGeometryUtils::FTwoTubeIntersectionData TubeIntersection;

	//DoobGeometryUtils::FTubeData TubeData;
	//DoobGeometryUtils::FTubeData LateralTubeData;

	LateralStartPosition = StartPosition + 350.0f * FVector(0, 0, 1);

	/*TArray<FVector> TubeVertices;
	TArray<FVector> LateralTubeVertices;*/

	//DoobGeometryUtils::FIntersectionRingData TestIntersection;

	//TArray<int32> TubeTriangles;
	//TArray<int32> LateralTubeTriangles;
	//// initialize vertex index
	//int32 BaseIndex = 0;
	//int32 LateralBaseIndex = 0;

	// changed for now 
	/*DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, FVector(0, 0, 1), FVector(0, 1, 0), MainTubeProfile.Points.Num(), 10, 700.0f, TubeData, true);
	DoobGeometryUtils::ConstructTubeFromProfile(LateralTubeProfile, LateralStartPosition, FVector(0, 1, 0.25f), FVector(0, 0, 1), LateralTubeProfile.Points.Num(), 10, 400.0f, LateralTubeData, true);*/

	DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, FVector(0, 0, 1), FVector(0, 1, 0), MainTubeProfile.Points.Num(), 10, 700.0f, TubeIntersection.MainTube, true);
	DoobGeometryUtils::ConstructTubeFromProfile(LateralTubeProfile, LateralStartPosition, FVector(0, 1, 0.25f), FVector(0, 0, 1), LateralTubeProfile.Points.Num(), 10, 400.0f, TubeIntersection.LateralTube, true);

	// Get the World reference
	UWorld* World = GetWorld();

	//DoobGeometryUtils::GenerateIntersectionRing(TubeData, LateralTubeData, TestIntersection);
	//DoobGeometryUtils::GenerateHalfIntersectionRing(LateralTubeData, TubeData, TestIntersection.LateralTubeVertices);

	DoobGeometryUtils::GenerateIntersectionRing(TubeIntersection.MainTube, TubeIntersection.LateralTube, TubeIntersection.IntersectionRing);

	/*TArray<FVector> TestRemoveVerts;

	for (int32 MainRingIndex = 0; MainRingIndex < TubeIntersection.MainTube.Rings.Num() - 1; ++MainRingIndex) {
		DoobGeometryUtils::FRingData CurrentMainRing = TubeIntersection.MainTube.Rings[MainRingIndex];
		DoobGeometryUtils::FRingData NextMainRing = TubeIntersection.MainTube.Rings[MainRingIndex + 1];

		for (int32 LateralRingIndex = 0; LateralRingIndex < TubeIntersection.LateralTube.Rings.Num(); ++LateralRingIndex) {
			DoobGeometryUtils::FRingData CurrentLateralRing = TubeIntersection.LateralTube.Rings[LateralRingIndex];

			for (int32 VertexIndex = 0; VertexIndex < CurrentLateralRing.Vertices.Num(); ++VertexIndex) {
				FVector CurrentVertex = CurrentLateralRing.Vertices[VertexIndex];
				if (DoobGeometryUtils::IsPointInsideFrustum(CurrentMainRing, NextMainRing, CurrentVertex)) {
					TestRemoveVerts.Add(CurrentVertex);
				}
			}
		}
	}

	for (int32 i = 0; i < TestRemoveVerts.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TestRemoveVerts[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TestRemoveVerts[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);
	}*/

	// DoobGeometryUtils::RemoveInternalVertices(TubeIntersection.MainTube, TubeIntersection.LateralTube);

	/*for (int32 i = 0; i < TubeIntersection.LateralTube.Rings.Num(); ++i) {
		for (int32 j = 0; j < TubeIntersection.LateralTube.Rings[i].Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTube.Rings[i].Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTube.Rings[i].Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

			UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
		}
	}

	UE_LOG(LogTemp, Log, TEXT("Test intersection vertices: %d"), TubeIntersection.IntersectionRing.LateralTubeVertices.Num());*/

	/*for (int32 i = 0; i < TubeIntersection.IntersectionRing.CombinedVertices.Num() - 1; ++i) {
		// Convert local space to world space
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);


		UE_LOG(LogTemp, Log, TEXT("Intersection Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	DoobGeometryUtils::FindIntersectionRingCardinalPoints(TubeIntersection.IntersectionRing, StartPosition, TubeIntersection.MainTube.EndPosition);

	DoobGeometryUtils::GenerateAboveBelowIntersectionRings(TubeIntersection);

	DoobGeometryUtils::GenerateLateralTubeIntersectionRings(TubeIntersection);

	/*for (int32 i = 0; i < TubeIntersection.LateralTubeIntersectionRings.Rings.Num(); ++i) {
		for (int32 j = 0; j < TubeIntersection.LateralTubeIntersectionRings.Rings[i].Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeIntersectionRings.Rings[i].Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeIntersectionRings.Rings[i].Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

			UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
		}
	}*/

	/*for (int32 j = 0; j < TubeIntersection.LateralTubeFirstFullRing.Vertices.Num() - 1; ++j) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeFirstFullRing.Vertices[j]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeFirstFullRing.Vertices[j + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), j, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.LateralTube.Rings.Num(); ++i) {
		for (int32 j = 0; j < TubeIntersection.LateralTube.Rings[i].Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTube.Rings[i].Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTube.Rings[i].Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

			UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
		}
	}

	UE_LOG(LogTemp, Log, TEXT("Test intersection vertices: %d"), TubeIntersection.IntersectionRing.LateralTubeVertices.Num());*/

	/*int32 LowestIndex;
	FVector LowestVertex;
	DoobGeometryUtils::FindClosestVertex(StartPosition, TubeIntersection.IntersectionRing.CombinedVertices, LowestVertex, LowestIndex);
	int32 HighestIndex;
	FVector HighestVertex;
	DoobGeometryUtils::FindClosestVertex(TubeIntersection.MainTube.EndPosition, TubeIntersection.IntersectionRing.CombinedVertices, HighestVertex, HighestIndex);

	int32 MidIndex = (HighestIndex + LowestIndex) / 2;
	int32 MidIndexB = LowestIndex + MidIndex; */

	//UE_LOG(LogTemp, Log, TEXT("MidIndex: %d, LowestIndex: %d, HighestIndex: %d, number: %d"), MidIndex, LowestIndex, HighestIndex, TubeIntersection.IntersectionRing.CombinedVertices.Num());

	/*DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[0]), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f);
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[1]), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f);
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[2]), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f); // left
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[3]), 5.0f, 12, FColor::Purple, true, 0.0f, 0, 1.0f); // right

	DoobGeometryUtils::FRingData LowVertStartRing;
	DoobGeometryUtils::FRingData LowVertEndRing;
	DoobGeometryUtils::FRingData HighVertStartRing;
	DoobGeometryUtils::FRingData HighVertEndRing;*/

	//UE_LOG(LogTemp, Log, TEXT("---------------------------------------------Lowest Vertex: (X:%f/Y:%f/Z:%f)"), LowestVertex.X, LowestVertex.Y, LowestVertex.Z);

	/*DoobGeometryUtils::FindSegmentForPoint(TubeIntersection.IntersectionRing.CardinalVertices[2], TubeIntersection.MainTube, LowVertStartRing, LowVertEndRing);
	DoobGeometryUtils::FindSegmentForPoint(TubeIntersection.IntersectionRing.CardinalVertices[0], TubeIntersection.MainTube, HighVertStartRing, HighVertEndRing);

	UE_LOG(LogTemp, Log, TEXT("---------------------------------------------Lowest Vertex Start Ring Center: (X:%f/Y:%f/Z:%f)"), LowVertStartRing.Center.X, LowVertStartRing.Center.Y, LowVertStartRing.Center.Z);
	UE_LOG(LogTemp, Log, TEXT("---------------------------------------------Lowest Vertex End Ring Center: (X:%f/Y:%f/Z:%f)"), LowVertEndRing.Center.X, LowVertEndRing.Center.Y, LowVertEndRing.Center.Z);

	FVector LowestVertexCenter = DoobGeometryUtils::CalculateCenterLinePoint(TubeIntersection.IntersectionRing.CardinalVertices[2], LowVertStartRing, LowVertEndRing);
	FVector HighestVertexCenter = DoobGeometryUtils::CalculateCenterLinePoint(TubeIntersection.IntersectionRing.CardinalVertices[0], HighVertStartRing, HighVertEndRing);*/

	/*DrawDebugSphere(World, GetActorTransform().TransformPosition(LowestVertexCenter), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f);
	DrawDebugLine(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[2]), GetActorTransform().TransformPosition(LowestVertexCenter), FColor::Red, true, 0.0f, 0, 10.0f);*/

	//UE_LOG(LogTemp, Log, TEXT("---------------------------------------------Lowest Vertex Center: (X:%f/Y:%f/Z:%f)"), LowestVertexCenter.X, LowestVertexCenter.Y, LowestVertexCenter.Z);

	// --------------------------------------------------------------- Temp to find intersection of lateral tube ---------------------------------------------------------------//

	/*DoobGeometryUtils::FRingData TestLateralIntersectionLowRing;
	DoobGeometryUtils::FRingData TestLateralIntersectionHighRing;
	DoobGeometryUtils::FRingData TestLateralIntersectionRing;
	FVector TestHighestVertex;
	int32 TestHighestIndex = 0;

	TArray<FVector> TestNorthSouth = { TubeIntersection.IntersectionRing.CardinalVertices[0], TubeIntersection.IntersectionRing.CardinalVertices[2] };
	DoobGeometryUtils::FindClosestVertex(TubeIntersection.LateralTube.EndPosition, TestNorthSouth, TestHighestVertex, TestHighestIndex);
	if (TestHighestIndex == 0) {

	}*/

	// --------------------------------------------------------------- Temp to find intersection of lateral tube end -----------------------------------------------------------//

	//DoobGeometryUtils::FRingData LowestVertexRing;
	//DoobGeometryUtils::FRingData HighestVertexRing;

	// DoobGeometryUtils::GenerateRing(LowestVertexCenter, FVector(0, 0, 1), FVector(0, 1, 0), 200, 10, LowestVertexRing);
	// DoobGeometryUtils::GenerateRing(HighestVertexCenter, FVector(0, 0, 1), FVector(0, 1, 0), 200, 10, HighestVertexRing);

	//DoobGeometryUtils::GenerateRing(LowestVertexCenter, FVector(0, 0, 1), FVector(0, 1, 0), 200, 10, TubeIntersection.MTBelowIntersectionRing);
	//DoobGeometryUtils::GenerateRing(HighestVertexCenter, FVector(0, 0, 1), FVector(0, 1, 0), 200, 10, TubeIntersection.MTAboveIntersectionRing);

	/*for (int32 i = 0; i < TubeIntersection.MTBelowIntersectionRing.Vertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MTBelowIntersectionRing.Vertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MTBelowIntersectionRing.Vertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.MTAboveIntersectionRing.Vertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRing.Vertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRing.Vertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	//TArray<FVector> TestSquare;
	//DoobGeometryUtils::GenerateSquareAroundIntersection(TubeIntersection.MTBelowIntersectionRing, TubeIntersection.MTAboveIntersectionRing, TubeIntersection.IntersectionRing, MidIndexB, MidIndex, TestSquare);
	DoobGeometryUtils::GenerateSquareAroundIntersection(
		TubeIntersection.MTBelowIntersectionRing, 
		TubeIntersection.MTAboveIntersectionRing, 
		TubeIntersection.IntersectionRing, 
		TubeIntersection.IntersectionRing.CardinalIndices[3], 
		TubeIntersection.IntersectionRing.CardinalIndices[1], 
		TubeIntersection.IntersectionSquare.Corners
	);

	/*DrawDebugLine(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[0]), GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[3]), FColor::Red, true, 0.0f, 0, 10.0f);
	DrawDebugLine(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[1]), GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[2]), FColor::Red, true, 0.0f, 0, 10.0f);*/

	DoobGeometryUtils::RemoveVerticesByInterpolatedDirections(TubeIntersection, TubeIntersection.IntersectionSquare);

	/*for (int32 i = 0; i < TubeIntersection.MTAboveIntersectionRingPartial.Vertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRingPartial.Vertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRingPartial.Vertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Blue, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	//DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRingPartial.Vertices[0]), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f);

	for (int32 i = 0; i < TubeIntersection.MTBelowIntersectionRingPartial.Vertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MTBelowIntersectionRingPartial.Vertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MTBelowIntersectionRingPartial.Vertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Blue, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	DoobGeometryUtils::OrderSquareIntersectionConnections(TubeIntersection);

	/*for (int32 i = 0; i < TubeIntersection.IntersectionRing.TopRightVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.TopRightVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.TopRightVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Yellow, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionRing.BottomRightVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.BottomRightVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.BottomRightVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Orange, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionRing.TopLeftVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.TopLeftVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.TopLeftVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Orange, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionRing.BottomLeftVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.BottomLeftVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.BottomLeftVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Yellow, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionSquare.TopRightVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.TopRightVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.TopRightVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Yellow, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionSquare.BottomRightVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.BottomRightVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.BottomRightVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Orange, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionSquare.TopLeftVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.TopLeftVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.TopLeftVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Orange, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.IntersectionSquare.BottomLeftVertices.Num() - 1; ++i) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.BottomLeftVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.BottomLeftVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Yellow, true, 0.0f, 0, 10.0f);

		UE_LOG(LogTemp, Log, TEXT("Lowest Vertex Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	DoobGeometryUtils::ConnectTwoTubeIntersection(TubeIntersection);

	//DoobGeometryUtils::ConnectRingArray(TubeIntersection.MainTube.Rings, TubeVertices, TubeTriangles, BaseIndex);
	//DoobGeometryUtils::ConnectRingArray(LateralTubeData.Rings, TubeVertices, TubeTriangles, BaseIndex);

	TArray<FVector> Normals;
	//Normals.Init(FVector::UpVector, TubeVertices.Num());
	Normals.Init(FVector::UpVector, TubeIntersection.AllVertices.Num());

	//DoobMeshUtils::RemoveDegenerateTriangles(TubeVertices, TubeTriangles);

	TArray<FVector2d> UV0;
	for (const FVector& Vertex : TubeIntersection.AllVertices)
	{
		UV0.Add(FVector2D(Vertex.X, Vertex.Y));
	}

	TArray<FColor> VertexColors;
	VertexColors.Init(FColor::White, TubeIntersection.AllVertices.Num());

	TArray<FProcMeshTangent> Tangents;
	Tangents.Init(FProcMeshTangent(), TubeIntersection.AllVertices.Num());

	UE_LOG(LogTemp, Log, TEXT("Vertices: %d, Normals: %d, UV0: %d, VertexColors: %d, Tangents: %d"),
		TubeIntersection.AllVertices.Num(), Normals.Num(), UV0.Num(), VertexColors.Num(), Tangents.Num());

	/*// temp logs
	UE_LOG(LogTemp, Log, TEXT("AllVertices Count: %d"), TubeIntersection.AllVertices.Num());
	for (int32 i = 0; i < TubeIntersection.AllVertices.Num(); i++)
	{
		UE_LOG(LogTemp, Log, TEXT("Vertex[%d]: %s"), i, *TubeIntersection.AllVertices[i].ToString());
	}

	UE_LOG(LogTemp, Log, TEXT("Triangles Count: %d"), TubeIntersection.Triangles.Num());
	for (int32 i = 0; i < TubeIntersection.Triangles.Num(); i += 3)
	{
		UE_LOG(LogTemp, Log, TEXT("Triangle[%d]: %d, %d, %d"), i / 3, TubeIntersection.Triangles[i], TubeIntersection.Triangles[i + 1], TubeIntersection.Triangles[i + 2]);
	}

	UE_LOG(LogTemp, Log, TEXT("Normals Count: %d"), Normals.Num());
	for (int32 i = 0; i < Normals.Num(); i++)
	{
		UE_LOG(LogTemp, Log, TEXT("Normal[%d]: %s"), i, *Normals[i].ToString());
	}

	UE_LOG(LogTemp, Log, TEXT("UV0 Count: %d"), UV0.Num());
	for (int32 i = 0; i < UV0.Num(); i++)
	{
		UE_LOG(LogTemp, Log, TEXT("UV0[%d]: (%f, %f)"), i, UV0[i].X, UV0[i].Y);
	}

	UE_LOG(LogTemp, Log, TEXT("VertexColors Count: %d"), VertexColors.Num());
	for (int32 i = 0; i < VertexColors.Num(); i++)
	{
		UE_LOG(LogTemp, Log, TEXT("VertexColor[%d]: (%d, %d, %d, %d)"),
			i,
			VertexColors[i].R,
			VertexColors[i].G,
			VertexColors[i].B,
			VertexColors[i].A);
	}

	UE_LOG(LogTemp, Log, TEXT("Tangents Count: %d"), Tangents.Num());
	for (int32 i = 0; i < Tangents.Num(); i++)
	{
		UE_LOG(LogTemp, Log, TEXT("Tangent[%d]: (%f, %f, %f)"), i, Tangents[i].TangentX.X, Tangents[i].TangentX.Y, Tangents[i].TangentX.Z);
	}*/

	// Create the node mesh section
	TubularJointNodeMesh->CreateMeshSection(0, TubeIntersection.AllVertices, TubeIntersection.Triangles, Normals, UV0, VertexColors, Tangents, true);

	/*UE_LOG(LogTemp, Log, TEXT("Creating Mesh Section..."));
	UE_LOG(LogTemp, Log, TEXT("Vertices: %d, Triangles: %d, Normals: %d, UV0: %d, VertexColors: %d, Tangents: %d"),
		TubeIntersection.AllVertices.Num(),
		TubeIntersection.Triangles.Num(),
		Normals.Num(),
		UV0.Num(),
		VertexColors.Num(),
		Tangents.Num());

	if (TubularJointNodeMesh)
	{
		UE_LOG(LogTemp, Log, TEXT("Procedural Mesh Component Found."));
		UE_LOG(LogTemp, Log, TEXT("Mesh Visibility: %s"), TubularJointNodeMesh->IsVisible() ? TEXT("True") : TEXT("False"));
		FBoxSphereBounds Bounds = TubularJointNodeMesh->CalcBounds(FTransform());
		UE_LOG(LogTemp, Log, TEXT("Mesh Bounds: Origin=%s, BoxExtent=%s, SphereRadius=%f"),
			*Bounds.Origin.ToString(),
			*Bounds.BoxExtent.ToString(),
			Bounds.SphereRadius);
	}
	else
	{
		UE_LOG(LogTemp, Error, TEXT("Procedural Mesh Component is NULL."));
	}*/
}

//void AProceduralTubularJointNode::CreateProceduralMesh() {
//
//}