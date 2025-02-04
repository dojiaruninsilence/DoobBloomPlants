// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralTubularJointNode.h"
#include "DrawDebugHelpers.h"

#include "DoobMathUtils.h"

// Sets default values
AProceduralTubularJointNode::AProceduralTubularJointNode()
{
	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	// create the procedural mesh component
	TubularJointNodeMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("ProceduralMesh"));
	RootComponent = TubularJointNodeMesh;

	// intialize properties
	MainTubeLength = 700.0f;
	LateralTubeLength = 400.0f;
	MainTubeSegments = 20;
	LateralTubeSegments = 20;
	MainTubeDirection = FVector(0, 0, 1);
	MainTubeStartRadius = 0.0f;
	MainTubeMidRadius = 300.0f;
	MainTubeEndRadius = 100.0f;
	LateralTubeRadius = 25.0f;
	MainTubeNumSides = 10;
	LateralTubeNumSides = 10;
	bIsMainTubeClosed = true;
	bIsLateralTubeClosed = true;

	MainTubeProfile = DoobProfileUtils::GenerateEggShapedCylinderProfile(MainTubeSegments, MainTubeStartRadius, MainTubeMidRadius, MainTubeEndRadius, 0.35f);

	LateralTubeProfile = DoobProfileUtils::GenerateEggShapedCylinderProfile(MainTubeSegments, 0.0f, 200.0f, 50.0f, 0.75f);

	// Generate a default circular profile
	//MainTubeProfile = DoobProfileUtils::GenerateLinearProfile(MainTubeSegments, 200.0f, 100.0f);
	MainTubeTransform = FTransform::Identity;

	//LateralTubeProfile = DoobProfileUtils::GenerateLinearProfile(LateralTubeSegments, 100.0f, 50.0f);
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

	DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, MainTubeDirection, FVector(0, 1, 0), MainTubeProfile.Points.Num(), MainTubeNumSides, MainTubeLength, TubeIntersection.MainTube, true);

	LateralStartPosition = StartPosition + (MainTubeLength / 2) * MainTubeDirection;
	// ------------------------------------ either this doesnt work or directions need some work. need to check
	//LateralTubeDirection = DoobMathUtils::GetRandomDirectionBetween(DoobMathUtils::GenerateRandomPerpendicularVector(MainTubeDirection), MainTubeDirection, 0.2f, 0.9f);
	LateralTubeDirection = FVector(0, 1, 0.75f);

	DoobGeometryUtils::ConstructTubeFromProfile(LateralTubeProfile, LateralStartPosition, LateralTubeDirection, FVector(0, 0, 1), LateralTubeProfile.Points.Num(), LateralTubeNumSides, LateralTubeLength, TubeIntersection.LateralTube, true);

	// Get the World reference
	UWorld* World = GetWorld();

	DoobGeometryUtils::GenerateIntersectionRing(TubeIntersection.MainTube, TubeIntersection.LateralTube, TubeIntersection.IntersectionRing);

	for (int32 j = 0; j < TubeIntersection.IntersectionRing.CombinedVertices.Num() - 1; ++j) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[j]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[j + 1]);
		DrawDebugSphere(World, CurrentVertex, 5.0f, 12, FColor::Black, true, 0.0f, 0, 1.0f);
		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Cyan, true, 0.0f, 0, 1.0f);
		UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), j, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 i = 0; i < TubeIntersection.LateralTube.Rings.Num(); ++i) {
		DoobGeometryUtils::FRingData CurrentRing = TubeIntersection.LateralTube.Rings[i];
		for (int32 j = 0; j < CurrentRing.Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Green, true, 0.0f, 0, 1.0f);
		}
	}

	for (int32 i = 0; i < TubeIntersection.MainTube.Rings.Num(); ++i) {
		DoobGeometryUtils::FRingData CurrentRing = TubeIntersection.MainTube.Rings[i];
		DrawDebugSphere(World, GetActorTransform().TransformPosition(CurrentRing.Vertices[0]), 5.0f, 12, FColor::Black, true, 0.0f, 0, 1.0f);
		DrawDebugSphere(World, GetActorTransform().TransformPosition(CurrentRing.Vertices[1]), 5.0f, 12, FColor::White, true, 0.0f, 0, 1.0f);
		for (int32 j = 0; j < CurrentRing.Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Green, true, 0.0f, 0, 1.0f);
		}
	}

	DoobGeometryUtils::FindIntersectionRingCardinalPoints(TubeIntersection.IntersectionRing, StartPosition, TubeIntersection.MainTube.EndPosition);

	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[0]), 5.0f, 12, FColor::Black, true, 0.0f, 0, 1.0f);

	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[0]), 5.0f, 12, FColor::Yellow, true, 0.0f, 0, 1.0f);
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[1]), 5.0f, 12, FColor::White, true, 0.0f, 0, 1.0f);
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[2]), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f); // left
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[3]), 5.0f, 12, FColor::Purple, true, 0.0f, 0, 1.0f); // right

	DoobGeometryUtils::GenerateAboveBelowIntersectionRings(TubeIntersection);

	TubeIntersection.LateralTubeRemovedVertices = TubeIntersection.LateralTube;
	DoobGeometryUtils::RemoveInternalVertices(TubeIntersection.MainTube, TubeIntersection.LateralTubeRemovedVertices);

	DoobGeometryUtils::GenerateLateralTubeIntersectionRings(TubeIntersection);

	DoobGeometryUtils::GenerateSquareAroundIntersection(
		TubeIntersection.MTBelowIntersectionRing, 
		TubeIntersection.MTAboveIntersectionRing, 
		TubeIntersection.IntersectionRing, 
		TubeIntersection.IntersectionRing.CardinalIndices[3], 
		TubeIntersection.IntersectionRing.CardinalIndices[1], 
		TubeIntersection.IntersectionSquare.Corners
	);

	DoobGeometryUtils::OrderSquareIntersectionConnections(TubeIntersection);

	//DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.TestVerts[0]), 5.0f, 12, FColor::Purple, true, 0.0f, 0, 1.0f);

	/*for (int32 j = 0; j < TubeIntersection.TestVerts.Num() - 1; ++j) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.TestVerts[j]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.TestVerts[j + 1]);
		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Purple, true, 0.0f, 0, 1.0f);
		UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), j, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	for (int32 i = 0; i < TubeIntersection.IntersectionSquare.TopRightPartialRings.Num(); ++i) {
		TArray<FVector> CurrentRing = TubeIntersection.IntersectionSquare.TopRightPartialRings[i];
		for (int32 j = 0; j < CurrentRing.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(CurrentRing[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(CurrentRing[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 1.0f);
		}
	}

	DoobGeometryUtils::RemoveVerticesByInterpolatedDirections(TubeIntersection, TubeIntersection.IntersectionSquare);

	DoobGeometryUtils::ConnectTwoTubeIntersection(TubeIntersection);

	/*for (int32 i = 0; i < TubeIntersection.MainTubePartialRings.Num(); ++i) {
		DoobGeometryUtils::FRingData CurrentRing = TubeIntersection.MainTubePartialRings[i];
		for (int32 j = 0; j < CurrentRing.Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 1.0f);
		}
	}*/

	/*for (int32 i = 0; i < TubeIntersection.LateralTubeIntersectionRings.Rings.Num(); ++i) {
		DoobGeometryUtils::FRingData CurrentRing = TubeIntersection.LateralTubeIntersectionRings.Rings[i];
		for (int32 j = 0; j < CurrentRing.Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(CurrentRing.Vertices[j + 1]);

			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 1.0f);
		}
	}*/

	DoobMeshUtils::RemoveDegenerateTriangles(TubeIntersection.AllVertices, TubeIntersection.Triangles);

	DoobMeshUtils::RemoveInwardFacingTriangles(TubeIntersection.AllVertices, TubeIntersection.Triangles);

	TArray<FVector> Normals;
	Normals.Init(FVector::UpVector, TubeIntersection.AllVertices.Num());

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

	// Create the node mesh section
	TubularJointNodeMesh->CreateMeshSection(0, TubeIntersection.AllVertices, TubeIntersection.Triangles, Normals, UV0, VertexColors, Tangents, true);
}