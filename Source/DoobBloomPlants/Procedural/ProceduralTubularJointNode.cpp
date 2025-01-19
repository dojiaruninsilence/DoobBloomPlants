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
	MainTubeSegments = 20;
	LateralTubeSegments = 20;
	MainTubeRadius = 50.0f;
	LateralTubeRadius = 25.0f;
	bIsMainTubeClosed = true;
	bIsLateralTubeClosed = true;

	//MainTubeProfile = DoobProfileUtils::GenerateEggShapedCylinderProfile(MainTubeSegments, 0.0f, 300.0f, 100.0f, 0.35f);

	// Generate a default circular profile
	MainTubeProfile = DoobProfileUtils::GenerateLinearProfile(MainTubeSegments, 200.0f, 100.0f);
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

	LateralStartPosition = StartPosition + 350.0f * FVector(0, 0, 1);

	DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, FVector(0, 0, 1), FVector(0, 1, 0), MainTubeProfile.Points.Num(), 10, 700.0f, TubeIntersection.MainTube, true);
	DoobGeometryUtils::ConstructTubeFromProfile(LateralTubeProfile, LateralStartPosition, FVector(0, 1, 0.25f), FVector(0, 0, 1), LateralTubeProfile.Points.Num(), 10, 400.0f, TubeIntersection.LateralTube, true);

	// Get the World reference
	UWorld* World = GetWorld();

	for (int32 i = 0; i < TubeIntersection.MainTube.Rings.Num(); ++i) {
		for (int32 j = 0; j < TubeIntersection.MainTube.Rings[i].Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MainTube.Rings[i].Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MainTube.Rings[i].Vertices[j + 1]);
			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Blue, true, 0.0f, 0, 1.0f);
			//UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
		}
	}

	DoobGeometryUtils::GenerateIntersectionRing(TubeIntersection.MainTube, TubeIntersection.LateralTube, TubeIntersection.IntersectionRing);

	/*for (int32 j = 0; j < TubeIntersection.IntersectionRing.CombinedVertices.Num() - 1; ++j) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[j]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CombinedVertices[j + 1]);
		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Cyan, true, 0.0f, 0, 10.0f);
		UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), j, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	DoobGeometryUtils::FindIntersectionRingCardinalPoints(TubeIntersection.IntersectionRing, StartPosition, TubeIntersection.MainTube.EndPosition);

	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[0]), 5.0f, 12, FColor::Yellow, true, 0.0f, 0, 1.0f);
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[1]), 5.0f, 12, FColor::Red, true, 0.0f, 0, 1.0f);
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[2]), 5.0f, 12, FColor::Green, true, 0.0f, 0, 1.0f); // left
	DrawDebugSphere(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionRing.CardinalVertices[3]), 5.0f, 12, FColor::Purple, true, 0.0f, 0, 1.0f); // right

	DoobGeometryUtils::GenerateAboveBelowIntersectionRings(TubeIntersection);

	DoobGeometryUtils::GenerateLateralTubeIntersectionRings(TubeIntersection);

	TubeIntersection.LateralTubeRemovedVertices = TubeIntersection.LateralTube;
	DoobGeometryUtils::RemoveInternalVertices(TubeIntersection.MainTube, TubeIntersection.LateralTubeRemovedVertices);

	//for (int32 i = 0; i < TubeIntersection.LateralTubeRemovedVertices.Rings.Num(); ++i) {
	//	for (int32 j = 0; j < TubeIntersection.LateralTubeRemovedVertices.Rings[i].Vertices.Num() - 1; ++j) {
	//		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeRemovedVertices.Rings[i].Vertices[j]);
	//		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeRemovedVertices.Rings[i].Vertices[j + 1]);
	//		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 1.0f);
	//		//UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	//	}
	//}

	for (int32 i = 0; i < TubeIntersection.LateralTubeIntersectionRings.Rings.Num(); ++i) {
		for (int32 j = 0; j < TubeIntersection.LateralTubeIntersectionRings.Rings[i].Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeIntersectionRings.Rings[i].Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.LateralTubeIntersectionRings.Rings[i].Vertices[j + 1]);
			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Green, true, 0.0f, 0, 1.0f);
			//UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
		}
	}

	/*for (int32 j = 0; j < TubeIntersection.MTAboveIntersectionRing.Vertices.Num() - 1; ++j) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRing.Vertices[j]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MTAboveIntersectionRing.Vertices[j + 1]);
		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Cyan, true, 0.0f, 0, 10.0f);
		UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), j, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	for (int32 j = 0; j < TubeIntersection.MTBelowIntersectionRing.Vertices.Num() - 1; ++j) {
		FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MTBelowIntersectionRing.Vertices[j]);
		FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MTBelowIntersectionRing.Vertices[j + 1]);
		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Cyan, true, 0.0f, 0, 10.0f);
		UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), j, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}*/

	DoobGeometryUtils::GenerateSquareAroundIntersection(
		TubeIntersection.MTBelowIntersectionRing, 
		TubeIntersection.MTAboveIntersectionRing, 
		TubeIntersection.IntersectionRing, 
		TubeIntersection.IntersectionRing.CardinalIndices[3], 
		TubeIntersection.IntersectionRing.CardinalIndices[1], 
		TubeIntersection.IntersectionSquare.Corners
	);

	DrawDebugLine(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[0]), GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[3]), FColor::Cyan, true, 0.0f, 0, 10.0f);
	DrawDebugLine(World, GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[1]), GetActorTransform().TransformPosition(TubeIntersection.IntersectionSquare.Corners[2]), FColor::Cyan, true, 0.0f, 0, 10.0f);

	DoobGeometryUtils::RemoveVerticesByInterpolatedDirections(TubeIntersection, TubeIntersection.IntersectionSquare);

	DoobGeometryUtils::OrderSquareIntersectionConnections(TubeIntersection);

	DoobGeometryUtils::ConnectTwoTubeIntersection(TubeIntersection);

	for (int32 i = 0; i < TubeIntersection.MainTubePartialRings.Num(); ++i) {
		for (int32 j = 0; j < TubeIntersection.MainTubePartialRings[i].Vertices.Num() - 1; ++j) {
			FVector CurrentVertex = GetActorTransform().TransformPosition(TubeIntersection.MainTubePartialRings[i].Vertices[j]);
			FVector NextVertex = GetActorTransform().TransformPosition(TubeIntersection.MainTubePartialRings[i].Vertices[j + 1]);
			DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 1.0f);
			//UE_LOG(LogTemp, Log, TEXT("Intersection Remaining Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
		}
	}

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