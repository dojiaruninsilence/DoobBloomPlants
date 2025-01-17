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

	LateralStartPosition = StartPosition + 350.0f * FVector(0, 0, 1);

	DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, FVector(0, 0, 1), FVector(0, 1, 0), MainTubeProfile.Points.Num(), 10, 700.0f, TubeIntersection.MainTube, true);
	DoobGeometryUtils::ConstructTubeFromProfile(LateralTubeProfile, LateralStartPosition, FVector(0, 1, 0.25f), FVector(0, 0, 1), LateralTubeProfile.Points.Num(), 10, 400.0f, TubeIntersection.LateralTube, true);

	// Get the World reference
	UWorld* World = GetWorld();

	DoobGeometryUtils::GenerateIntersectionRing(TubeIntersection.MainTube, TubeIntersection.LateralTube, TubeIntersection.IntersectionRing);

	DoobGeometryUtils::FindIntersectionRingCardinalPoints(TubeIntersection.IntersectionRing, StartPosition, TubeIntersection.MainTube.EndPosition);

	DoobGeometryUtils::GenerateAboveBelowIntersectionRings(TubeIntersection);

	DoobGeometryUtils::GenerateLateralTubeIntersectionRings(TubeIntersection);

	DoobGeometryUtils::GenerateSquareAroundIntersection(
		TubeIntersection.MTBelowIntersectionRing, 
		TubeIntersection.MTAboveIntersectionRing, 
		TubeIntersection.IntersectionRing, 
		TubeIntersection.IntersectionRing.CardinalIndices[3], 
		TubeIntersection.IntersectionRing.CardinalIndices[1], 
		TubeIntersection.IntersectionSquare.Corners
	);

	DoobGeometryUtils::RemoveVerticesByInterpolatedDirections(TubeIntersection, TubeIntersection.IntersectionSquare);

	DoobGeometryUtils::OrderSquareIntersectionConnections(TubeIntersection);

	DoobGeometryUtils::ConnectTwoTubeIntersection(TubeIntersection);

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