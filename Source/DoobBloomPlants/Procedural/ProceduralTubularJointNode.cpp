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
	DoobGeometryUtils::FTubeData TubeData;
	DoobGeometryUtils::FTubeData LateralTubeData;

	LateralStartPosition = StartPosition + 350.0f * FVector(0, 0, 1);

	TArray<FVector> TubeVertices;
	TArray<FVector> LateralTubeVertices;

	DoobGeometryUtils::FIntersectionRingData TestIntersection;

	TArray<int32> TubeTriangles;
	TArray<int32> LateralTubeTriangles;
	// initialize vertex index
	int32 BaseIndex = 0;
	int32 LateralBaseIndex = 0;

	DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, FVector(0, 0, 1), FVector(0, 1, 0), MainTubeProfile.Points.Num(), 10, 700.0f, TubeData, true);
	DoobGeometryUtils::ConstructTubeFromProfile(LateralTubeProfile, LateralStartPosition, FVector(0, 1, 0.25f), FVector(0, 0, 1), LateralTubeProfile.Points.Num(), 10, 400.0f, LateralTubeData, true);

	// Get the World reference
	UWorld* World = GetWorld();

	DoobGeometryUtils::GenerateIntersectionRing(TubeData, LateralTubeData, TestIntersection);
	//DoobGeometryUtils::GenerateHalfIntersectionRing(LateralTubeData, TubeData, TestIntersection.LateralTubeVertices);

	UE_LOG(LogTemp, Log, TEXT("Test intersection vertices: %d"), TestIntersection.LateralTubeVertices.Num());

	for (int32 i = 0; i < TestIntersection.CombinedVertices.Num() - 1; ++i) {
		// Convert local space to world space
		FVector CurrentVertex = GetActorTransform().TransformPosition(TestIntersection.CombinedVertices[i]);
		FVector NextVertex = GetActorTransform().TransformPosition(TestIntersection.CombinedVertices[i + 1]);

		DrawDebugLine(World, CurrentVertex, NextVertex, FColor::Red, true, 0.0f, 0, 10.0f);

		
		UE_LOG(LogTemp, Log, TEXT("Intersection Ring Debug iteration: %d, Current Vertex: (X:%f/Y:%f/Z:%f), Next Vertex: (X:%f/Y:%f/Z:%f)"), i, CurrentVertex.X, CurrentVertex.Y, CurrentVertex.Z, NextVertex.X, NextVertex.Y, NextVertex.Z);
	}

	DoobGeometryUtils::ConnectRingArray(TubeData.Rings, TubeVertices, TubeTriangles, BaseIndex);
	DoobGeometryUtils::ConnectRingArray(LateralTubeData.Rings, TubeVertices, TubeTriangles, BaseIndex);

	TArray<FVector> Normals;
	Normals.Init(FVector::UpVector, TubeVertices.Num());

	//DoobMeshUtils::RemoveDegenerateTriangles(TubeVertices, TubeTriangles);

	TArray<FVector2d> UV0;
	for (const FVector& Vertex : TubeVertices)
	{
		UV0.Add(FVector2D(Vertex.X, Vertex.Y));
	}

	TArray<FColor> VertexColors;
	VertexColors.Init(FColor::White, TubeVertices.Num());

	TArray<FProcMeshTangent> Tangents;
	Tangents.Init(FProcMeshTangent(), TubeVertices.Num());

	UE_LOG(LogTemp, Log, TEXT("Vertices: %d, Normals: %d, UV0: %d, VertexColors: %d, Tangents: %d"),
		TubeVertices.Num(), Normals.Num(), UV0.Num(), VertexColors.Num(), Tangents.Num());


	// Create the node mesh section
	TubularJointNodeMesh->CreateMeshSection(0, TubeVertices, TubeTriangles, Normals, UV0, VertexColors, Tangents, true);
}

//void AProceduralTubularJointNode::CreateProceduralMesh() {
//
//}