// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralTubularJointNode.h"

// Sets default values
AProceduralTubularJointNode::AProceduralTubularJointNode()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	// create the procedural mesh component
	TubularJointNodeMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("ProceduralMesh"));
	RootComponent = TubularJointNodeMesh;

	// intialize properties
	MainTubeSegments = 32;
	MainTubeRadius = 50.0f;
	bIsMainTubeClosed = true;

	MainTubeProfile = DoobProfileUtils::GenerateEggShapedCylinderProfile(MainTubeSegments, 0.0f, 300.0f, 0.0f, 0.35f);

	// Generate a default circular profile
	//MainTubeProfile = DoobProfileUtils::GenerateRegularCylinderProfile(MainTubeSegments, MainTubeRadius);
	MainTubeTransform = FTransform::Identity;
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
	TArray<TArray<FVector>> Rings;

	TArray<FVector> TubeVertices;

	TArray<int32> TubeTriangles;
	// initialize vertex index
	int32 BaseIndex = 0;


	DoobGeometryUtils::ConstructTubeFromProfile(MainTubeProfile, StartPosition, EndPosition, FVector(0, 0, 1), FVector(0, 1, 0), MainTubeProfile.Points.Num(), 10, 700.0f, Rings);

	DoobGeometryUtils::ConnectRingArray(Rings, TubeVertices, TubeTriangles, BaseIndex);

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