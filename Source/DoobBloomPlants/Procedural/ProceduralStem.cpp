// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralStem.h"

// Sets default values
AProceduralStem::AProceduralStem()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	//PrimaryActorTick.bCanEverTick = true;

	// Create the Procedural Mesh Component
	StemMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("StemMesh"));
	RootComponent = StemMesh;

	// Enable Collision
	StemMesh->bUseAsyncCooking = true;

}

// Called when the game starts or when spawned
void AProceduralStem::BeginPlay()
{
	Super::BeginPlay();
	GenerateStem();
	
}

// Called every frame
void AProceduralStem::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

// Generate the full stem
void AProceduralStem::GenerateStem()
{
	TArray<FVector> Vertices;
	TArray<int32> Triangles;

	FVector CurrentPosition = FVector::ZeroVector;
	FVector CurrentDirection = FVector(0, 0, 1); // Initial growth direction up

	for (int32 i = 0; i < NumSegments; ++i)
	{
		// Determine the next segment's end pos
		FVector NextPosition = CurrentPosition + (CurrentDirection * SegmentLength);

		// Generate the cylinder for this segment
		GenerateCylinder(CurrentPosition, NextPosition, SegmentRadius, Vertices, Triangles);
		
		// Update the current position and apply a random rotation
		CurrentPosition = NextPosition;
		FRotator RandomDirection = FRotator(FMath::RandRange(-30.0f, 30.0f), FMath::RandRange(-30.0f, 30.0f), 0);
		CurrentDirection = RandomDirection.RotateVector(CurrentDirection);
		CurrentDirection.Normalize();
	}

	TArray<FVector> Normals;
	Normals.Init(FVector::UpVector, Vertices.Num());

	TArray<FVector2d> UV0;
	for (const FVector& Vertex : Vertices)
	{
		UV0.Add(FVector2D(Vertex.X, Vertex.Y));
	}

	TArray<FColor> VertexColors;
	VertexColors.Init(FColor::White, Vertices.Num());

	TArray<FProcMeshTangent> Tangents;
	Tangents.Init(FProcMeshTangent(), Vertices.Num());

	// create the mesh section
	StemMesh->CreateMeshSection(0, Vertices, Triangles, Normals, UV0, VertexColors,Tangents, true);
}

// Generate a single cylinder
void AProceduralStem::GenerateCylinder(FVector Start, FVector End, float Radius, TArray<FVector>& Vertices, TArray<int32>& Triangles)
{
	const int NumSides = 12; // number of sides for the cylinder
	const float AngleStep = 360.0f / NumSides;

	FVector Up = (End - Start).GetSafeNormal();
	FVector Right = FVector::CrossProduct(Up, FVector(0, 1, 0)).GetSafeNormal();
	FVector Forward = FVector::CrossProduct(Right, Up).GetSafeNormal();

	int32 BaseIndex = Vertices.Num();
	
	for (int32 i = 0; i <= NumSides; ++i)
	{
		float Angle = FMath::DegreesToRadians(i * AngleStep);
		FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;

		Vertices.Add(Start + Offset);
		Vertices.Add(End + Offset);

		if (i < NumSides)
		{
			int32 Current = BaseIndex + (i * 2);
			int32 Next = BaseIndex + ((i + 1) * 2);

			// Triangle 1
			Triangles.Add(Current);
			Triangles.Add(Next);
			Triangles.Add(Current + 1);

			// Triangle 2
			Triangles.Add(Next);
			Triangles.Add(Next + 1);
			Triangles.Add(Current + 1);
		}
	}
}

