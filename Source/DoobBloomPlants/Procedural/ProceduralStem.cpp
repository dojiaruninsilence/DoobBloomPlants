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

	// Initialize vertex index
	int32 BaseIndex = 0;

	FVector CurrentPosition = FVector::ZeroVector;
	FVector CurrentDirection = FVector(0, 0, 1); // Initial growth direction up

	FVector RightVector = FVector(1, 0, 0);
	FVector UpVector = FVector(0, 1, 0);

	TArray<FVector> LastRingVertices;

	for (int32 i = 0; i < NumSegments; ++i)
	{
		// interpolate radius for this segment
		float CurrentRadius = FMath::Lerp(BaseRadius, TopRadius, static_cast<float>(i) / NumSegments);
		// float CurrentRadius = BaseRadius * FMath::Pow(TopRadius / BaseRadius, static_cast<float>(i) / NumSegments);

		// Determine the next segment's end pos
		FVector NextPosition = CurrentPosition + (CurrentDirection * SegmentLength);

		// Generate the ring at the start of this segment
		TArray<FVector> StartRingVertices;
		GenerateRing(CurrentPosition, CurrentDirection, UpVector, CurrentRadius, StartRingVertices);

		// Generate the ring at the end of this segment
		TArray<FVector> EndRingVertices;
		GenerateRing(NextPosition, CurrentDirection, UpVector, CurrentRadius, EndRingVertices);

		// connect the rings with triangles
		if (i > 0) // skip connecting for the first segment
		{
			ConnectRings(LastRingVertices, StartRingVertices, Vertices, Triangles, BaseIndex);
		}

		// Generate the cylinder for this segment
		ConnectRings(StartRingVertices, EndRingVertices, Vertices, Triangles, BaseIndex);
		
		// Update the current position and apply a random rotation
		LastRingVertices = EndRingVertices;
		CurrentPosition = NextPosition + (CurrentDirection * SegmentGapLength);

		// Apply random tilt
		FRotator RandomTilt = FRotator(FMath::RandRange(-30.0f, 30.0f), FMath::RandRange(-30.0f, 30.0f), 0);
		FQuat TiltQuat = FQuat(RandomTilt);

		CurrentDirection = TiltQuat.RotateVector(CurrentDirection).GetSafeNormal();
		RightVector = TiltQuat.RotateVector(RightVector).GetSafeNormal();
		UpVector = FVector::CrossProduct(CurrentDirection, RightVector).GetSafeNormal();

		// Ensure RightVector stays perpendicular to the new direction
		RightVector = FVector::CrossProduct(UpVector, CurrentDirection).GetSafeNormal();

		/*FRotator RandomDirection = FRotator(FMath::RandRange(-30.0f, 30.0f), FMath::RandRange(-30.0f, 30.0f), 0);
		CurrentDirection = RandomDirection.RotateVector(CurrentDirection).GetSafeNormal();*/
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
void AProceduralStem::GenerateCylinder(FVector Start, FVector End, float Radius, int32 LocalNumSides, TArray<FVector>& Vertices, TArray<int32>& Triangles)
{
	const float AngleStep = 360.0f / LocalNumSides;

	FVector Up = (End - Start).GetSafeNormal();
	FVector Right = FVector::CrossProduct(Up, FVector(0, 1, 0)).GetSafeNormal();
	FVector Forward = FVector::CrossProduct(Right, Up).GetSafeNormal();

	int32 BaseIndex = Vertices.Num();
	
	for (int32 i = 0; i <= LocalNumSides; ++i)
	{
		float Angle = FMath::DegreesToRadians(i * AngleStep);
		FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;

		Vertices.Add(Start + Offset);
		Vertices.Add(End + Offset);

		if (i < LocalNumSides)
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

void AProceduralStem::GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, TArray<FVector>& RingVertices)
{
	const float AngleStep = 360.0f / StemNumSides;
	FVector Right = FVector::CrossProduct(Direction, UpVector).GetSafeNormal();
	FVector Forward = FVector::CrossProduct(Right, Direction).GetSafeNormal();

	for (int i = 0; i <= StemNumSides; ++i)
	{
		float Angle = FMath::DegreesToRadians(i * AngleStep);
		FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;
		RingVertices.Add(Center + Offset);
	}
}

void AProceduralStem::ConnectRings(
	const TArray<FVector>& RingA,
	const TArray<FVector>& RingB,
	TArray<FVector>& Vertices,
	TArray<int32>& Triangles, 
	int32& BaseIndex
)
{
	for (int32 i = 0; i < RingA.Num(); ++i)
	{
		int32 CurrentA = BaseIndex + i;
		int32 NextA = BaseIndex + (i + 1) % RingA.Num();

		int32 CurrentB = BaseIndex + RingA.Num() + i;
		int32 NextB = BaseIndex + RingA.Num() + ((i + 1) % RingB.Num());

		// Triangle 1
		Triangles.Add(CurrentA);
		Triangles.Add(NextA);
		Triangles.Add(CurrentB);

		// Triangle 2
		Triangles.Add(NextA);
		Triangles.Add(NextB);
		Triangles.Add(CurrentB);
	}
	Vertices.Append(RingA);
	Vertices.Append(RingB);
	BaseIndex += RingA.Num() + RingB.Num();
}