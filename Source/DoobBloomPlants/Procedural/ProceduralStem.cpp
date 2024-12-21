// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralStem.h"
#include "ProceduralStemNode.h"

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
	//GenerateStem();
	
}

// Called every frame
void AProceduralStem::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

void AProceduralStem::SetParentNode(AProceduralStemNode* NewParentNode)
{
	if (NewParentNode)
	{
		ParentNode = NewParentNode;
	}
	else
	{
		UE_LOG(LogTemp, Warning, TEXT("Attempted to set a null parent node."));
	}
}

void AProceduralStem::SetEndNode(AProceduralStemNode* NewEndNode)
{
	if (NewEndNode)
	{
		// if there was a previous node, detach it
		if (EndNode)
		{
			EndNode->DetachFromStem(this);
		}

		// attach new end node
		EndNode = NewEndNode;
		EndNode->AttachToStem(this); // need to add to node class first

		// update the procedural mesh if needed
		//GenerateStem(); // Rebuild Stem to reflect new endpoint connection

		UE_LOG(LogTemp, Log, TEXT("End node successfully set to: %s"), *NewEndNode->GetName());
	}
	else
	{
		UE_LOG(LogTemp, Warning, TEXT("Attempted to set a null end node."));
	}
}

AProceduralStemNode* AProceduralStem::GetParentNode()
{
	return ParentNode;
}

AProceduralStemNode* AProceduralStem::GetEndNode()
{
	return EndNode;
}

// Generate the full stem
void AProceduralStem::GenerateStem()
{
	TArray<FVector> Vertices;
	TArray<int32> Triangles;

	// Initialize vertex index
	int32 BaseIndex = 0;

	FVector CurrentPosition = FVector::ZeroVector;
	FVector NextPosition = FVector::ZeroVector;
	FVector CurrentDirection = FVector(0, 0, 1); // Initial growth direction up

	FVector RightVector = FVector(1, 0, 0);
	FVector UpVector = FVector(0, 1, 0);

	FVector PerpVector = GenerateRandomPerpendicularVector(TargetPoint);

	TArray<FVector> LastRingVertices;

	// store the start position and direction
	StartPosition = CurrentPosition;
	StartDirection = CurrentDirection;

	for (int32 i = 0; i < NumSegments; ++i)
	{
		// interpolate radius for this segment
		float CurrentRadius;
		if (TaperType == 0)
		{
			CurrentRadius = FMath::Lerp(BaseRadius, TopRadius, static_cast<float>(i) / NumSegments);
		}
		else if (TaperType == 1)
		{
			CurrentRadius = BaseRadius * FMath::Pow(TopRadius / BaseRadius, static_cast<float>(i) / NumSegments);
		}
		else
		{
			CurrentRadius = BaseRadius;
		}

		// Determine the next segment's end pos
		NextPosition = CurrentPosition + (CurrentDirection * SegmentLength);

		UE_LOG(LogTemp, Log, TEXT("Generate Stem iteration num: %d, The current position is: %s, The next position is: %s, the current direction is: %s"), i, *CurrentPosition.ToString(), *NextPosition.ToString(), *CurrentDirection.ToString());

		EndDirection = CurrentDirection;
		EndRadius = CurrentRadius;
		EndUpVector = UpVector;

		// Generate the ring at the start of this segment
		TArray<FVector> StartRingVertices;
		GenerateRing(CurrentPosition, CurrentDirection, UpVector, CurrentRadius, StemNumSides, StartRingVertices);

		// Generate the ring at the end of this segment
		TArray<FVector> EndRingVertices;
		GenerateRing(NextPosition, CurrentDirection, UpVector, CurrentRadius, StemNumSides, EndRingVertices);

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
		// probability check to determine growth behavior
		float RandomValue = FMath::FRand(); // rand val between 0.0 and 1.0
		
		//Determine the bias direction based on probabilities
		FVector BiasDirection;
		if (RandomValue < GrowTowardProbability)
		{
			// grow toward the target point
			BiasDirection = FMath::Lerp(CurrentDirection, TargetPoint, GrowTowardAmount).GetSafeNormal();
		}
		else if (RandomValue < GrowTowardProbability + GrowAwayProbability)
		{			
			if (GrowCurveType == 1 && i > (NumSegments / 3) * 2)
			{
				BiasDirection = FMath::Lerp(CurrentDirection, FVector(0, 0, -1), GrowAwayAmount).GetSafeNormal();
			}

			else if (GrowCurveType == 2 && i > NumSegments / 2)
			{
				BiasDirection = FMath::Lerp(CurrentDirection, FVector(0, 0, -1), GrowAwayAmount).GetSafeNormal();
			}

			else
			{
				//  Grow away from the target point
				BiasDirection = FMath::Lerp(CurrentDirection, PerpVector, GrowAwayAmount).GetSafeNormal();
			}			
		}
		else
		{
			// completely random growth
			FVector TempRandomDirection = FVector(FMath::FRandRange(-1.0f , 1.0f), FMath::FRandRange(-1.0f, 1.0f), FMath::FRandRange(-1.0f, 1.0f)).GetSafeNormal();
			BiasDirection = FMath::Lerp(CurrentDirection, TempRandomDirection, Randomness).GetSafeNormal();
		}

		// apply randomness to the bias direction
		FRotator RandomTilt = FRotator(
			FMath::RandRange(-Randomness * 30.0f, Randomness * 30.0f),
			FMath::RandRange(-Randomness * 30.0f, Randomness * 30.0f),
			0
		);
		FQuat TiltQuat = FQuat(RandomTilt);

		// combine the bias direction and randomness
		CurrentDirection = TiltQuat.RotateVector(BiasDirection).GetSafeNormal();

		// update right and up vectors
		RightVector = FVector::CrossProduct(UpVector, CurrentDirection).GetSafeNormal();
		UpVector = FVector::CrossProduct(CurrentDirection, RightVector).GetSafeNormal();

	}


	// store the end position and direction
	EndPosition = NextPosition;

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

FVector AProceduralStem::GetStartPosition()
{
	return StartPosition;
}

FVector AProceduralStem::GetStartDirection()
{
	return StartDirection;
}

FVector AProceduralStem::GetEndPosition()
{
	return EndPosition;
}

FVector AProceduralStem::GetEndDirection()
{
	return EndDirection;
}

FVector AProceduralStem::GetEndUpVector()
{
	return EndUpVector;
}

float AProceduralStem::GetEndRadius()
{
	return EndRadius;
}


FVector AProceduralStem::GenerateRandomPerpendicularVector(const FVector& BaseVector)
{
	// ensure non zero length
	if (!BaseVector.IsNearlyZero())
	{
		// choose an arbitrary vector dynamically
		FVector ArbitraryVector;

		// avoid parallel vectors by choosing one with non-zero cross product
		if (FMath::Abs(BaseVector.X) <= FMath::Abs(BaseVector.Y) && FMath::Abs(BaseVector.X) <= FMath::Abs(BaseVector.Z))
		{
			ArbitraryVector = FVector(1, 0, 0); // use x axis if base vector is predominantly in x
		} 
		else if (FMath::Abs(BaseVector.Y) <= FMath::Abs(BaseVector.Z))
		{
			ArbitraryVector = FVector(0, 1, 0); // use y axis if mostly in y
		}
		else
		{
			ArbitraryVector = FVector(0, 0, 1); // otherwise use z
		}

		// calculate the perpendicular vector using the cross product
		FVector PerpendicularVector = FVector::CrossProduct(BaseVector, ArbitraryVector).GetSafeNormal();

		// add randomness by rotating around basevactor
		float RandomAngle = FMath::FRandRange(0.0f, 360.0f); // random angle in degrees
		FQuat Rotation = FQuat(BaseVector.GetSafeNormal(), FMath::DegreesToRadians(RandomAngle));
		FVector RandomizedPerpendicular = Rotation.RotateVector(PerpendicularVector);

		return RandomizedPerpendicular.GetSafeNormal();
	}

	// Fallback to prevent errors if BaseVector is invalid
	return FVector::ZeroVector;
}

void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, TArray<FVector>& RingVertices)
{
	const float AngleStep = 360.0f / NumSides;
	FVector Right = FVector::CrossProduct(Direction, UpVector).GetSafeNormal();
	FVector Forward = FVector::CrossProduct(Right, Direction).GetSafeNormal();

	for (int i = 0; i <= NumSides; ++i)
	{
		float Angle = FMath::DegreesToRadians(i * AngleStep);
		FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;
		RingVertices.Add(Center + Offset);
	}
}

void ConnectRings(
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