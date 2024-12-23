// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralStem.h"
#include "ProceduralStemNode.h"

#include "GeometryUtilities.h"
#include "MathUtililities.h"
#include "DoobMeshUtils.h"

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
	TArray<TArray<FVector>> Rings;

	TArray<FVector> Vertices;
	TArray<int32> Triangles;

	// Initialize vertex index
	int32 BaseIndex = 0;

	FVector CurrentPosition = StartPosition;
	FVector NextPosition = StartPosition;

	FVector CurrentDirection = StartDirection;
	FVector LastDirection = StartDirection;

	FVector RightVector = FVector(1, 0, 0);
	FVector UpVector = StartUpVector;

	TArray<FVector> LastRingVertices;

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

		// generate joint ring
		TArray<FVector> BendRingVertices;
		if (i > 0) {
			FVector BendPosition = (CurrentPosition + NextPosition) * 0.5f;
			FVector BendDirection = (CurrentDirection + LastDirection).GetSafeNormal();
			float BendRadius = CurrentRadius + (CurrentRadius * 0.05f);
			GeometryUtilities::GenerateRing(BendPosition, BendDirection, UpVector, BendRadius, StemNumSides, BendRingVertices);
			Rings.Add(BendRingVertices);
		}

		LastDirection = CurrentDirection;

		// Determine the next segment's end pos
		NextPosition = CurrentPosition + (CurrentDirection * SegmentLength);

		UE_LOG(LogTemp, Log, TEXT("Generate Stem iteration num: %d, The current position is: %s, The next position is: %s, the current direction is: %s"), i, *CurrentPosition.ToString(), *NextPosition.ToString(), *CurrentDirection.ToString());

		EndDirection = CurrentDirection;
		EndRadius = CurrentRadius;
		EndUpVector = UpVector;

		// Generate the ring at the start of this segment
		TArray<FVector> StartRingVertices;
		GeometryUtilities::GenerateRing(CurrentPosition, CurrentDirection, UpVector, CurrentRadius, StemNumSides, StartRingVertices);
		Rings.Add(StartRingVertices);

		// Generate the ring at the end of this segment
		TArray<FVector> EndRingVertices;
		GeometryUtilities::GenerateRing(NextPosition, CurrentDirection, UpVector, CurrentRadius, StemNumSides, EndRingVertices);
		Rings.Add(EndRingVertices);

		// connect the rings with triangles
		if (i > 0) // skip connecting for the first segment
		{
			//GeometryUtilities::ConnectRings(LastRingVertices, BendRingVertices, Vertices, Triangles, BaseIndex);
			//GeometryUtilities::ConnectRings(BendRingVertices, StartRingVertices, Vertices, Triangles, BaseIndex);
		}

		// Generate the cylinder for this segment
		//GeometryUtilities::ConnectRings(StartRingVertices, EndRingVertices, Vertices, Triangles, BaseIndex);
		
		// Update the current position and apply a random rotation
		LastRingVertices = EndRingVertices;
		//CurrentPosition = NextPosition + (CurrentDirection * SegmentGapLength);

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

		CurrentPosition = NextPosition + (CurrentDirection * SegmentGapLength);

		// update right and up vectors
		RightVector = FVector::CrossProduct(UpVector, CurrentDirection).GetSafeNormal();
		UpVector = FVector::CrossProduct(CurrentDirection, RightVector).GetSafeNormal();

	}


	// store the end position and direction
	EndPosition = NextPosition;

	

	// smooth vertices
	DoobMeshUtils::SmoothMeshVertices(Vertices, Triangles, 1);

	GeometryUtilities::ConnectRingArray(Rings, Vertices, Triangles, BaseIndex);

	TArray<FVector> Normals;
	Normals.Init(FVector::UpVector, Vertices.Num());

	DoobMeshUtils::RemoveDegenerateTriangles(Vertices, Triangles);

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

FVector AProceduralStem::GetStartDirection()
{
	return StartDirection;
}

void AProceduralStem::SetStartDirection(FVector Direcetion) {
	StartDirection = Direcetion;
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

