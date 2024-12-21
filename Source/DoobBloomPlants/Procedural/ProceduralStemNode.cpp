// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralStemNode.h"
#include "ProceduralStem.h"

// Sets default values
AProceduralStemNode::AProceduralStemNode()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	//PrimaryActorTick.bCanEverTick = true;

	// Create the Procedural Mesh Component
	NodeMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("StemMesh"));
	RootComponent = NodeMesh;

	// Enable Collision
	NodeMesh->bUseAsyncCooking = true;
}

// Called when the game starts or when spawned
void AProceduralStemNode::BeginPlay()
{
	Super::BeginPlay();
	//GenerateNode()
	
}

// Called every frame
void AProceduralStemNode::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

void AProceduralStemNode::AttachToStem(AProceduralStem* Stem)
{
	if (Stem)
	{
		// store a reference to the connected stem
		ConnectedStem = Stem;
		UE_LOG(LogTemp, Log, TEXT("Node %s attached to stem %s"), *GetName(), *Stem->GetName());
	}
}

void AProceduralStemNode::DetachFromStem(AProceduralStem* Stem)
{
	if (ConnectedStem == Stem)
	{
		ConnectedStem = nullptr;
		UE_LOG(LogTemp, Log, TEXT("Node %s detached from stem %s"), *GetName(), *Stem->GetName());
	}
}

void AProceduralStemNode::GenerateNode()
{
	// validate node mesh
	if (!NodeMesh)
	{
		UE_LOG(LogTemp, Error, TEXT("NodeMesh is not initialized in GenerateNode!"));
		return;
	}

	// ensure non zero base radius
	if (StartRadius <= 0)
	{
		UE_LOG(LogTemp, Error, TEXT("Invalid BaseRadius: %f"), StartRadius);
		return;
	}



	TArray<FVector> Vertices;
	TArray<int32> Triangles;

	// initialize vertex index
	int32 BaseIndex = 0;

	FVector CurrentPosition = StartPosition;
	FVector NextPosition = StartPosition;
	FVector CurrentDirection = StartDirection.IsNormalized() ? StartDirection : StartDirection.GetSafeNormal();

	FVector RightVector = FVector(1, 0, 0);
	FVector UpVector = StartUpVector;

	TArray<FVector> LastRingVertices;

	float CurrentRadius = StartRadius;
	float NextRadius = StartRadius;
	float RadiusModified = StartRadius + RadiusChangeAmount;

	// number of segments
	int32 NodeSegments = 20;

	for (int32 i = 0; i <= NodeSegments; ++i)
	{
		// interpolate radius to form a spherical bulge
		float Fraction = static_cast<float>(i) / NodeSegments;
		CurrentRadius = NextRadius;

		CurrentPosition = NextPosition;

		// generate the ring at this position
		TArray<FVector> CurrentRingVertices;
		GenerateRing(CurrentPosition, CurrentDirection, UpVector, CurrentRadius, StemNodeNumSides, CurrentRingVertices);

		// Connect the Previous ring to the current ring
		if (i > 0)
		{
			ConnectRings(LastRingVertices, CurrentRingVertices, Vertices, Triangles, BaseIndex);
		}

		// update variables for the next iteration
		LastRingVertices = CurrentRingVertices;

		// calc the position for this ring
		NextPosition = CurrentPosition + (CurrentDirection * (StartRadius * 0.5f));

		int32 curveType = 0;
		int32 steps = NodeSegments / 2;
		if (i <= steps)
		{			
			NextRadius = RadiusChange(StartRadius, RadiusModified, CurrentRadius, steps, i, curveType);
		}
		/*else if (i == NodeSegments - 1)
		{
			NextRadius = EndRadius;
		}*/
		else
		{
			int32 step = i - steps;
			NextRadius = RadiusChange(RadiusModified, EndRadius, CurrentRadius, steps, step, curveType+1);
			UE_LOG(LogTemp, Log, TEXT("Node Segment step %d"), step);
		}

		UE_LOG(LogTemp, Log, TEXT("Node Segment: %d, Fraction: %f, CurrentRadius: %f"), i, Fraction, CurrentRadius);

	}

	if (Vertices.Num() == 0 || Triangles.Num() == 0)
	{
		UE_LOG(LogTemp, Error, TEXT("No vertices or triangles generated."));
		return;
	}

	// Finalize the mesh
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

	UE_LOG(LogTemp, Log, TEXT("Vertices: %d, Normals: %d, UV0: %d, VertexColors: %d, Tangents: %d"),
		Vertices.Num(), Normals.Num(), UV0.Num(), VertexColors.Num(), Tangents.Num());


	// Create the node mesh section
	NodeMesh->CreateMeshSection(0, Vertices, Triangles, Normals, UV0, VertexColors, Tangents, true);
}

float EaseOutIncrement(int32 TotalSteps, int32 CurrentStep)
{
	// normalize the step (from 0 to 1)
	float t = (float)CurrentStep / (float)TotalSteps;
	// ease out curve for increment scaling
	return FMath::Lerp(1.0f, 0.0f, t);
}

float EaseInIncrement(int32 TotalSteps, int32 CurrentStep)
{
	// normalize the current step / total steps
	float normalizedStep = FMath::Clamp((float)CurrentStep / (float)TotalSteps, 0.0f, 1.0f);

	return normalizedStep * normalizedStep;
}

float RadiusChange(float StartRadius, float EndRadius, float CurrentRadius, int32 TotalSteps, int32 CurrentStep, int32 CurveType)
{
	float ReturnRadius = CurrentRadius;

	// calc the scale
	float scale = 0.0f;
	float ChangeAmount = 0.0f;

	if (CurveType == 0)
	{
		scale = EaseOutIncrement(TotalSteps, CurrentStep);
		// calc the change
		ChangeAmount = (EndRadius - StartRadius) * scale / TotalSteps;
	}
	else if (CurveType == 1)
	{
		scale = EaseInIncrement(TotalSteps, CurrentStep);
		// calc the change
		ChangeAmount = (EndRadius - StartRadius) * scale;
		UE_LOG(LogTemp, Log, TEXT("Node Segment scale %f"), scale);
	}

	// get new radius
	ReturnRadius += ChangeAmount;
	UE_LOG(LogTemp, Log, TEXT("Node Segment start radius: %f, end radius: %f, change: %f, return: %f"), StartRadius, EndRadius, ChangeAmount, ReturnRadius);

	// make sure not over the end

	if (ChangeAmount < 0)
	{
		ReturnRadius = FMath::Clamp(ReturnRadius, EndRadius, StartRadius);
	}
	else
	{
		ReturnRadius = FMath::Clamp(ReturnRadius, StartRadius, EndRadius);
	}

	return ReturnRadius;

}