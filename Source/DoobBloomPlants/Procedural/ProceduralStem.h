// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"
#include "ProceduralStem.generated.h"

UCLASS()
class DOOBBLOOMPLANTS_API AProceduralStem : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AProceduralStem();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	// Procedural mesh component for the stem
	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Mesh")
	UProceduralMeshComponent* StemMesh;

	// Params for stem generation
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	int32 NumSegments = 10;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float SegmentLength = 100.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float SegmentRadius = 10.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float BaseRadius = 20.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float TopRadius = 5.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	int32 StemNumSides = 12; // Default to 12

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float SegmentGapLength = 20.0f;

	// Generate the stem geometry
	void GenerateStem();
	void GenerateCylinder(FVector Start, FVector End, float Radius, int32 NumSides, TArray<FVector>& Vertices, TArray<int32>& Triangles);

protected:
	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, TArray<FVector>& RingVertices);
	void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex);
	void CreateTransitionGeometry(const TArray<FVector>& LastRingVertices, const TArray<FVector>& NextRingVertices, TArray<FVector>& Vertices, TArray<int32>& Triangles);
};
