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
	float SegmentRadius = 10.f;

	// Generate the stem geometry
	void GenerateStem();
	void GenerateCylinder(FVector Start, FVector End, float Radus, TArray<FVector>& Vertices, TArray<int32>& Triangles);

};
