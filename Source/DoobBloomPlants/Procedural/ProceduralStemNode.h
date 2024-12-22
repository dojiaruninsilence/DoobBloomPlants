// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"

#include "ProceduralStemNode.generated.h"

// forward declarations
class AProceduralStem;

UCLASS()
class DOOBBLOOMPLANTS_API AProceduralStemNode : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AProceduralStemNode();	

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Mesh")
	UProceduralMeshComponent* NodeMesh;

	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	AProceduralStem* ConnectedStem;

	// Start positioning
	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	FVector StartPosition;

	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	FVector StartDirection;

	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	FVector StartUpVector;

	// ring parameters
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem Node")
	float StartRadius = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem Node")
	float EndRadius = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem Node")
	float RadiusChangeAmount = 30.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem Node")
	int32 StemNodeNumSides = 12;

	// length parameters
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem Node")
	int32 NodeSegments = 20;

	void AttachToStem(AProceduralStem* Stem);
	void DetachFromStem(AProceduralStem* Stem);

	void GenerateNode();

	FVector GetEndPosition();
	FVector GetEndDirection();
	FVector GetEndUpVector();

protected:
	FVector EndPosition;
	FVector EndDirection;
	FVector EndUpVector;
};

