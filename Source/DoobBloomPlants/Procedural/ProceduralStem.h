// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"

#include "ProceduralStem.generated.h"

// forward declarations
class AProceduralStemNode;

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
	// position params
	UPROPERTY(BlueprintReadWrite, Category = "Stem")
	FVector StartPosition = FVector::ZeroVector;

	// length params
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	int32 NumSegments = 10;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float SegmentLength = 100.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float SegmentGapLength = 20.0f;

	// Radius params
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	int32 TaperType = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float BaseRadius = 20.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float TopRadius = 5.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	int32 StemNumSides = 12; // Default to 12

	// Direction params
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float GrowTowardProbability = 0.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float GrowAwayProbability = 0.2f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float GrowAwayAmount = 0.25;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float GrowTowardAmount = 0.25;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	int32 GrowCurveType = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	float Randomness = 0.3f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Stem")
	FVector TargetPoint = FVector(0, 0, 1);

	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	FVector StartUpVector = FVector(0, 1, 0);

	// For Connecting to Nodes
	UPROPERTY(BlueprintReadWrite, Category = "Stem")
	AProceduralStemNode* ParentNode;

	UPROPERTY(BlueprintReadWrite, Category = "Stem")
	AProceduralStemNode* EndNode;

	UFUNCTION(BlueprintCallable, Category = "Stem")
	void SetParentNode(AProceduralStemNode* NewParentNode);

	UFUNCTION(BlueprintCallable, Category = "Stem")
	void SetEndNode(AProceduralStemNode* NewEndNode);

	UFUNCTION(BlueprintCallable, Category = "Stem")
	AProceduralStemNode* GetParentNode();

	UFUNCTION(BlueprintCallable, Category = "Stem")
	AProceduralStemNode* GetEndNode();

	// Generate the stem geometry
	UFUNCTION(BlueprintCallable, Category = "Stem")
	void GenerateStem();

	FVector GetStartDirection();
	void SetStartDirection(FVector Direcetion);
	FVector GetEndPosition();
	FVector GetEndDirection();
	FVector GetEndUpVector();
	float GetEndRadius();

private:

protected:
	FVector StartDirection = FVector(0, 0, 1);
	FVector EndPosition;
	FVector EndDirection;
	FVector EndUpVector;
	float EndRadius;
};
