// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMainStem.generated.h"

class AProceduralStem;
class AProceduralStemNode;

UCLASS()
class DOOBBLOOMPLANTS_API AProceduralMainStem : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AProceduralMainStem();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;


public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	// length params
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	int32 NumStems = 3;

	// radius params
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float BaseRadius = 20.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float TopRadius = 0.0f;

	// Direction params
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float GrowTowardProbability = 0.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float GrowAwayProbability = 0.2f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float GrowAwayAmount = 0.25;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float GrowTowardAmount = 0.25;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	int32 GrowCurveType = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	float Randomness = 0.3f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	FVector TargetPoint = FVector(0, 0, 1);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Main Stem")
	FVector PerpVector = FVector(0, 1, 0);

	UPROPERTY(BlueprintReadWrite, Category = "Main Stem")
	FVector StartUpVector = FVector(0, 1, 0);

	void GenerateMainStem();

	void GenerateMainStemChain();

private:

	// Root component for transform handling
	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Components", meta = (AllowPrivateAccess = "true"))
	USceneComponent* Root;

	UPROPERTY()
	AProceduralStem* MainStem;

	UPROPERTY()
	AProceduralStemNode* MainStemNode;

	UPROPERTY()
	TArray<AProceduralStem*> Stems;

	UPROPERTY()
	TArray<AProceduralStemNode*> Nodes;
};
