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

	void GenerateMainStem();

private:

	// Root component for transform handling
	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Components", meta = (AllowPrivateAccess = "true"))
	USceneComponent* Root;

	UPROPERTY()
	AProceduralStem* MainStem;

	UPROPERTY()
	AProceduralStemNode* MainStemNode;
};
