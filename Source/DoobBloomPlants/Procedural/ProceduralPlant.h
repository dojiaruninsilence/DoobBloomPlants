// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralPlant.generated.h"

class AProceduralMainStem;

UCLASS()
class DOOBBLOOMPLANTS_API AProceduralPlant : public AActor
{
	GENERATED_BODY()

public:
	// Sets default values for this actor's properties
	AProceduralPlant();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	void GeneratePlant();

private:
	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Components", meta = (AllowPrivateAccess = "true"))
	USceneComponent* Root;

	UPROPERTY()
	AProceduralMainStem* MainStem;

	UPROPERTY()
	TArray<AProceduralMainStem*> Stems;

};