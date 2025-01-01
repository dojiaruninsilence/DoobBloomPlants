// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"

#include "DoobProfileUtils.h"
#include "DoobGeometryUtils.h"
#include "DoobMeshUtils.h"

#include "ProceduralTubularJointNode.generated.h"

UCLASS()
class DOOBBLOOMPLANTS_API AProceduralTubularJointNode : public AActor
{
	GENERATED_BODY()

public:
	// Sets default values for this actor's properties
	AProceduralTubularJointNode();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	// function to build and render the main tube
	void BuildMainTube();

	// Start positioning
	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	FVector StartPosition = FVector::ZeroVector;

private:
	UPROPERTY(VisibleAnywhere)
	UProceduralMeshComponent* TubularJointNodeMesh;

	//UPROPERTY(EditAnywhere)
	DoobProfileUtils::F2DProfile MainTubeProfile;

	UPROPERTY(EditAnywhere)
	FTransform MainTubeTransform;

	UPROPERTY(EditAnywhere)
	int32 MainTubeSegments;

	UPROPERTY(EditAnywhere)
	float MainTubeRadius;

	UPROPERTY(EditAnywhere)
	bool bIsMainTubeClosed;

	//void CreateProceduralMesh();

	FVector EndPosition;
};