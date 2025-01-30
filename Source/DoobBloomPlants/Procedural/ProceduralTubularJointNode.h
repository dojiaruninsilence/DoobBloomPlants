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

	UPROPERTY(BlueprintReadWrite, Category = "Stem Node")
	FVector LateralStartPosition = FVector::ZeroVector;

	// direction params
	FVector MainTubeDirection;

	FVector LateralTubeDirection;

	// length params
	float MainTubeLength;
	float LateralTubeLength;

	UPROPERTY(EditAnywhere)
	int32 MainTubeSegments;

	UPROPERTY(EditAnywhere)
	int32 LateralTubeSegments;

	// radius params
	UPROPERTY(EditAnywhere)
	float MainTubeStartRadius;

	UPROPERTY(EditAnywhere)
	float MainTubeMidRadius;

	UPROPERTY(EditAnywhere)
	float MainTubeEndRadius;

	UPROPERTY(EditAnywhere)
	float LateralTubeRadius;

	int32 MainTubeNumSides;
	int32 LateralTubeNumSides;

private:
	UPROPERTY(VisibleAnywhere)
	UProceduralMeshComponent* TubularJointNodeMesh;

	//UPROPERTY(EditAnywhere)
	DoobProfileUtils::F2DProfile MainTubeProfile;
	DoobProfileUtils::F2DProfile LateralTubeProfile;

	UPROPERTY(EditAnywhere)
	FTransform MainTubeTransform;

	UPROPERTY(EditAnywhere)
	FTransform LateralTubeTransform;	

	UPROPERTY(EditAnywhere)
	bool bIsMainTubeClosed;

	UPROPERTY(EditAnywhere)
	bool bIsLateralTubeClosed;

	//void CreateProceduralMesh();

	FVector EndPosition;
	FVector LateralEndPosition;
};