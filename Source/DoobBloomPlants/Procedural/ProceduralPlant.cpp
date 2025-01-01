// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralPlant.h"

#include "ProceduralMainStem.h"

// Sets default values
AProceduralPlant::AProceduralPlant()
{
	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	// Create and set the root component
	Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	SetRootComponent(Root);
}

// Called when the game starts or when spawned
void AProceduralPlant::BeginPlay()
{
	Super::BeginPlay();

	GeneratePlant();

}

// Called every frame
void AProceduralPlant::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}


void AProceduralPlant::GeneratePlant() {
	MainStem = GetWorld()->SpawnActor<AProceduralMainStem>(AProceduralMainStem::StaticClass());

	if (MainStem) {
		// attach the stem to this actor
		MainStem->AttachToComponent(Root, FAttachmentTransformRules::KeepRelativeTransform);

		// generate the stem geometry
		MainStem->GenerateMainStemChain();
	}
}