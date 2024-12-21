// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralMainStem.h"
#include "ProceduralStem.h"
#include "ProceduralStemNode.h"

// Sets default values
AProceduralMainStem::AProceduralMainStem()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	// Create and set the root component
	Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	SetRootComponent(Root);

	// initialize member variables
	MainStem = nullptr;
	MainStemNode = nullptr;
}

// Called when the game starts or when spawned
void AProceduralMainStem::BeginPlay()
{
	Super::BeginPlay();

	GenerateMainStem();
	
}

// Called every frame
void AProceduralMainStem::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

void AProceduralMainStem::GenerateMainStem()
{
	// spawn stem as child actor
	MainStem = GetWorld()->SpawnActor<AProceduralStem>(AProceduralStem::StaticClass());
	if (MainStem)
	{
		// attach the stem to this actor
		MainStem->AttachToComponent(Root, FAttachmentTransformRules::KeepRelativeTransform);

		// generate the stem geometry
		MainStem->GenerateStem();

		// retrieve the end position and direction of the stem
		FVector StemEndPosition = MainStem->GetEndPosition();
		FVector StemEndDirection = MainStem->GetEndDirection();

		// Convert the end position from local space to world space (relative to the root)
		FVector StemEndPositionWorld = RootComponent->GetComponentTransform().TransformPosition(StemEndPosition);

		UE_LOG(LogTemp, Log, TEXT("Generate Main Stem The End position is: %s, The end direction is: %s, the world end position is: %s"), *StemEndPosition.ToString(), *StemEndDirection.ToString(), *StemEndPositionWorld.ToString());

		// spawn the node at the end of the stem
		MainStemNode = GetWorld()->SpawnActor<AProceduralStemNode>(AProceduralStemNode::StaticClass(), StemEndPosition, FRotator::ZeroRotator);
		if (MainStemNode)
		{
			// attach the node to the stem
			MainStemNode->AttachToComponent(Root, FAttachmentTransformRules::KeepRelativeTransform);

			// set node params
			MainStemNode->StartPosition = FVector::ZeroVector;
			MainStemNode->StartDirection = StemEndDirection;
			MainStemNode->StartRadius = MainStem->GetEndRadius();
			MainStemNode->StartUpVector = MainStem->GetEndUpVector();

			// generate the node geometry
			MainStemNode->GenerateNode();

			// Link the node back to the stem
			MainStem->SetEndNode(MainStemNode);
		}
	}
}