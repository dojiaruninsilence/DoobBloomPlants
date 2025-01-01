// Fill out your copyright notice in the Description page of Project Settings.

using System.IO;
using UnrealBuildTool;

public class DoobBloomPlants : ModuleRules
{
	public DoobBloomPlants(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;
	
		PublicDependencyModuleNames.AddRange(new string[] { "Core", "CoreUObject", "Engine", "InputCore", "ProceduralMeshComponent", "AutomationController" });

		PrivateDependencyModuleNames.AddRange(new string[] {  });

        // Add include paths for Utilities and other custom folders
        PublicIncludePaths.AddRange(new string[]
        {
            Path.Combine(ModuleDirectory, "DoobUtils"),
            Path.Combine(ModuleDirectory, "Procedural"),
            Path.Combine(ModuleDirectory, "Tests")
        });

        // Uncomment if you are using Slate UI
        // PrivateDependencyModuleNames.AddRange(new string[] { "Slate", "SlateCore" });

        // Uncomment if you are using online features
        // PrivateDependencyModuleNames.Add("OnlineSubsystem");

        // To include OnlineSubsystemSteam, add it to the plugins section in your uproject file with the Enabled attribute set to true
    }
}
