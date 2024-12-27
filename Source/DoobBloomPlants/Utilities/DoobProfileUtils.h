#pragma once

#include "CoreMinimal.h"
#include "Math/Vector2D.h"
#include "Containers/Array.h"

namespace DoobProfileUtils {

	struct F2DProfile {
		TArray<FVector2D> Points; // points defining the profile
		bool bIsClosed; // whether the profile is closed, looped

		F2DProfile();
		F2DProfile(TArray<FVector2D> InPoints, bool InClosed = true);

		// scales the profile by a given factor
		void Scale(float Factor);

		// translate the profile by a given offset
		void Translate(const FVector2D& Offset);

		// rotates the profile around the origin by a given angle in radians
		void Rotate(float Angle);
	};

	F2DProfile GenerateLinearProfile(int32 NumSegments, float StartHeight, float EndHeight);
	F2DProfile GenerateEggShapedCylinderProfile(int32 NumSegments, float StartHeight, float MidHeight, float EndHeight, float BulgePosition);

	FVector2D EvaluateProfile(const F2DProfile& Profile, float t);

	float SumYValuesInProfile(const F2DProfile& Profile);
	float MaxYValueInProfile(const F2DProfile& Profile);
}