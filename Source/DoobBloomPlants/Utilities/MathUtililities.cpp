// Fill out your copyright notice in the Description page of Project Settings.


#include "MathUtililities.h"

namespace MathUtilities {
	FVector GenerateRandomPerpendicularVector(const FVector& BaseVector)
	{
		// ensure non zero length
		if (!BaseVector.IsNearlyZero())
		{
			// choose an arbitrary vector dynamically
			FVector ArbitraryVector;

			// avoid parallel vectors by choosing one with non-zero cross product
			if (FMath::Abs(BaseVector.X) <= FMath::Abs(BaseVector.Y) && FMath::Abs(BaseVector.X) <= FMath::Abs(BaseVector.Z))
			{
				ArbitraryVector = FVector(1, 0, 0); // use x axis if base vector is predominantly in x
			}
			else if (FMath::Abs(BaseVector.Y) <= FMath::Abs(BaseVector.Z))
			{
				ArbitraryVector = FVector(0, 1, 0); // use y axis if mostly in y
			}
			else
			{
				ArbitraryVector = FVector(0, 0, 1); // otherwise use z
			}

			// calculate the perpendicular vector using the cross product
			FVector PerpendicularVector = FVector::CrossProduct(BaseVector, ArbitraryVector).GetSafeNormal();

			// add randomness by rotating around basevactor
			float RandomAngle = FMath::FRandRange(0.0f, 360.0f); // random angle in degrees
			FQuat Rotation = FQuat(BaseVector.GetSafeNormal(), FMath::DegreesToRadians(RandomAngle));
			FVector RandomizedPerpendicular = Rotation.RotateVector(PerpendicularVector);

			return RandomizedPerpendicular.GetSafeNormal();
		}

		// Fallback to prevent errors if BaseVector is invalid
		return FVector::ZeroVector;
	}

	float EaseOutIncrement(int32 TotalSteps, int32 CurrentStep)
	{
		// normalize the step (from 0 to 1)
		float t = (float)CurrentStep / (float)TotalSteps;
		// ease out curve for increment scaling
		return FMath::Lerp(1.0f, 0.0f, t);
	}

	float EaseInIncrement(int32 TotalSteps, int32 CurrentStep)
	{
		// normalize the current step / total steps
		float normalizedStep = FMath::Clamp((float)CurrentStep / (float)TotalSteps, 0.0f, 1.0f);

		return normalizedStep * normalizedStep;
	}

	float FloatChangeWithCurve(float StartRadius, float EndRadius, float CurrentRadius, int32 TotalSteps, int32 CurrentStep, int32 CurveType)
	{
		float ReturnRadius = CurrentRadius;

		// calc the scale
		float scale = 0.0f;
		float ChangeAmount = 0.0f;

		if (CurveType == 0)
		{
			scale = EaseOutIncrement(TotalSteps, CurrentStep);
			// calc the change
			ChangeAmount = (EndRadius - StartRadius) * scale / TotalSteps;
		}
		else if (CurveType == 1)
		{
			scale = EaseInIncrement(TotalSteps, CurrentStep);
			// calc the change
			ChangeAmount = (EndRadius - StartRadius) * scale;
		}

		// get new radius
		ReturnRadius += ChangeAmount;

		// this is super strange, cant remove with out messing up code, need to rework everything soon
		UE_LOG(LogTemp, Log, TEXT("Need to figure out why i cant remove this"));

		// make sure not over the end

		if (ChangeAmount < 0)
		{
			ReturnRadius = FMath::Clamp(ReturnRadius, EndRadius, StartRadius);
		}
		else
		{
			ReturnRadius = FMath::Clamp(ReturnRadius, StartRadius, EndRadius);
		}

		return ReturnRadius;

	}
}