// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"

namespace DoobMathUtils {
	FVector GenerateRandomPerpendicularVector(const FVector& BaseVector);

	float EaseOutIncrement(int32 TotalSteps, int32 CurrentStep);
	float EaseInIncrement(int32 TotalSteps, int32 CurrentStep);

	float FloatChangeWithCurve(float StartRadius, float EndRadius, float CurrentRadius, int32 TotalSteps, int32 CurrentStep, int32 CurveType);

	void InterpolateBetweenBounds(float StartBound, float EndBound, int32 NumSegments, int32 SegmentIndex, float& StartValue, float& EndValue);
}