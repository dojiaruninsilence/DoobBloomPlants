// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"

/**
 * @namespace DoobMathUtils
 * Contains utility functions and data structures for mathamatecal operations.
 */
namespace DoobMathUtils {
	/**
	 * Generates a random vector perpendicular to the given base vector.
	 * @param BaseVector The input vector to generate a perpendicular vector.
	 * @return A normalized vector perpendicular to the BaseVector.
	 */
	FVector GenerateRandomPerpendicularVector(const FVector& BaseVector);

	/**
	 * Calculates an ease-out increment value based on the current step.
	 * @param TotalSteps The total number of steps in the process.
	 * @param CurrentStep The current step of the process.
	 * @return The ease-out increment value.
	 */
	float EaseOutIncrement(int32 TotalSteps, int32 CurrentStep);

	/**
	 * Calculates an ease-in increment value based on the current step.
	 * @param TotalSteps The total number of steps in the process.
	 * @param CurrentStep The current step of the process.
	 * @return The ease-in increment value.
	 */
	float EaseInIncrement(int32 TotalSteps, int32 CurrentStep);

	/**
	 * Adjusts a floating-point value using an easing curve (ease-in or ease-out).
	 * @param StartRadius The starting value of the range.
	 * @param EndRadius The ending value of the range.
	 * @param CurrentRadius The current radius being adjusted.
	 * @param TotalSteps The total number of steps in the process.
	 * @param CurrentStep The current step of the process.
	 * @param CurveType The type of easing curve (0 for ease-out, 1 for ease-in).
	 * @return The adjusted floating-point value.
	 */
	float FloatChangeWithCurve(float StartRadius, float EndRadius, float CurrentRadius, int32 TotalSteps, int32 CurrentStep, int32 CurveType);

	/**
	 * Interpolates between two bounds, dividing the range into equal segments.
	 * @param StartBound The lower bound of the range.
	 * @param EndBound The upper bound of the range.
	 * @param NumSegments The number of segments to divide the range into.
	 * @param SegmentIndex The index of the segment for which to calculate the bounds.
	 * @param StartValue Output parameter for the start of the segment.
	 * @param EndValue Output parameter for the end of the segment.
	 */
	void InterpolateBetweenBounds(float StartBound, float EndBound, int32 NumSegments, int32 SegmentIndex, float& StartValue, float& EndValue);
}