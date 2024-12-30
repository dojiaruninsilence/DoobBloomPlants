#pragma once

#include "CoreMinimal.h"
#include "Math/Vector2D.h"
#include "Containers/Array.h"

/**
 * @namespace DoobProfileUtils
 * Utility functions and classes for generating and manipulating 2D profiles.
 */
namespace DoobProfileUtils {

	/**
	 * @struct F2DProfile
	 * Represents a 2D profile defined by a set of points.
	 */
	struct F2DProfile {
		TArray<FVector2D> Points; ///< Points defining the profile.
		bool bIsClosed; ///< Whether the profile is closed (looped).

		/**
		 * Default constructor.
		 */
		F2DProfile();

		/**
		 * Constructs an F2DProfile with specified points.
		 * @param InPoints Points defining the profile.
		 * @param InClosed Whether the profile is closed (default: true).
		 */
		F2DProfile(TArray<FVector2D> InPoints, bool InClosed = true);

		/**
		 * Scales the profile by a given factor.
		 * @param Factor Scaling factor.
		 */
		void Scale(float Factor);

		/**
		 * Translates the profile by a given offset.
		 * @param Offset Translation vector.
		 */
		void Translate(const FVector2D& Offset);

		/**
		 * Rotates the profile around the origin by a given angle in radians.
		 * @param Angle Angle in radians.
		 */
		void Rotate(float Angle);
	};

	// ----------------------------------------------------------- Generate Profiles ----------------------------------------------------------------------- // 

	/**
	 * Generates a linear profile.
	 * @param NumSegments Number of segments in the profile.
	 * @param StartHeight Starting height of the profile.
	 * @param EndHeight Ending height of the profile.
	 * @return A linear 2D profile.
	 */
	F2DProfile GenerateLinearProfile(int32 NumSegments, float StartHeight, float EndHeight);

	/**
	 * Generates an egg-shaped cylindrical profile with a bulge.
	 * @param NumSegments Number of segments in the profile.
	 * @param StartHeight Starting height of the profile.
	 * @param MidHeight Maximum height at the bulge.
	 * @param EndHeight Ending height of the profile.
	 * @param BulgePosition Position of the bulge (0 to 1).
	 * @return An egg-shaped cylindrical profile.
	 */
	F2DProfile GenerateEggShapedCylinderProfile(int32 NumSegments, float StartHeight, float MidHeight, float EndHeight, float BulgePosition);

	// ----------------------------------------------------------- Retrieve Data from Profiles ----------------------------------------------------------------------- //

	/**
	 * Evaluates a profile at a given parameter t.
	 * @param Profile The profile to evaluate.
	 * @param t A parameter (0 to 1) indicating the position along the profile.
	 * @return The evaluated point on the profile.
	 */
	FVector2D EvaluateProfile(const F2DProfile& Profile, float t);

	/**
	 * Calculates the sum of Y-values in a profile.
	 * @param Profile The profile to process.
	 * @return The sum of all Y-values.
	 */
	float SumYValuesInProfile(const F2DProfile& Profile);

	/**
	 * Finds the maximum Y-value in a profile.
	 * @param Profile The profile to process.
	 * @return The maximum Y-value.
	 */
	float MaxYValueInProfile(const F2DProfile& Profile);
}