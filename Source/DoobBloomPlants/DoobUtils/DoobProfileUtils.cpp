#include "DoobProfileUtils.h"

namespace DoobProfileUtils {

	F2DProfile::F2DProfile() : Points(), bIsClosed(true) {

	}

	F2DProfile::F2DProfile(TArray<FVector2D> InPoints, bool InClosed)
		: Points(InPoints), bIsClosed(InClosed) {
	}

	// scale the profile by a given factor
	void F2DProfile::Scale(float Factor) {
		for (FVector2D& Point : Points) {
			Point *= Factor;
		}
	}

	// translate the profile by given offset
	void F2DProfile::Translate(const FVector2D& Offset) {
		for (FVector2D& Point : Points) {
			Point += Offset;
		}
	}

	// rotate the profile arount the origin by a given angle in radians
	void F2DProfile::Rotate(float Angle) {
		float CosAngle = FMath::Cos(Angle);
		float SinAngle = FMath::Sin(Angle);

		for (FVector2D& Point : Points) {
			float X = Point.X * CosAngle - Point.Y * SinAngle;
			float Y = Point.X * SinAngle + Point.Y * CosAngle;
			Point = FVector2D(X, Y);
		}
	}


	F2DProfile GenerateLinearProfile(int32 NumSegments, float StartHeight, float EndHeight) {
		TArray<FVector2D> Points;

		float Slope = (EndHeight - StartHeight) / NumSegments;

		for (int32 i = 0; i < NumSegments; ++i) {
			float t = static_cast<float>(i) / (NumSegments - 1);
			float X = t;
			float Y = FMath::Lerp(StartHeight, EndHeight, t);
			Points.Add(FVector2D(X, Y));
		}

		return F2DProfile(Points, true);
	}

	F2DProfile GenerateEggShapedCylinderProfile(int32 NumSegments, float StartHeight, float MidHeight, float EndHeight, float BulgePosition) {
		TArray<FVector2D> Points;

		float step = 0.02f;
		float currentStep = 0.0f;

		for (int32 i = 0; i < NumSegments; ++i) {
			float t = static_cast<float>(i) / (NumSegments - 1);
			float bulgeTAdjust = (NumSegments - 1) * BulgePosition;
			float bulgeT = static_cast<float>(i) / bulgeTAdjust;
			float startT = FMath::Pow(bulgeT, 0.4f);
			float startTStep = FMath::Clamp(startT + currentStep, 0, 1);
			float endT = FMath::Clamp((static_cast<float>(i) - bulgeTAdjust) / ((NumSegments - 1) * (1.0f - BulgePosition)), 0, 1);

			// calc the y using eased params
			float Y = FMath::Lerp(
				FMath::Lerp(StartHeight, MidHeight, startTStep),
				FMath::Lerp(MidHeight, EndHeight, endT * endT * endT),
				FMath::Clamp(bulgeT, 0, 1)
			);

			currentStep += step;
			step = FMath::Pow(step, 1.3f);

			// pass t to x
			float X = t;

			Points.Add(FVector2D(X, Y));
		}

		return F2DProfile(Points, true);
	}

	FVector2D EvaluateProfile(const F2DProfile& Profile, float t) {
		if (Profile.Points.Num() == 0) {
			// handle empty profile
			return FVector2D::ZeroVector;
		}

		int32 NumPoints = Profile.Points.Num();
		t = FMath::Clamp(t, 0.0f, 1.0f); // ensure t in bounds

		float ScaledT = t * (Profile.bIsClosed ? NumPoints : (NumPoints - 1));
		int32 IndexA = FMath::FloorToInt(ScaledT); // Base index
		int32 IndexB = (IndexA + 1) % NumPoints; // next index, wraps if closed

		if (!Profile.bIsClosed && IndexB >= NumPoints) {
			// if profile is open and at the last point
			return Profile.Points.Last();
		}

		float LocalT = ScaledT - IndexA; // fraction between index a and index b
		return FMath::Lerp(Profile.Points[IndexA], Profile.Points[IndexB], LocalT);
	}

	float SumYValuesInProfile(const F2DProfile& Profile) {
		float SumY = 0.0f;

		for (const FVector2D& Point : Profile.Points) {
			SumY += Point.Y;
		}

		return SumY;
	}

	float MaxYValueInProfile(const F2DProfile& Profile) {
		float Max = 0.0f;

		for (const FVector2D& Point : Profile.Points) {
			if (Point.Y > Max) {
				Max = Point.Y;
			}
		}

		return Max;
	}
}