#include "DoobGeometryUtils.h"

#include "Math/Vector.h"
#include "Math/Transform.h"
#include "Containers/Array.h"

//#include "DoobProfileUtils.h"

namespace DoobGeometryUtils {
	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, TArray<FVector>& RingVertices) {
		const float AngleStep = 360.0f / NumSides;
		FVector Right = FVector::CrossProduct(Direction, UpVector).GetSafeNormal();
		FVector Forward = FVector::CrossProduct(Right, Direction).GetSafeNormal();

		for (int i = 0; i <= NumSides; ++i) {
			float Angle = FMath::DegreesToRadians(i * AngleStep);
			FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;
			RingVertices.Add(Center + Offset);
		}
	}

	void ConnectRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		for (int32 i = 0; i < RingA.Num(); ++i) {
			int32 CurrentA = BaseIndex + i;
			int32 NextA = BaseIndex + (i + 1) % RingA.Num();

			int32 CurrentB = BaseIndex + RingA.Num() + i;
			int32 NextB = BaseIndex + RingA.Num() + ((i + 1) % RingB.Num());

			// Triangle 1
			Triangles.Add(CurrentA);
			Triangles.Add(NextA);
			Triangles.Add(CurrentB);

			// Triangle 2
			Triangles.Add(NextA);
			Triangles.Add(NextB);
			Triangles.Add(CurrentB);
		}
		Vertices.Append(RingA);
		Vertices.Append(RingB);
		BaseIndex += RingA.Num() + RingB.Num();
	}

	void ConnectRingArray(const TArray<TArray<FVector>>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		if (Rings.Num() < 2) return;

		// iterate through the array of rings and connect consecutive rings
		for (int32 i = 0; i < Rings.Num() - 1; ++i) {
			const TArray<FVector>& CurrentRing = Rings[i];
			const TArray<FVector>& NextRing = Rings[i + 1];

			ConnectRings(CurrentRing, NextRing, Vertices, Triangles, BaseIndex);
		}

	}

	TArray<FVector> GenerateIntersectionCurve(
		const DoobProfileUtils::F2DProfile& MainProfile,
		const FTransform& MainTransform,
		const DoobProfileUtils::F2DProfile& LateralProfile,
		const FTransform& LateralTransform,
		int32 MainSegments,
		int32 LateralSegments,
		float Threshold
	) {
		TArray<FVector> IntersectionCurve;

		for (int32 i = 0; i < MainSegments; ++i) {
			// parametric value for the main cylinder
			float tMain = static_cast<float>(i) / MainSegments;

			// evaluate the main profile
			FVector2D Main2DPoint = DoobProfileUtils::EvaluateProfile(MainProfile, tMain);
			FVector MainPoint = MainTransform.TransformPosition(FVector(0.0f, Main2DPoint.X, Main2DPoint.Y));

			for (int32 j = 0; j < LateralSegments; ++j) {
				// parametric value for the lateral cylinder
				float tLateral = static_cast<float>(j) / LateralSegments;

				// Evaluate the lateral profile
				FVector2D Lateral2DPoint = DoobProfileUtils::EvaluateProfile(LateralProfile, tLateral);
				FVector LateralPoint = LateralTransform.TransformPosition(FVector(Lateral2DPoint.Y, 0.0f, Lateral2DPoint.X));

				// compute the distance between the points
				float Distance = FVector::Dist(MainPoint, LateralPoint);

				// if the points are close enough, add to the intersection curve
				if (Distance < Threshold) {
					IntersectionCurve.Add((MainPoint + LateralPoint) / 2); // avg the mid pt
				}
			}
		}

		return IntersectionCurve;
	}

	void ConstructTubeFromProfile(
		const DoobProfileUtils::F2DProfile& Profile,
		const FVector& StartPosition,
		const FVector& EndPosition,
		const FVector& Direction,
		const FVector& UpVector,
		int32 NumSegments,
		int32 NumSides,
		float Height,
		TArray<TArray<FVector>>& OutRingVerticesArrays,
		bool ApplyBothAxis
	) {
		OutRingVerticesArrays.Empty();

		FVector CurrentPosition = StartPosition;

		// Normalize the input direction to avoid errors
		FVector NormalizedDirection = Direction.GetSafeNormal();

		float SegmentLengthCurveTotal = 0.0f;
		float SegmentLengthMax = 0.0f;
		bool MissingSegmentPassed = false;

		if (ApplyBothAxis) {
			SegmentLengthCurveTotal = DoobProfileUtils::SumYValuesInProfile(Profile);
			SegmentLengthMax = DoobProfileUtils::MaxYValueInProfile(Profile);
		}

		// Iterate over height segments to generate rings
		for (int32 i = 0; i < NumSegments; ++i) {
			float t = static_cast<float>(i) / NumSegments;

			float segmentLength = 0.0f;

			if (ApplyBothAxis) {
				float segPercent = Profile.Points[i].Y / SegmentLengthCurveTotal;

				if (i == 1) {
					segPercent = Profile.Points[0].Y / SegmentLengthCurveTotal;
				}
				else if (Profile.Points[i].Y >= SegmentLengthMax - 1.0f && !MissingSegmentPassed) {
					segPercent = (Profile.Points[i].Y / SegmentLengthCurveTotal) + (Profile.Points[1].Y / SegmentLengthCurveTotal);
					MissingSegmentPassed = true;
				}

				segmentLength = segPercent * Height;

				CurrentPosition = CurrentPosition + segmentLength * NormalizedDirection;

				UE_LOG(LogTemp, Log, TEXT("iteration: %d, segPercent: %f, segmentLength: %f"), i, segPercent, segmentLength);
			}
			else {
				segmentLength = t * Height;
				CurrentPosition = StartPosition + segmentLength * NormalizedDirection;
			}

			// Calculate the position along the height
			//CurrentPosition = StartPosition + segmentLength * NormalizedDirection;

			// Generate a ring at this position
			TArray<FVector> RingVertices;
			DoobGeometryUtils::GenerateRing(CurrentPosition, NormalizedDirection, UpVector, Profile.Points[i].Y, NumSides, RingVertices);

			// Store the ring vertices in the output array
			OutRingVerticesArrays.Add(RingVertices);
		}		
	}
}