#include "DoobGeometryUtils.h"

#include "Math/Vector.h"
#include "Math/Transform.h"
#include "Containers/Array.h"

//#include "DoobProfileUtils.h"

namespace DoobGeometryUtils {
	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, FRingData& RingData) {

		// <------------------------------------- Need to figure out if this is necassary and if so why it ruins geometry ----------------------------------------------> //
		// Ensure the UpVector is not parallel to the Direction to avoid degenerate cross-product results
		//FVector SafeUpVector = FMath::IsNearlyZero(Direction | UpVector) ? FVector::UpVector : UpVector;
		// <---------------- Temp Remove all this logic if above cant be fixed ------------------------> // 
		FVector SafeUpVector = UpVector;

		const float AngleStep = 360.0f / NumSides;
		FVector Right = FVector::CrossProduct(Direction, SafeUpVector).GetSafeNormal();
		FVector Forward = FVector::CrossProduct(Right, Direction).GetSafeNormal();

		FVector FirstVertex;

		for (int i = 0; i <= NumSides; ++i) {
			float Angle = FMath::DegreesToRadians(i * AngleStep);
			FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;

			FVector Vertex = Center + Offset;

			if (i == 0) {
				FirstVertex = Vertex; // Store the first vertex
			}

			RingData.Vertices.Add(Vertex);
		}

		// close the ring by adding the first vertex again
		RingData.Vertices.Add(FirstVertex);

		// calc the normal of the ring using right and forward vectors
		RingData.Normal = FVector::CrossProduct(Right, Forward).GetSafeNormal();

		RingData.Center = Center;
		RingData.Direction = Direction;
		RingData.Radius = Radius;
		RingData.UpVector = UpVector;
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

	void ConnectRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		if (Rings.Num() < 2) return;

		// iterate through the array of rings and connect consecutive rings
		for (int32 i = 0; i < Rings.Num() - 1; ++i) {
			const TArray<FVector>& CurrentRing = Rings[i].Vertices;
			const TArray<FVector>& NextRing = Rings[i + 1].Vertices;

			ConnectRings(CurrentRing, NextRing, Vertices, Triangles, BaseIndex);
		}

	}

	void ConstructTubeFromProfile(
		const DoobProfileUtils::F2DProfile& Profile,
		const FVector& StartPosition,
		const FVector& Direction,
		const FVector& UpVector,
		int32 NumSegments,
		int32 NumSides,
		float Length,
		FTubeData& OutTube,
		bool ApplyBothAxis
	) {
		OutTube.Rings.Empty();

		FVector CurrentPosition = StartPosition;

		// Normalize the input direction to avoid errors
		FVector NormalizedDirection = Direction.GetSafeNormal();

		float SegmentLengthCurveTotal = 0.0f;
		float SegmentLengthMax = 0.0f;
		bool MissingSegmentPassed = false;

		// if true apply the curve to both axis
		if (ApplyBothAxis) {
			SegmentLengthCurveTotal = DoobProfileUtils::SumYValuesInProfile(Profile);
			SegmentLengthMax = DoobProfileUtils::MaxYValueInProfile(Profile);
		}

		// Iterate over height segments to generate rings
		for (int32 i = 0; i < NumSegments; ++i) {
			float t = static_cast<float>(i) / NumSegments;

			float segmentLength = 0.0f;

			// Generate a ring at this position
			DoobGeometryUtils::FRingData RingData;
			DoobGeometryUtils::GenerateRing(CurrentPosition, NormalizedDirection, UpVector, Profile.Points[i].Y, NumSides, RingData);

			// Store the ring vertices in the output array
			OutTube.Rings.Add(RingData);

			// if true apply the curve to both axis
			if (ApplyBothAxis) {
				float segPercent = Profile.Points[i].Y / SegmentLengthCurveTotal;

				if (i == 1) {
					segPercent = Profile.Points[0].Y / SegmentLengthCurveTotal;
				}
				else if (Profile.Points[i].Y >= SegmentLengthMax - 1.0f && !MissingSegmentPassed) {
					segPercent = (Profile.Points[i].Y / SegmentLengthCurveTotal) + (Profile.Points[1].Y / SegmentLengthCurveTotal);
					MissingSegmentPassed = true;
				}

				segmentLength = segPercent * Length;

				CurrentPosition = CurrentPosition + segmentLength * NormalizedDirection;

				UE_LOG(LogTemp, Log, TEXT("iteration: %d, segPercent: %f, segmentLength: %f"), i, segPercent, segmentLength);
			}
			else {
				segmentLength = t * Length;
				CurrentPosition = StartPosition + segmentLength * NormalizedDirection;
			}
		}

		OutTube.Direction = Direction;
		OutTube.EndPosition = CurrentPosition;
		OutTube.Length = Length;
		OutTube.NumSegments = NumSegments;
		OutTube.NumSides = NumSides;
		OutTube.Profile = Profile;
		OutTube.StartPosition = StartPosition;
	}

	FPlaneEquation::FPlaneEquation()
		: Normal(), D() {}

	FPlaneEquation::FPlaneEquation(const FVector& InNormal, float InD)
		: Normal(InNormal), D(InD) { }

	FPlaneEquation CalculatePlaneFromRings(const FRingData& RingA, const FRingData& RingB) {
		// calc teh direction vector between the centers of the 2 rings
		FVector DirectionBetweenRings = (RingB.Center - RingA.Center).GetSafeNormal();

		// use the normal of RingA directly
		FVector NormalA = RingA.Normal;

		// compute the plane normal by taking the cross product of DirectionBetweenRings and RingA's normal
		FVector PlaneNormal = FVector::CrossProduct(DirectionBetweenRings, NormalA).GetSafeNormal();

		// calculate the D value for the plane equation
		float D = FVector::DotProduct(PlaneNormal, RingA.Center);

		return FPlaneEquation(PlaneNormal, D);
	}

	bool CalculateIntersectionWithPlane(
		const FPlaneEquation& Plane,
		const FVector& RayOrigin,
		const FVector& RayDirection,
		FVector& OutIntersectionPoint
	) {
		// calc the dot product of the plane's normal and the ray direction
		float NDotD = FVector::DotProduct(Plane.Normal, RayDirection);

		// check if the ray is parallel to the plane
		if (FMath::IsNearlyZero(NDotD)) {
			return false; // no intersection, ray parallel to the plane
		}

		// calc t, the parameter of the ray at the intersection point
		float t = -(FVector::DotProduct(Plane.Normal, RayOrigin) + Plane.D) / NDotD;

		// if t < 0, intersection is behind the ray origin
		if (t < 0.0f) {
			return false; // no intersect the ray does not reach the plane
		}

		// compute the intersection point using the ray equation
		OutIntersectionPoint = RayOrigin + t * RayDirection;

		return true; // successfull calculated
	}

	bool PointInsideCircle(const FVector& Point, const FVector& Center, float Radius) {
		// calc the squared distance between the point and the circle's center
		float DistanceSquared = FVector::DistSquared(Point, Center);

		// compare squared distance to squared radius for efficiency
		return DistanceSquared <= FMath::Square(Radius);
	}

	void FilterPlanesAndLines(
		const FTubeData& TubeA,
		const FTubeData& TubeB,
		TArray<int32>& OutPlaneIndicesA,
		TArray<int32>& OutPlaneIndicesB,
		TArray<int32>& OutLineIndicesA,
		TArray<int32>& OutLineIndicesB
	) {
		UE_LOG(LogTemp, Log, TEXT("Starting FilterPlanesAndLines"));
		UE_LOG(LogTemp, Log, TEXT("TubeA has %d rings, TubeB has %d rings"), TubeA.Rings.Num(), TubeB.Rings.Num());

		for (int32 i = 0; i < TubeA.Rings.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("Processing TubeA Ring %d: Center %s, Radius %f"), i, *TubeA.Rings[i].Center.ToString(), TubeA.Rings[i].Radius);

			for (int32 j = 0; j < TubeB.Rings.Num(); ++j) {
				UE_LOG(LogTemp, Log, TEXT("Processing TubeB Ring %d: Center %s, Radius %f"), j, *TubeB.Rings[j].Center.ToString(), TubeB.Rings[j].Radius);

				float Distance = FVector::Dist(TubeA.Rings[i].Center, TubeB.Rings[j].Center);
				float CombinedRadius = TubeA.Rings[i].Radius + TubeB.Rings[j].Radius;

				UE_LOG(LogTemp, Log, TEXT("Distance between TubeA Ring %d and TubeB Ring %d: %f"), i, j, Distance);
				UE_LOG(LogTemp, Log, TEXT("Combined radius: %f"), CombinedRadius);

				if (Distance <= CombinedRadius) {
					UE_LOG(LogTemp, Log, TEXT("Rings intersect: TubeA Ring %d and TubeB Ring %d"), i, j);

					// add rings that are close enough to possibly intersect
					OutPlaneIndicesA.Add(i);
					OutPlaneIndicesB.Add(j);

					// add connected lines for TubeA
					OutLineIndicesA.Add(i);
					OutLineIndicesA.Add((i + 1) % TubeA.Rings.Num());

					// add connected lines for TubeB
					OutLineIndicesB.Add(j);
					OutLineIndicesB.Add((j + 1) % TubeB.Rings.Num());
				}
				else {
					UE_LOG(LogTemp, Log, TEXT("Rings do not intersect: TubeA Ring %d and TubeB Ring %d"), i, j);
				}
			}
		}

		UE_LOG(LogTemp, Log, TEXT("Filtering complete"));
		UE_LOG(LogTemp, Log, TEXT("OutPlaneIndicesA: %d, OutPlaneIndicesB: %d, OutLineIndicesA: %d, OutLineIndicesB: %d"),
			OutPlaneIndicesA.Num(), OutPlaneIndicesB.Num(), OutLineIndicesA.Num(), OutLineIndicesB.Num());

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
}