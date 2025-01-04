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

	void GenerateIntersectionRing(
		const FTubeData& MainTube,
		const FTubeData& LateralTube,
		FIntersectionRingData& OutRingData,
		float Precision
	) {
		TArray<FVector> CombinedVertices;

		GenerateHalfIntersectionRing(MainTube, LateralTube, OutRingData.MainTubeVertices);
		GenerateHalfIntersectionRing(LateralTube, MainTube, OutRingData.LateralTubeVertices);

		CombinedVertices = OutRingData.MainTubeVertices;
		CombinedVertices.Append(OutRingData.LateralTubeVertices);

		OutRingData.CombinedVertices = OrderRingVertices(CombinedVertices);

		FVector FirstVertex = OutRingData.CombinedVertices[0];
		OutRingData.CombinedVertices.Add(FirstVertex);
	}

	void GenerateHalfIntersectionRing(
		const FTubeData& TubeA,
		const FTubeData& TubeB,
		TArray<FVector>& OutRingVertices,
		float Precision
	) {
		for (int32 LineIndexA = 0; LineIndexA < TubeA.Rings.Num() - 1; ++LineIndexA) {
			const FRingData& CurrentRingA = TubeA.Rings[LineIndexA];
			const FRingData& NextRingA = TubeA.Rings[LineIndexA + 1];
			int32 NumVerticesA = CurrentRingA.Vertices.Num();

			UE_LOG(LogTemp, Log, TEXT("NumVerticesA: %d"), NumVerticesA);

			for (int32 VertexIndexA = 0; VertexIndexA < NumVerticesA; ++VertexIndexA) {
				FVector LineStartA = CurrentRingA.Vertices[VertexIndexA];
				FVector LineEndA = NextRingA.Vertices[VertexIndexA];
				
				// loop through TubeB rectangles
				for (int32 RingIndexB = 0; RingIndexB < TubeB.Rings.Num() - 1; ++RingIndexB) {
					const FRingData& CurrentRingB = TubeB.Rings[RingIndexB];
					const FRingData& NextRingB = TubeB.Rings[RingIndexB + 1];
					int32 NumVerticesB = CurrentRingB.Vertices.Num();

					for (int32 VertexIndexB = 0; VertexIndexB < NumVerticesB; ++VertexIndexB) {
						FVector V0 = CurrentRingB.Vertices[VertexIndexB];
						FVector V1 = CurrentRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];
						FVector V2 = NextRingB.Vertices[VertexIndexB];
						FVector V3 = NextRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];

						// check for intersection with triangles
						FVector IntersectionPoint;
						if (LineSegmentIntersectsTriangle(LineStartA, LineEndA, V0, V1, V2, IntersectionPoint) ||
							LineSegmentIntersectsTriangle(LineStartA, LineEndA, V1, V2, V3, IntersectionPoint)) {
							OutRingVertices.Add(IntersectionPoint);
						}
					}
				}
			}
		}
	}

	bool LineSegmentIntersectsTriangle(
		const FVector& LineStart,
		const FVector& LineEnd,
		const FVector& V0,
		const FVector& V1,
		const FVector& V2,
		FVector& OutIntersectionPoint
	) {
		FVector LineDir = LineEnd - LineStart;
		FVector Edge1 = V1 - V0;
		FVector Edge2 = V2 - V0;
		FVector PVec = FVector::CrossProduct(LineDir, Edge2);
		float Det = FVector::DotProduct(Edge1, PVec);

		// parallel check
		if (FMath::Abs(Det) < KINDA_SMALL_NUMBER) {
			return false;
		}

		float InvDet = 1.0f / Det;
		FVector TVec = LineStart - V0;
		float U = FVector::DotProduct(TVec, PVec) * InvDet;

		if (U < 0.0f || U > 1.0f) {
			return false;
		}

		FVector QVec = FVector::CrossProduct(TVec, Edge1);
		float V = FVector::DotProduct(LineDir, QVec) * InvDet;

		if (V < 0.0f || (U + V) > 1.0f) {
			return false;
		}

		float T = FVector::DotProduct(Edge2, QVec) * InvDet;
		
		if (T < 0.0f || T > 1.0f) {
			return false;
		}

		// compute intersection point
		OutIntersectionPoint = LineStart + T * LineDir;
		return true;
	}

	FVector ComputeCentroid(const TArray<FVector>& Vertices) {
		FVector Centroid = FVector::ZeroVector;
		for (const FVector& Vertex : Vertices) {
			Centroid += Vertex;
		}
		return Centroid / Vertices.Num();
	}

	TArray<FVector> OrderRingVertices(const TArray<FVector>& InputVertices) {
		FVector Centroid = ComputeCentroid(InputVertices);

		// array to store vertices with angles
		TArray<TPair<float, FVector>> VerticesWithAngles;

		// define the plane normal (approx, using 1st 2 vectors)
		//FVector PlaneNormal = FVector::CrossProduct(InputVertices[1] - Centroid, InputVertices[0] - Centroid).GetSafeNormal();

		// Compute a better plane normal based on averaging cross products of centroid and vertices
		FVector AverageNormal = FVector::ZeroVector;
		for (int32 i = 0; i < InputVertices.Num(); i++) {
			FVector NextVertex = InputVertices[(i + 1) % InputVertices.Num()];
			FVector Edge1 = InputVertices[i] - Centroid;
			FVector Edge2 = NextVertex - Centroid;
			AverageNormal += FVector::CrossProduct(Edge1, Edge2);
		}
		AverageNormal = AverageNormal.GetSafeNormal(); // Normalize the average normal

		// project vertices onto a plane
		for (const FVector& Vertex : InputVertices) {
			FVector Offset = Vertex - Centroid;
			FVector Projected = Offset - FVector::DotProduct(Offset, AverageNormal) * AverageNormal;

			// compute angle in projected space
			float Angle = FMath::Atan2(Projected.Y, Projected.X);
			VerticesWithAngles.Add(TPair<float, FVector>(Angle, Vertex));
		}

		// sort vertices by angle
		VerticesWithAngles.Sort([](const TPair<float, FVector>& A, const TPair<float, FVector>& B) {
			return A.Key < B.Key; // sort by angle
		});

		// extract sorted vertices
		TArray<FVector> SortedVertices;
		for (const TPair<float, FVector>& Pair : VerticesWithAngles) {
			SortedVertices.Add(Pair.Value);
		}

		return SortedVertices;
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
		: Normal(), D() {
	}

	FPlaneEquation::FPlaneEquation(const FVector& InNormal, float InD)
		: Normal(InNormal), D(InD) {
	}

	FPlaneEquation CalculatePlaneFromRings(const FRingData& RingA, const FRingData& RingB) {
		// calc teh direction vector between the centers of the 2 rings
		FVector DirectionBetweenRings = (RingB.Center - RingA.Center).GetSafeNormal();

		// Compute the slope adjustment vector based on the radius difference
		float RadiusDifference = RingB.Radius - RingA.Radius;
		FVector SlopeAdjustment = RadiusDifference * DirectionBetweenRings;

		// adjust the normal of RingA to account for the slope
		FVector AdjustedNormalA = (RingA.Normal + SlopeAdjustment).GetSafeNormal();

		// compute the plane normal by taking the cross product of DirectionBetweenRings and RingA's normal
		FVector PlaneNormal = FVector::CrossProduct(DirectionBetweenRings, AdjustedNormalA).GetSafeNormal();

		// calculate the D value for the plane equation, adjust for shifted center
		FVector AdjustedCenter = RingA.Center + SlopeAdjustment;
		float D = FVector::DotProduct(PlaneNormal, AdjustedCenter);

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

		float Numerator = -(FVector::DotProduct(Plane.Normal, RayOrigin) + Plane.D);

		// check if the ray is parallel to the plane
		if (FMath::IsNearlyZero(NDotD)) {
			float DistanceFromPlane = -Numerator;

			// check if the ray origin lies on the plane
			if (FMath::IsNearlyZero(DistanceFromPlane)) {
				OutIntersectionPoint = RayOrigin; // orign is an intersection point
				return true;
			}

			return false; // no intersection, ray parallel to the plane
		}

		// calc t, the parameter of the ray at the intersection point
		float t = Numerator / NDotD;

		// if t < 0, intersection is behind the ray origin
		if (t < 0.0f) {
			return false; // no intersect the ray does not reach the plane
		}

		// compute the intersection point using the ray equation
		OutIntersectionPoint = RayOrigin + t * RayDirection;

		return true; // successfully calculated an intersection point 
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
		TArray<TPair<int32, int32>>& OutPlaneIndicesA,
		TArray<TPair<int32, int32>>& OutPlaneIndicesB,
		TArray<int32>& OutLineIndicesA,
		TArray<int32>& OutLineIndicesB
	) {
		// process TubeA planes and TubeB lines
		for (int32 RingIndexA = 0; RingIndexA < TubeA.Rings.Num(); ++RingIndexA) {
			const FRingData& RingA = TubeA.Rings[RingIndexA];
			float RadiusA = RingA.Radius;
			int32 NumSidesA = RingA.Vertices.Num();

			for (int32 SideIndexA = 0; SideIndexA < NumSidesA; ++SideIndexA) {
				FVector SideStartA = RingA.Vertices[SideIndexA];
				FVector SideEndA = RingA.Vertices[(SideIndexA + 1) % NumSidesA];

				// loop through TubeB lines
				for (int32 LineIndexB = 0; LineIndexB < TubeB.Rings.Num() - 1; ++LineIndexB) {
					const FRingData& CurrentRingB = TubeB.Rings[LineIndexB];
					const FRingData& NextRingB = TubeB.Rings[LineIndexB + 1];
					int32 NumVerticesB = CurrentRingB.Vertices.Num();

					for (int32 VertexIndexB = 0; VertexIndexB < NumVerticesB; ++VertexIndexB) {
						if (!CurrentRingB.Vertices.IsValidIndex(VertexIndexB) ||
							!NextRingB.Vertices.IsValidIndex(VertexIndexB)) {
							continue;
						}

						FVector LineStartB = CurrentRingB.Vertices[VertexIndexB];
						FVector LineEndB = NextRingB.Vertices[VertexIndexB];

						FVector PlaneNormal = (SideEndA - SideStartA).GetSafeNormal();
						float PlaneD = FVector::DotProduct(SideStartA, PlaneNormal);

						FVector LineDir = LineEndB - LineStartB;
						float Denom = FVector::DotProduct(LineDir, PlaneNormal);

						// Check for parallel lines
						if (FMath::Abs(Denom) < KINDA_SMALL_NUMBER) {
							continue;
						}

						// Calculate T (intersection point parameter)
						float T = (PlaneD - FVector::DotProduct(LineStartB, PlaneNormal)) / Denom;

						// Ensure T lies within [0, 1]
						if (T < 0.0f || T > 1.0f) {
							continue;
						}

						FVector IntersectionPoint = LineStartB + T * LineDir;
						float DistanceToCenter = FVector::Dist(IntersectionPoint, SideStartA);

						// Check if intersection point is within the ring's radius
						if (DistanceToCenter <= RadiusA) {
							TPair<int32, int32> PlaneIndexA(RingIndexA, SideIndexA);
							if (!OutPlaneIndicesA.Contains(PlaneIndexA)) {
								OutPlaneIndicesA.Add(PlaneIndexA);
							}
							if (!OutLineIndicesB.Contains(LineIndexB)) {
								OutLineIndicesB.Add(LineIndexB);
							}
						}
					}
				}
			}
		}

		// process TubeB planes and TubeA lines
		for (int32 RingIndexB = 0; RingIndexB < TubeB.Rings.Num(); ++RingIndexB) {
			const FRingData& RingB = TubeB.Rings[RingIndexB];
			float RadiusB = RingB.Radius;
			int32 NumSidesB = RingB.Vertices.Num();

			for (int32 SideIndexB = 0; SideIndexB < NumSidesB; ++SideIndexB) {
				FVector SideStartB = RingB.Vertices[SideIndexB];
				FVector SideEndB = RingB.Vertices[(SideIndexB + 1) % NumSidesB];

				// loop through TubeA lines
				for (int32 LineIndexA = 0; LineIndexA < TubeA.Rings.Num(); ++LineIndexA) {
					const FRingData& RingA = TubeA.Rings[LineIndexA];
					int32 NumVerticesA = RingA.Vertices.Num();

					for (int32 VertexIndexA = 0; VertexIndexA < NumVerticesA; ++VertexIndexA) {
						FVector LineStartA = RingA.Vertices[VertexIndexA];
						FVector LineEndA = RingA.Vertices[(VertexIndexA + 1) % NumVerticesA];

						float DistanceStart = FVector::PointPlaneDist(LineStartA, SideStartB, SideEndB - SideStartB);
						float DistanceEnd = FVector::PointPlaneDist(LineEndA, SideStartB, SideEndB - SideStartB);


						if (FMath::Abs(DistanceStart) <= RadiusB || FMath::Abs(DistanceStart) <= RadiusB) {
							TPair<int32, int32> PlaneIndexB(RingIndexB, SideIndexB);
							if (!OutPlaneIndicesB.Contains(PlaneIndexB)) {
								OutPlaneIndicesB.Add(PlaneIndexB);
							}
							if (!OutLineIndicesA.Contains(LineIndexA)) {
								OutLineIndicesA.Add(LineIndexA);
							}
						}
					}
				}
			}
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
}