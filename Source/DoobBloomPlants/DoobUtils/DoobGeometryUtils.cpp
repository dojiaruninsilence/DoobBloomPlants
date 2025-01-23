#include "DoobGeometryUtils.h"

#include "Math/Vector.h"
#include "Math/Transform.h"
#include "Containers/Array.h"
#include "DrawDebugHelpers.h"

//#include "DoobProfileUtils.h"
#include "DoobContainerUtils.h"

namespace DoobGeometryUtils {
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                                 2. Ring Operations                                                                   //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

	void GenerateRing(FVector Center, FVector Direction, FVector UpVector, float Radius, int32 NumSides, FRingData& RingData) {
		// Ensure the UpVector is valid and not parallel to the Direction to avoid degenerate results
		// For now, we assume the provided UpVector is valid, but this logic could be reintroduced if needed.
		FVector SafeUpVector = UpVector;

		// Calculate the angle step for creating vertices around the ring
		const float AngleStep = 360.0f / NumSides;

		// Calculate the initial basis vectors for constructing the ring
		FVector Right = FVector::CrossProduct(Direction, SafeUpVector).GetSafeNormal();
		FVector Forward = FVector::CrossProduct(Right, Direction).GetSafeNormal();

		FVector FirstVertex; // Variable to store the first vertex for potential use later

		// Generate vertices around the ring
		for (int i = 0; i <= NumSides; ++i) {
			// Calculate the angle for the current vertex
			float Angle = FMath::DegreesToRadians(i * AngleStep);

			// Calculate the offset for the vertex in the ring's plane
			FVector Offset = (FMath::Cos(Angle) * Right + FMath::Sin(Angle) * Forward) * Radius;

			// Compute the vertex position by adding the offset to the center
			FVector Vertex = Center + Offset;

			// Store the first vertex for potential use (e.g., for closed rings)
			if (i == 0) {
				FirstVertex = Vertex; // Store the first vertex
			}

			// Add the vertex to the ring's vertex array
			RingData.Vertices.Add(Vertex);
		}

		// calc the normal of the ring using right and forward vectors
		RingData.Normal = FVector::CrossProduct(Right, Forward).GetSafeNormal();

		// Store additional properties in the FRingData structure
		RingData.Center = Center;       // The center point of the ring
		RingData.Direction = Direction; // The direction vector of the ring
		RingData.Radius = Radius;       // The radius of the ring
		RingData.UpVector = UpVector;   // The up vector used for orientation
		RingData.NumVertices = RingData.Vertices.Num(); // Num verts for convenience
	}

	TArray<FVector> OrderRingVertices(const TArray<FVector>& InputVertices) {
		FVector Centroid = ComputeCentroid(InputVertices);

		// array to store vertices with angles
		TArray<TPair<float, FVector>> VerticesWithAngles;

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

	void FindIntersectionRingCardinalPoints(FIntersectionRingData& IntersectionRing, const FVector& StartCenter, const FVector& EndCenter) {
		int32 ArraySize = IntersectionRing.CombinedVertices.Num();

		if (ArraySize == 0) {
			UE_LOG(LogTemp, Error, TEXT("CombinedVertices array is empty!"));
			return;
		}

		// Find highest (north) and lowest (south) points
		int32 LowestIndex, HighestIndex;
		FVector LowestVertex, HighestVertex;
		DoobGeometryUtils::FindClosestVertex(StartCenter, IntersectionRing.CombinedVertices, LowestVertex, LowestIndex);
		DoobGeometryUtils::FindClosestVertex(EndCenter, IntersectionRing.CombinedVertices, HighestVertex, HighestIndex);

		// Calculate the ring's center
		FVector RingCenter = ComputeCentroid(IntersectionRing.CombinedVertices);

		// Define the ring's axis (from StartCenter to EndCenter)
		FVector RingAxis = (EndCenter - StartCenter).GetSafeNormal();

		// Define a perpendicular axis
		FVector UpVector = FVector::UpVector; // Can be replaced with another consistent vector
		FVector PerpendicularAxis = FVector::CrossProduct(UpVector, RingAxis).GetSafeNormal();

		// Ensure the PerpendicularAxis is valid
		if (PerpendicularAxis.IsZero()) {
			PerpendicularAxis = FVector::CrossProduct(FVector::RightVector, RingAxis).GetSafeNormal();
		}

		// Find left and right points based on the perpendicular axis
		float MaxProjection = -FLT_MAX;
		float MinProjection = FLT_MAX;
		FVector LeftVertex, RightVertex;
		int32 LeftIndex = -1, RightIndex = -1;

		for (int32 i = 0; i < ArraySize; i++) {
			float Projection = FVector::DotProduct(IntersectionRing.CombinedVertices[i] - RingCenter, PerpendicularAxis);

			if (Projection > MaxProjection) {
				MaxProjection = Projection;
				RightVertex = IntersectionRing.CombinedVertices[i];
				RightIndex = i;
			}

			if (Projection < MinProjection) {
				MinProjection = Projection;
				LeftVertex = IntersectionRing.CombinedVertices[i];
				LeftIndex = i;
			}
		}

		// Assign cardinal points
		IntersectionRing.CardinalVertices = { HighestVertex, RightVertex, LowestVertex, LeftVertex };
		IntersectionRing.CardinalIndices = { HighestIndex, RightIndex, LowestIndex, LeftIndex };
	}

	int32 FindRingIndexByCenter(const TArray<FRingData>& Rings, const FVector& Center, float Tolerance)
	{
		for (int32 Index = 0; Index < Rings.Num(); ++Index)
		{
			if (FVector::DistSquared(Rings[Index].Center, Center) <= FMath::Square(Tolerance))
			{
				return Index; // Return the index if the center matches within the tolerance
			}
		}
		return INDEX_NONE; // Return -1 if no match is found
	}

	void RemoveVerticesByInterpolatedDirections(
		FTwoTubeIntersectionData& TubeIntersectionData,
		FIntersectionSquareData& IntersectionSquare
	) {
		if (IntersectionSquare.Corners.Num() != 4) {
			UE_LOG(LogTemp, Warning, TEXT("Square must have exactly 4 vertices"));
			return;
		}

		// calc center of square
		CalculateSquareCenter(IntersectionSquare);

		// get perpendicular directions for top and bottom rings
		TArray<FVector> TopDirections, BottomDirections;
		TopDirections.Add((IntersectionSquare.Corners[3] - TubeIntersectionData.MTAboveIntersectionRing.Center).GetSafeNormal());
		TopDirections.Add((IntersectionSquare.Corners[2] - TubeIntersectionData.MTAboveIntersectionRing.Center).GetSafeNormal());

		BottomDirections.Add((IntersectionSquare.Corners[0] - TubeIntersectionData.MTBelowIntersectionRing.Center).GetSafeNormal());
		BottomDirections.Add((IntersectionSquare.Corners[1] - TubeIntersectionData.MTBelowIntersectionRing.Center).GetSafeNormal());

		FRingData TempTopRing = TubeIntersectionData.MTAboveIntersectionRing;
		TempTopRing.Vertices.Empty();

		// Flags to track corner insertion
		bool Corner2Added = false;
		bool Corner3Added = false;
		bool Corner23Added = false;

		int32 indexCount = 0;
		int32 indexCountFinal = 0;

		for (const FVector& Vertex : TubeIntersectionData.MTAboveIntersectionRing.Vertices) {
			FVector VertexDirection = (Vertex - TubeIntersectionData.MTAboveIntersectionRing.Center).GetSafeNormal();

			// check if the vertex direction is between the two interpolated directions
			FVector Normal = FVector::CrossProduct(TopDirections[0], TopDirections[1]).GetSafeNormal(); // Normal of the plane
			bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(TopDirections[0], VertexDirection)) >= 0 &&
				FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, TopDirections[1])) >= 0;

			if (!IsBetween) {
				if (!Corner3Added) {
					Corner3Added = true;
				}
				TempTopRing.Vertices.Add(Vertex);
				indexCount++;
			}
			else if (Corner3Added) {
				TempTopRing.Vertices.Add(IntersectionSquare.Corners[3]);
				Corner3Added = false;
				Corner2Added = true;
				indexCount++;
			}
			else if (Corner2Added) {
				TempTopRing.Vertices.Add(IntersectionSquare.Corners[2]);
				Corner2Added = false;
				Corner23Added = true;
				indexCountFinal = indexCount;
			}
		}

		TArray<FVector> ReorderedTopVertices;

		if (!Corner23Added) {
			FVector TempTestVertex = FindIntersectionOnRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TopDirections[1], TubeIntersectionData.MTAboveIntersectionRing.Center);

			for (int32 i = 0; i < TubeIntersectionData.MTAboveIntersectionRing.Vertices.Num(); ++i) {
				int32 TempNextVertexIndex = (i + 1) % TubeIntersectionData.MTAboveIntersectionRing.Vertices.Num();
				FVector TempCurrentVertex = TubeIntersectionData.MTAboveIntersectionRing.Vertices[i];
				FVector TempNextVertex = TubeIntersectionData.MTAboveIntersectionRing.Vertices[TempNextVertexIndex];

				if (IsPointOnLine(TempCurrentVertex, TempNextVertex, TempTestVertex)) {
					indexCountFinal = i;
					break;
				}
			}

			TArray<FVector> TempModifiedVertices = DoobContainerUtils::ReorderedArray(TempTopRing.Vertices, indexCountFinal + 1);
			FVector TempLastVertex = FindIntersectionOnRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TopDirections[0], TubeIntersectionData.MTAboveIntersectionRing.Center);

			ReorderedTopVertices.Add(TempTestVertex);
			ReorderedTopVertices.Append(TempModifiedVertices);
			ReorderedTopVertices.Add(TempLastVertex);
		}
		else {
			ReorderedTopVertices = DoobContainerUtils::ReorderedArray(TempTopRing.Vertices, indexCountFinal);
		}

		TempTopRing.Vertices = ReorderedTopVertices;

		FRingData TempBottomRing = TubeIntersectionData.MTBelowIntersectionRing;
		TempBottomRing.Vertices.Empty();

		indexCount = 0;
		indexCountFinal = 0;

		// Flags to track corner insertion
		bool Corner0Added = false;
		bool Corner1Added = false;
		bool Corner01Added = false;
		for (const FVector& Vertex : TubeIntersectionData.MTBelowIntersectionRing.Vertices) {
			FVector VertexDirection = (Vertex - TubeIntersectionData.MTBelowIntersectionRing.Center).GetSafeNormal();

			// check if the vertex direction is between the two interpolated directions
			FVector Normal = FVector::CrossProduct(BottomDirections[0], BottomDirections[1]).GetSafeNormal(); // Normal of the plane
			bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(BottomDirections[0], VertexDirection)) >= 0 &&
				FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, BottomDirections[1])) >= 0;

			if (!IsBetween) {
				if (!Corner0Added) {
					Corner0Added = true;
				}
				TempBottomRing.Vertices.Add(Vertex);
				indexCount++;
			}
			else if (Corner0Added) {
				TempBottomRing.Vertices.Add(IntersectionSquare.Corners[0]);
				Corner0Added = false;
				Corner1Added = true;
				indexCount++;
			}
			else if (Corner1Added) {
				TempBottomRing.Vertices.Add(IntersectionSquare.Corners[1]);
				Corner1Added = false;
				Corner01Added = true;
				indexCountFinal = indexCount;
			}
		}

		TArray<FVector> ReorderedBottomVertices;

		if (!Corner01Added) {
			FVector TempTestVertex = FindIntersectionOnRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, BottomDirections[1], TubeIntersectionData.MTBelowIntersectionRing.Center);

			for (int32 i = 0; i < TubeIntersectionData.MTBelowIntersectionRing.Vertices.Num(); ++i) {
				int32 TempNextVertexIndex = (i + 1) % TubeIntersectionData.MTBelowIntersectionRing.Vertices.Num();
				FVector TempCurrentVertex = TubeIntersectionData.MTBelowIntersectionRing.Vertices[i];
				FVector TempNextVertex = TubeIntersectionData.MTBelowIntersectionRing.Vertices[TempNextVertexIndex];

				if (IsPointOnLine(TempCurrentVertex, TempNextVertex, TempTestVertex)) {
					indexCountFinal = i;
					break;
				}
			}

			TArray<FVector> TempModifiedVertices = DoobContainerUtils::ReorderedArray(TempBottomRing.Vertices, indexCountFinal + 1);
			FVector TempLastVertex = FindIntersectionOnRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, BottomDirections[0], TubeIntersectionData.MTBelowIntersectionRing.Center);

			ReorderedBottomVertices.Add(TempTestVertex);
			ReorderedBottomVertices.Append(TempModifiedVertices);
			ReorderedBottomVertices.Add(TempLastVertex);
		}
		else {
			ReorderedBottomVertices = DoobContainerUtils::ReorderedArray(TempBottomRing.Vertices, indexCountFinal);
		}

		TempBottomRing.Vertices = ReorderedBottomVertices;

		TempBottomRing.bIsClosed = false;
		TempBottomRing.bIsComplete = false;
		TempTopRing.bIsClosed = false;
		TempTopRing.bIsComplete = false;

		TubeIntersectionData.MTAboveIntersectionRingPartial = TempTopRing;
		TubeIntersectionData.MTBelowIntersectionRingPartial = TempBottomRing;

		// iterate through rings between top and bottom
		FRingData MainTubeTopRing = TubeIntersectionData.MTAboveIntersectionRing;
		FRingData MainTubeBottomRing = TubeIntersectionData.MTBelowIntersectionRing;
		FVector TubeDirection = (MainTubeTopRing.Center - MainTubeBottomRing.Center).GetSafeNormal();

		float SegmentLength = FVector::Dist(MainTubeBottomRing.Center, MainTubeTopRing.Center);

		TArray<FRingData> TempTube;

		for (const FRingData& Ring : TubeIntersectionData.MainTube.Rings) {
			float RingProjection = FVector::DotProduct(Ring.Center - MainTubeBottomRing.Center, TubeDirection);
			float TopProjection = FVector::DotProduct(MainTubeTopRing.Center - MainTubeBottomRing.Center, TubeDirection);
			float BottomProjection = FVector::DotProduct(MainTubeBottomRing.Center - MainTubeTopRing.Center, TubeDirection);

			FVector PointToCurrentRing = Ring.Center - MainTubeBottomRing.Center;
			float PointLength = FVector::DotProduct(PointToCurrentRing, (MainTubeTopRing.Center - MainTubeBottomRing.Center).GetSafeNormal());

			bool RingInBounds = PointLength >= 0 && PointLength <= SegmentLength;

			if (!RingInBounds) {
				TempTube.Add(Ring);
				continue;
			}

			// blend factor for this ring
			float BlendFactor = (RingProjection - 0.0f) / (TopProjection - 0.0f);

			// interpolate directions
			FVector InterpolatedDirection1 = FMath::Lerp(BottomDirections[0], TopDirections[0], BlendFactor).GetSafeNormal();
			FVector InterpolatedDirection2 = FMath::Lerp(BottomDirections[1], TopDirections[1], BlendFactor).GetSafeNormal();

			FRingData ModifiedRing = Ring;
			ModifiedRing.Vertices.Empty();

			indexCount = 0;
			indexCountFinal = 0;

			bool LeftRingConnection = false;
			bool RightRingConnection = false;
			bool BothConnection = false;

			// check and remove vertices
			for (const FVector& Vertex : Ring.Vertices) {
				FVector VertexDirection = (Vertex - Ring.Center).GetSafeNormal();

				// check if the vertex direction is between the two interpolated directions
				FVector Normal = FVector::CrossProduct(InterpolatedDirection1, InterpolatedDirection2).GetSafeNormal(); // Normal of the plane
				bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(InterpolatedDirection1, VertexDirection)) >= 0 &&
					FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, InterpolatedDirection2)) >= 0;

				if (!IsBetween) {
					if (!LeftRingConnection) {
						LeftRingConnection = true;
					}
					ModifiedRing.Vertices.Add(Vertex);
					indexCount++;
				}
				else if (LeftRingConnection) {
					FVector LeftRingConnectionVertex = FindIntersectionOnRing(Ring.Vertices, InterpolatedDirection1, Ring.Center);
					ModifiedRing.Vertices.Add(LeftRingConnectionVertex);
					LeftRingConnection = false;
					RightRingConnection = true;
					indexCount++;
				}
				else if (RightRingConnection) {
					FVector RightRingConnectionVertex = FindIntersectionOnRing(Ring.Vertices, InterpolatedDirection2, Ring.Center);
					ModifiedRing.Vertices.Add(RightRingConnectionVertex);
					RightRingConnection = false;
					BothConnection = true;
					indexCountFinal = indexCount;
				}
			}

			TArray<FVector> ReorderedModifiedVertices;

			if (!BothConnection) {
				FVector TempTestVertex = FindIntersectionOnRing(Ring.Vertices, InterpolatedDirection2, Ring.Center);

				for (int32 i = 0; i < Ring.Vertices.Num(); ++i) {
					int32 TempNextVertexIndex = (i + 1) % Ring.Vertices.Num();
					FVector TempCurrentVertex = Ring.Vertices[i];
					FVector TempNextVertex = Ring.Vertices[TempNextVertexIndex];

					if (IsPointOnLine(TempCurrentVertex, TempNextVertex, TempTestVertex)) {
						indexCountFinal = i;
						break;
					}
				}

				TArray<FVector> TempModifiedVertices = DoobContainerUtils::ReorderedArray(ModifiedRing.Vertices, indexCountFinal + 1);
				FVector TempLastVertex = FindIntersectionOnRing(Ring.Vertices, InterpolatedDirection1, Ring.Center);

				ReorderedModifiedVertices.Add(TempTestVertex);
				ReorderedModifiedVertices.Append(TempModifiedVertices);
				ReorderedModifiedVertices.Add(TempLastVertex);
			}
			else {
				ReorderedModifiedVertices = DoobContainerUtils::ReorderedArray(ModifiedRing.Vertices, indexCountFinal);
			}

			ModifiedRing.Vertices = ReorderedModifiedVertices;
			ModifiedRing.bIsClosed = false;
			ModifiedRing.bIsComplete = false;

			TempTube.Add(ModifiedRing);
		}

		TubeIntersectionData.MainTube.Rings = TempTube;
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                                 3. Tube Operations                                                                   //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

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
			}
			else {
				segmentLength = t * Length;
				CurrentPosition = StartPosition + segmentLength * NormalizedDirection;
			}
		}

		OutTube.Direction = Direction;
		OutTube.UpVector = UpVector;
		OutTube.EndPosition = CurrentPosition;
		OutTube.Length = Length;
		OutTube.NumSegments = NumSegments;
		OutTube.NumSides = NumSides;
		OutTube.Profile = Profile;
		OutTube.StartPosition = StartPosition;
	}

	void RemoveInternalVertices(const FTubeData& TubeA, FTubeData& TubeB) {
		TArray<FRingData> TempTube;

		// Iterate through all rings in TubeB
		for (int32 TubeRingIndexB = 0; TubeRingIndexB < TubeB.Rings.Num(); ++TubeRingIndexB) {
			FRingData CurrentRingB = TubeB.Rings[TubeRingIndexB];
			FRingData TempRing;

			// Check each vertex in CurrentRingB against all frustums in TubeA
			for (int32 VertexIndex = 0; VertexIndex < CurrentRingB.Vertices.Num(); ++VertexIndex) {
				FVector CurrentVertex = CurrentRingB.Vertices[VertexIndex];
				bool bInsideAnyFrustum = false;

				// Test against all frustums defined by TubeA's consecutive rings
				for (int32 TubeRingIndexA = 0; TubeRingIndexA < TubeA.Rings.Num() - 1; ++TubeRingIndexA) {
					const FRingData& CurrentRingA = TubeA.Rings[TubeRingIndexA];
					const FRingData& NextRingA = TubeA.Rings[TubeRingIndexA + 1];

					if (IsPointInsideFrustum(CurrentRingA, NextRingA, CurrentVertex)) {
						bInsideAnyFrustum = true;
						break; // No need to check further if inside any frustum
					}
				}

				// Only keep vertices that are outside all frustums
				if (!bInsideAnyFrustum) {
					TempRing.Vertices.Add(CurrentVertex);
				}
			}

			// If TempRing has vertices, add it to the temporary tube
			if (TempRing.Vertices.Num() > 0) {
				TempRing.Center = CurrentRingB.Center;
				TempRing.Direction = CurrentRingB.Direction;
				TempRing.Normal = CurrentRingB.Normal;
				TempRing.Radius = CurrentRingB.Radius;
				TempRing.RingID = CurrentRingB.RingID;
				TempRing.UpVector = CurrentRingB.UpVector;
				TempRing.bIsClosed = (TempRing.Vertices.Num() == CurrentRingB.Vertices.Num());
				TempRing.bIsComplete = (TempRing.Vertices.Num() == CurrentRingB.Vertices.Num());

				TempTube.Add(TempRing);
			}
		}

		TubeB.Rings = TempTube;
	}

	void FindSegmentForPoint(
		const FVector& Point,
		const FTubeData& Tube,
		FRingData& StartRing,
		FRingData& EndRing
	) {
		// iterate through the rings to find the segments the point is between
		for (int32 i = 0; i < Tube.Rings.Num() - 1; ++i) {
			// get the current ring and the next ring
			FRingData CurrentRing = Tube.Rings[i];
			FRingData NextRing = Tube.Rings[i + 1];

			// calc the distance between the 2 rings
			float SegmentLength = FVector::Dist(CurrentRing.Center, NextRing.Center);

			// calc the length of the point relative to current ring - for interpolation
			FVector PointToCurrentRing = Point - CurrentRing.Center;
			float PointLength = FVector::DotProduct(PointToCurrentRing, (NextRing.Center - CurrentRing.Center).GetSafeNormal());

			// check if point is within segments
			if (PointLength >= 0 && PointLength <= SegmentLength) {
				StartRing = CurrentRing;
				EndRing = NextRing;
			}
		}
	}

	void RemoveVerticesInsideSquareByAngle(
		FTwoTubeIntersectionData& TubeIntersectionData,
		FIntersectionSquareData& IntersectionSquare
	) {
		// Ensure the square has 4 vertices
		if (IntersectionSquare.Corners.Num() != 4) {
			UE_LOG(LogTemp, Warning, TEXT("Square must have exactly 4 vertices"));
			return;
		}

		CalculateSquareCenter(IntersectionSquare);

		GetSquareAngles(IntersectionSquare);

		FRingData MainTubeTopRing = TubeIntersectionData.MTAboveIntersectionRing;
		FRingData MainTubeBottomRing = TubeIntersectionData.MTBelowIntersectionRing;

		TArray<FRingData> TempTube;

		// calc direction
		FVector TubeDirection = (MainTubeTopRing.Center - MainTubeBottomRing.Center).GetSafeNormal();

		// loop through main tube rings
		for (const FRingData& Ring : TubeIntersectionData.MainTube.Rings) {
			// project center of ring onto direction vector
			float RingProjection = FVector::DotProduct(Ring.Center - MainTubeBottomRing.Center, TubeDirection);
			float BottomProjection = FVector::DotProduct(MainTubeBottomRing.Center - MainTubeBottomRing.Center, TubeDirection);
			float TopProjection = FVector::DotProduct(MainTubeTopRing.Center - MainTubeBottomRing.Center, TubeDirection);

			// skip rings outside of intersection rings
			if (RingProjection < BottomProjection || RingProjection > TopProjection) {
				TempTube.Add(Ring);
				continue;
			}

			// create a modified version of the ring
			FRingData ModifiedRing = Ring;
			ModifiedRing.Vertices.Empty();

			// loop through the vertices in the current ring
			for (const FVector& Vertex : Ring.Vertices) {
				// skip vertices that lie inside the square's angular area
				if (!IsVertexInsideSquareAngle(Vertex, IntersectionSquare)) {
					// keep vertices outside the angular range
					ModifiedRing.Vertices.Add(Vertex);
				}
			}

			// add connection points to the square's sides
			const FVector& TopRightCorner = IntersectionSquare.Corners[2];
			const FVector& BottomRightCorner = IntersectionSquare.Corners[1];

			FVector RightSideDirection = (TopRightCorner - BottomRightCorner).GetSafeNormal();
			FVector RingToRightCorner = Ring.Center - BottomRightCorner;
			float RightT = FVector::DotProduct(RingToRightCorner, RightSideDirection);

			const FVector& TopLeftCorner = IntersectionSquare.Corners[3];
			const FVector& BottomLeftCorner = IntersectionSquare.Corners[0];

			FVector LeftSideDirection = (TopLeftCorner - BottomLeftCorner).GetSafeNormal();
			FVector RingToLeftCorner = Ring.Center - BottomLeftCorner;
			float LeftT = FVector::DotProduct(RingToLeftCorner, LeftSideDirection);

			// clamp t's to side lengths
			RightT = FMath::Clamp(RightT, 0.0f, (TopRightCorner - BottomRightCorner).Size());
			FVector RightIntersectionPoint = BottomRightCorner + RightSideDirection * RightT;

			LeftT = FMath::Clamp(LeftT, 0.0f, (TopLeftCorner - BottomLeftCorner).Size());
			FVector LeftIntersectionPoint = BottomLeftCorner + LeftSideDirection * LeftT;

			// add vertices to side of square
			IntersectionSquare.LeftSideRingConnections.Add(LeftIntersectionPoint);
			IntersectionSquare.RightSideRingConnections.Add(RightIntersectionPoint);

			// add the modified ring to the output
			TempTube.Add(Ring);
		}

		TubeIntersectionData.MainTube.Rings = TempTube;



		FRingData TempTopRing = TubeIntersectionData.MTAboveIntersectionRing;
		TempTopRing.Vertices.Empty();
		FRingData TempBottomRing = TubeIntersectionData.MTBelowIntersectionRing;
		TempBottomRing.Vertices.Empty();

		// loop through top ring
		for (const FVector& Vertex : TubeIntersectionData.MTAboveIntersectionRing.Vertices) {
			if (!IsVertexInsideSquareAngle(Vertex, IntersectionSquare)) {
				TempTopRing.Vertices.Add(Vertex);
			}
		}

		// loop through bottom ring
		for (const FVector& Vertex : TubeIntersectionData.MTBelowIntersectionRing.Vertices) {
			if (!IsVertexInsideSquareAngle(Vertex, IntersectionSquare)) {
				TempBottomRing.Vertices.Add(Vertex);
			}
		}

		TubeIntersectionData.MTAboveIntersectionRingPartial = TempTopRing;
		TubeIntersectionData.MTBelowIntersectionRingPartial = TempBottomRing;
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                            4. Intersection Calculations                                                              //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

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

		RemoveDuplicateVertices(CombinedVertices);

		OutRingData.CombinedVertices = OrderRingVertices(CombinedVertices);

		OutRingData.Centroid = ComputeCentroid(OutRingData.CombinedVertices);
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

			for (int32 VertexIndexA = 0; VertexIndexA < NumVerticesA; ++VertexIndexA) {
				FVector LineStartA = CurrentRingA.Vertices[VertexIndexA];
				FVector LineEndA = NextRingA.Vertices[VertexIndexA];

				FRingData StartRingA;
				FRingData EndRingA;
				FRingData StartRingB;
				FRingData EndRingB;

				FindSegmentForPoint(LineStartA, TubeB, StartRingA, EndRingA);
				FindSegmentForPoint(LineEndA, TubeB, StartRingB, EndRingB);

				TArray<FRingData> TempTubeB;
				TempTubeB.Add(StartRingA);
				TempTubeB.Add(EndRingA);
				TempTubeB.Add(StartRingB);
				TempTubeB.Add(EndRingB);

				// loop through TubeB rectangles
				for (int32 RingIndexB = 0; RingIndexB < TempTubeB.Num() - 1; ++RingIndexB) {
					const FRingData& CurrentRingB = TempTubeB[RingIndexB];
					const FRingData& NextRingB = TempTubeB[RingIndexB + 1];
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

	void GenerateSquareAroundIntersection(
		const FRingData& LowestRing,
		const FRingData& HighestRing,
		const FIntersectionRingData& IntersectionRing,
		int32 LeftIndex,
		int32 RightIndex,
		TArray<FVector>& OutSquareVertices
	) {
		// Validate input
		if (IntersectionRing.CombinedVertices.Num() == 0 || LowestRing.Vertices.Num() == 0 || HighestRing.Vertices.Num() == 0) {
			UE_LOG(LogTemp, Warning, TEXT("Invalid input data for square generation"));
			return;
		}

		// Get leftmost and rightmost vertices
		FVector LeftVertex = IntersectionRing.CombinedVertices[LeftIndex];
		FVector RightVertex = IntersectionRing.CombinedVertices[RightIndex];

		// calc center
		FVector PointCenter = CalculateCenterLinePoint(LeftVertex, LowestRing, HighestRing);

		// calc directions for intersections
		FVector LeftDirection = (LeftVertex - PointCenter).GetSafeNormal();
		FVector RightDirection = (RightVertex - PointCenter).GetSafeNormal();

		// find intersection on the lowest and highest rings
		FVector BottomLeft = FindIntersectionOnRing(LowestRing.Vertices, LeftDirection, LowestRing.Center);
		FVector BottomRight = FindIntersectionOnRing(LowestRing.Vertices, RightDirection, LowestRing.Center);
		FVector TopLeft = FindIntersectionOnRing(HighestRing.Vertices, LeftDirection, HighestRing.Center);
		FVector TopRight = FindIntersectionOnRing(HighestRing.Vertices, RightDirection, HighestRing.Center);

		// store vertices in output array
		OutSquareVertices = { BottomLeft, BottomRight, TopRight, TopLeft };
	}

	void GenerateAboveBelowIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData) {
		FRingData LowVertStartRing;
		FRingData LowVertEndRing;
		FRingData HighVertStartRing;
		FRingData HighVertEndRing;

		// Find the segments in the main tube for the lowest and highest intersection points
		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.CardinalVertices[2], TubeIntersectionData.MainTube, LowVertStartRing, LowVertEndRing);
		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.CardinalVertices[0], TubeIntersectionData.MainTube, HighVertStartRing, HighVertEndRing);

		// Calculate the center points for the lowest and highest intersection vertices
		FVector LowestVertexCenter = CalculateCenterLinePoint(TubeIntersectionData.IntersectionRing.CardinalVertices[2], LowVertStartRing, LowVertEndRing);
		FVector HighestVertexCenter = CalculateCenterLinePoint(TubeIntersectionData.IntersectionRing.CardinalVertices[0], HighVertStartRing, HighVertEndRing);

		// Interpolate the ring radius for both the lowest and highest vertex centers
		float LowRingRadius = InterpolatedRingRadius(LowestVertexCenter, LowVertStartRing, LowVertEndRing);
		float HighRingRadius = InterpolatedRingRadius(HighestVertexCenter, HighVertStartRing, HighVertEndRing);

		GenerateRing(
			LowestVertexCenter,
			TubeIntersectionData.MainTube.Direction,
			TubeIntersectionData.MainTube.UpVector,
			LowRingRadius,
			TubeIntersectionData.MainTube.NumSides,
			TubeIntersectionData.MTBelowIntersectionRing
		);

		GenerateRing(
			HighestVertexCenter,
			TubeIntersectionData.MainTube.Direction,
			TubeIntersectionData.MainTube.UpVector,
			HighRingRadius,
			TubeIntersectionData.MainTube.NumSides,
			TubeIntersectionData.MTAboveIntersectionRing
		);
	}

	void GenerateLateralTubeIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData) {
		TubeIntersectionData.LateralTubeIntersectionRings = TubeIntersectionData.LateralTube;
		TubeIntersectionData.LateralTubeIntersectionRings.Rings.Empty();

		// create temp ring to reorder ring to north point
		TArray<FVector> ReorderedIntersection = DoobContainerUtils::ReorderedArray(
			TubeIntersectionData.IntersectionRing.CombinedVertices,
			TubeIntersectionData.IntersectionRing.CardinalIndices[0]
		);

		int32 StartIndex = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[1]) - 1;
		int32 EndIndexAlt = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[3]) + 1;

		TArray<FVector> NorthSouthVerts = {
			TubeIntersectionData.IntersectionRing.CardinalVertices[0],
			TubeIntersectionData.IntersectionRing.CardinalVertices[2]
		};

		FVector ClosestVertToEnd;
		int32 ClosestVertToEndIndex;

		FindClosestVertex(TubeIntersectionData.LateralTube.EndPosition, NorthSouthVerts, ClosestVertToEnd, ClosestVertToEndIndex);

		if (ClosestVertToEndIndex > 0) {
			ReorderedIntersection.Empty();
			ReorderedIntersection = DoobContainerUtils::ReorderedArray(
				TubeIntersectionData.IntersectionRing.CombinedVertices,
				TubeIntersectionData.IntersectionRing.CardinalIndices[2]
			);
		}

		FRingData IntersectionRingStartRing;
		FRingData IntersectionRingEndRing;

		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.Centroid, TubeIntersectionData.LateralTube, IntersectionRingStartRing, IntersectionRingEndRing);

		float IntersectionRingRadius = InterpolatedRingRadius(TubeIntersectionData.IntersectionRing.Centroid, IntersectionRingStartRing, IntersectionRingEndRing);

		FRingData FirstLateralIntersectionRing;

		GenerateRing(
			TubeIntersectionData.IntersectionRing.Centroid,
			TubeIntersectionData.LateralTube.Direction,
			TubeIntersectionData.LateralTube.UpVector,
			IntersectionRingRadius,
			TubeIntersectionData.LateralTube.NumSides,
			FirstLateralIntersectionRing
		);

		FirstLateralIntersectionRing.Vertices.Empty();

		for (int32 i = StartIndex; i <= EndIndexAlt; ++i) {
			FVector VertexToAdd = ReorderedIntersection[i];

			FirstLateralIntersectionRing.Vertices.Add(VertexToAdd);

			FVector SuperTempCenter = CalculateCenterLinePoint(VertexToAdd, IntersectionRingStartRing, IntersectionRingEndRing);
			FVector SuperTempDirection = (VertexToAdd - SuperTempCenter).GetSafeNormal();
		}

		int32 EndIndex = 0;

		FRingData CurrentLateralIntersectionRing;
		FRingData PreviousLateralIntersectionRing;

		for (int32 i = StartIndex; i >= EndIndex; --i) {
			FRingData CurrentStartRing;
			FRingData CurrentEndRing;

			TArray<FVector> TempCurrentVertices;
			TArray<FVector> TempPreviousVertices;

			FindSegmentForPoint(ReorderedIntersection[i], TubeIntersectionData.LateralTube, CurrentStartRing, CurrentEndRing);

			FVector CurrentVertexCenterPoint = CalculateCenterLinePoint(ReorderedIntersection[i], CurrentStartRing, CurrentEndRing);

			float CurrentVertexRingRadius = InterpolatedRingRadius(CurrentVertexCenterPoint, CurrentStartRing, CurrentEndRing);

			if (i == StartIndex) {
				PreviousLateralIntersectionRing = FirstLateralIntersectionRing;
			}
			else {
				PreviousLateralIntersectionRing = CurrentLateralIntersectionRing;
				int32 ReorderedIntersectionNum = ReorderedIntersection.Num();
				TempPreviousVertices.Add(ReorderedIntersection[i]);
				TempPreviousVertices.Append(CurrentLateralIntersectionRing.Vertices);

				if (EndIndexAlt >= ReorderedIntersectionNum) {
					EndIndexAlt = 0;
				}

				TempPreviousVertices.Add(ReorderedIntersection[EndIndexAlt]);
				PreviousLateralIntersectionRing.Vertices = TempPreviousVertices;
			}

			EndIndexAlt++;

			CurrentLateralIntersectionRing.Vertices.Empty();

			GenerateRing(
				CurrentVertexCenterPoint,
				TubeIntersectionData.LateralTube.Direction,
				TubeIntersectionData.LateralTube.UpVector,
				CurrentVertexRingRadius,
				TubeIntersectionData.LateralTube.NumSides,
				CurrentLateralIntersectionRing
			);

			TubeIntersectionData.LateralTubeFirstFullRing = CurrentLateralIntersectionRing;

			FRingData TempNextStartRing;
			FRingData TempNextEndRing;
			FVector TempNextVertex;

			for (int32 j = 0; j < PreviousLateralIntersectionRing.Vertices.Num(); ++j) {
				FRingData TempCurrentStartRing;
				FRingData TempCurrentEndRing;
				FindSegmentForPoint(PreviousLateralIntersectionRing.Vertices[j], TubeIntersectionData.LateralTube, TempCurrentStartRing, TempCurrentEndRing);
				FVector TempCurrentCenterlinePoint = CalculateCenterLinePoint(PreviousLateralIntersectionRing.Vertices[j], TempCurrentStartRing, TempCurrentEndRing);
				FVector TempCurrentDirection = (PreviousLateralIntersectionRing.Vertices[j] - TempCurrentCenterlinePoint).GetSafeNormal();

				FVector TempCurrentVertex = FindIntersectionOnNewRing(CurrentLateralIntersectionRing.Vertices, TempCurrentDirection, CurrentLateralIntersectionRing.Center, PreviousLateralIntersectionRing.Vertices[j]);

				FRingData CurrentMainTubeStartRing;
				FRingData CurrentMainTubeEndRing;

				FindSegmentForPoint(TempCurrentVertex, TubeIntersectionData.MainTube, CurrentMainTubeStartRing, CurrentMainTubeEndRing);

				int32 TempNextIndex = (j + 1) % PreviousLateralIntersectionRing.Vertices.Num();
				int32 TempPrevIndex = (j - 1 + PreviousLateralIntersectionRing.Vertices.Num()) % PreviousLateralIntersectionRing.Vertices.Num();

				bool IsInsideFrustum = IsPointInsideFrustum(CurrentMainTubeStartRing, CurrentMainTubeEndRing, TempCurrentVertex);
				bool IsAlongLine = IsPointOnLine(PreviousLateralIntersectionRing.Vertices[j], PreviousLateralIntersectionRing.Vertices[TempNextIndex], TempCurrentVertex) ||
					IsPointOnLine(PreviousLateralIntersectionRing.Vertices[TempPrevIndex], PreviousLateralIntersectionRing.Vertices[j], TempCurrentVertex);

				if (IsInsideFrustum && !IsAlongLine) {
					TempCurrentVertices.Add(PreviousLateralIntersectionRing.Vertices[j]);
				}
				else {
					TempCurrentVertices.Add(TempCurrentVertex);
				}

				FindSegmentForPoint(PreviousLateralIntersectionRing.Vertices[TempNextIndex], TubeIntersectionData.LateralTube, TempNextStartRing, TempNextEndRing);
				FVector TempNextCenterlinePoint = CalculateCenterLinePoint(PreviousLateralIntersectionRing.Vertices[TempNextIndex], TempNextStartRing, TempNextEndRing);
				FVector TempNextDirection = (PreviousLateralIntersectionRing.Vertices[TempNextIndex] - TempNextCenterlinePoint).GetSafeNormal();

				TempNextVertex = FindIntersectionOnNewRing(CurrentLateralIntersectionRing.Vertices, TempNextDirection, CurrentLateralIntersectionRing.Center, PreviousLateralIntersectionRing.Vertices[TempNextIndex]);
			}

			CurrentLateralIntersectionRing.Vertices = TempCurrentVertices;
			CurrentLateralIntersectionRing.bIsClosed = false;
			CurrentLateralIntersectionRing.bIsComplete = false;

			PreviousLateralIntersectionRing.bIsClosed = false;
			PreviousLateralIntersectionRing.bIsComplete = false;

			TubeIntersectionData.LateralTubeIntersectionRings.Rings.Add(PreviousLateralIntersectionRing);
			TubeIntersectionData.LateralTubeIntersectionRings.Rings.Add(CurrentLateralIntersectionRing);
		}
	}

	FVector FindIntersectionOnRing(const TArray<FVector>& RingVertices, const FVector& Direction, const FVector& Center) {


		FVector NormalizedDirection = Direction.GetSafeNormal();
		FVector IntersectionPoint = FVector::ZeroVector;

		int32 VertexCount = RingVertices.Num();
		for (int32 i = 0; i < VertexCount; ++i) {
			// get the current line segment of the ring
			FVector Start = RingVertices[i];
			FVector End = RingVertices[(i + 1) % VertexCount];  // Wrap around to the start for the last segment

			// calc the line's vector
			FVector LineVector = End - Start;
			FVector LineToCenter = Center - Start;

			// calc the determinant for intersection test
			float Determinant = FVector::CrossProduct(LineVector, NormalizedDirection).Size();

			// if determinant is zero, the direction is parallel to the line
			if (FMath::IsNearlyZero(Determinant)) {
				continue;
			}

			// compute the intersection point using parametric line equation
			float T = FVector::CrossProduct(LineToCenter, NormalizedDirection).Size() / Determinant;
			float U = FVector::CrossProduct(LineToCenter, LineVector).Size() / Determinant;

			// check if the intersection point lies on the segment
			if (T >= 0.0f && T <= 1.0f && U >= 0.0f) {
				IntersectionPoint = Start + T * LineVector;
				break;
			}
		}

		return IntersectionPoint;
	}

	FVector FindIntersectionOnNewRing(const TArray<FVector>& RingVertices, const FVector& Direction, const FVector& Center, const FVector& TargetVertex) {
		FVector NormalizedDirection = Direction.GetSafeNormal();
		TArray<FVector> IntersectionPoints;

		int32 VertexCount = RingVertices.Num();

		for (int32 i = 0; i < VertexCount; ++i) {
			// Get the current line segment of the ring
			FVector Start = RingVertices[i];
			FVector End = RingVertices[(i + 1) % VertexCount]; // Wrap-around to start for the last segment

			// Calculate the line vector and vector from Start to Center
			FVector LineVector = End - Start;
			FVector LineToCenter = Center - Start;

			// Calculate the determinant to check parallelism
			float Determinant = FVector::CrossProduct(LineVector, NormalizedDirection).Size();
			if (FMath::IsNearlyZero(Determinant)) {
				continue; // Skip parallel lines
			}

			// Compute parametric values T and U
			float T = FVector::CrossProduct(LineToCenter, NormalizedDirection).Size() / Determinant;
			float U = FVector::CrossProduct(LineToCenter, LineVector).Size() / Determinant;

			// Check if the intersection lies within the segment bounds
			if (T >= 0.0f && T <= 1.0f && U >= 0.0f) {
				FVector IntersectionPoint = Start + T * LineVector;
				IntersectionPoints.Add(IntersectionPoint);
			}
		}

		// If no intersections were found, return zero vector
		if (IntersectionPoints.Num() == 0) {
			UE_LOG(LogTemp, Log, TEXT("No intersection found on the ring."));
			return FVector::ZeroVector;
		}

		// Find the intersection point closest to the target vertex
		FVector ClosestPoint = IntersectionPoints[0];
		float MinDistance = FVector::DistSquared(ClosestPoint, TargetVertex);

		for (const FVector& Point : IntersectionPoints) {
			float Distance = FVector::DistSquared(Point, TargetVertex);
			if (Distance < MinDistance) {
				ClosestPoint = Point;
				MinDistance = Distance;
			}
		}

		return ClosestPoint;
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

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                           5. Polygon and Shape Utilities                                                             //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

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

	void CalculateSquareCenter(FIntersectionSquareData& IntersectionSquare) {
		FVector Center(0.0f, 0.0f, 0.0f);
		for (const FVector& Vertex : IntersectionSquare.Corners) {
			Center += Vertex;
		}
		Center /= IntersectionSquare.Corners.Num();
		IntersectionSquare.Center = Center;
	}

	void GetSquareAngles(FIntersectionSquareData& IntersectionSquare) {
		TArray<float> Angles;

		// calc the square's normal vector
		FVector SquareNormal = FVector::CrossProduct(
			IntersectionSquare.Corners[1] - IntersectionSquare.Corners[0],
			IntersectionSquare.Corners[2] - IntersectionSquare.Corners[0]
		).GetSafeNormal();

		// create transformation matrix to local space
		FMatrix ToLocalSpace = FRotationMatrix::MakeFromZ(SquareNormal);

		// transform the square center to local space
		FVector LocalCenter = ToLocalSpace.TransformPosition(IntersectionSquare.Center);

		// transform corners to local space and calc angle
		for (const FVector& Corner : IntersectionSquare.Corners) {
			FVector LocalCorner = ToLocalSpace.TransformPosition(Corner);

			// calc the vector relative to the center in local space (ignoring Z-Axis)
			FVector2D Vector2D(LocalCorner.X - LocalCenter.X, LocalCorner.Y - LocalCenter.Y);

			// calc the angle
			float Angle = FMath::Atan2(Vector2D.Y, Vector2D.X);
			Angles.Add(Angle);
		}
		Angles.Sort();

		IntersectionSquare.Angles = Angles;
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                               6. Connection Utilities                                                                //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

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

	void ConnectPartialRings(const TArray<FVector>& RingA, const TArray<FVector>& RingB, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		for (int32 i = 0; i < RingA.Num() - 1; ++i) {
			int32 CurrentA = BaseIndex + i;
			int32 NextA = BaseIndex + (i + 1);

			int32 CurrentB = BaseIndex + RingA.Num() + i;
			int32 NextB = BaseIndex + RingA.Num() + (i + 1);

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

	void ConnectPartialRingArray(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		if (Rings.Num() < 2) return;

		// iterate through the array of rings and connect consecutive rings
		for (int32 i = 0; i < Rings.Num() - 1; ++i) {
			const TArray<FVector>& CurrentRing = Rings[i].Vertices;
			const TArray<FVector>& NextRing = Rings[i + 1].Vertices;

			ConnectPartialRings(CurrentRing, NextRing, Vertices, Triangles, BaseIndex);
		}

	}

	void ConnectPartialRingArrayPaired(const TArray<FRingData>& Rings, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		if (Rings.Num() < 2) return;

		// iterate through the array of rings and connect consecutive rings
		for (int32 i = 0; i < Rings.Num() - 1; i += 2) {
			const TArray<FVector>& CurrentRing = Rings[i].Vertices;
			const TArray<FVector>& NextRing = Rings[i + 1].Vertices;

			ConnectPartialRings(CurrentRing, NextRing, Vertices, Triangles, BaseIndex);
		}
	}

	void ConnectIntersectionCornerArrays(const TArray<TArray<FVector>>& CornerVertexArrays, TArray<FVector>& Vertices, TArray<int32>& Triangles, int32& BaseIndex) {
		if (CornerVertexArrays.Num() < 2) return;

		for (int32 i = 0; i < CornerVertexArrays.Num() - 1; ++i) {
			TArray<FVector> CurrentVertexArray = CornerVertexArrays[i];
			TArray<FVector> NextVertexArray = CornerVertexArrays[i + 1];

			for (int32 j = 0; j < CurrentVertexArray.Num() - 1; ++j) {
				int32 CurrentA = BaseIndex + j;
				int32 NextA = BaseIndex + (j + 1);

				int32 CurrentB = BaseIndex + CurrentVertexArray.Num() + j;
				int32 NextB = BaseIndex + CurrentVertexArray.Num() + (j + 1);

				// Triangle 1
				Triangles.Add(CurrentA);
				Triangles.Add(NextA);
				Triangles.Add(CurrentB);

				// Triangle 2
				Triangles.Add(CurrentB);
				Triangles.Add(NextA);
				Triangles.Add(NextB);
			}

			Vertices.Append(CurrentVertexArray);
			Vertices.Append(NextVertexArray);
			BaseIndex += CurrentVertexArray.Num() + NextVertexArray.Num();
		}
	}

	void ConnectIntersectionRingToSquare(FTwoTubeIntersectionData& TubeIntersectionData, int32& BaseIndex) {
		ConnectIntersectionCornerArrays(TubeIntersectionData.IntersectionSquare.TopRightPartialRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
		ConnectIntersectionCornerArrays(TubeIntersectionData.IntersectionSquare.BottomRightPartialRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
		ConnectIntersectionCornerArrays(TubeIntersectionData.IntersectionSquare.BottomLeftPartialRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
		ConnectIntersectionCornerArrays(TubeIntersectionData.IntersectionSquare.TopLeftPartialRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
	}

	void ConnectTwoTubeIntersection(FTwoTubeIntersectionData& TubeIntersectionData) {
		TArray<FRingData> BelowIntersectionRings;
		TArray<FRingData> PartialIntersectionRings;
		TArray<FRingData> AboveIntersectionRings;

		TArray<FRingData> FullLateralTubeRings;

		PartialIntersectionRings.Add(TubeIntersectionData.MTBelowIntersectionRingPartial);

		int32 MainTubeRingIndex = 0;
		int32 BaseIndex = 0;

		bool AboveBelowAdded = false;

		for (const FRingData Ring : TubeIntersectionData.MainTube.Rings) {

			if (Ring.bIsComplete && MainTubeRingIndex <= (TubeIntersectionData.MainTube.Rings.Num() / 2) - 1) {
				BelowIntersectionRings.Add(Ring);
			}
			else if (!Ring.bIsComplete) {
				PartialIntersectionRings.Add(Ring);
			}
			else if (Ring.bIsComplete && MainTubeRingIndex > (TubeIntersectionData.MainTube.Rings.Num() / 2) - 1) {
				if (!AboveBelowAdded) {
					BelowIntersectionRings.Add(TubeIntersectionData.MTBelowIntersectionRing);
					AboveIntersectionRings.Add(TubeIntersectionData.MTAboveIntersectionRing);
					AboveBelowAdded = true;
				}
				AboveIntersectionRings.Add(Ring);
			}

			MainTubeRingIndex++;
		}

		PartialIntersectionRings.Add(TubeIntersectionData.MTAboveIntersectionRingPartial);

		ConnectRingArray(BelowIntersectionRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);

		ConnectPartialRingArray(PartialIntersectionRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);

		ConnectIntersectionRingToSquare(TubeIntersectionData, BaseIndex);

		ConnectRingArray(AboveIntersectionRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);

		TArray<FRingData> ReversedLateralTubeIntersectionRings;

		for (int32 i = 0; i < TubeIntersectionData.LateralTubeIntersectionRings.Rings.Num(); ++i) {
			FRingData TempReverseRing = TubeIntersectionData.LateralTubeIntersectionRings.Rings[i];
			Algo::Reverse(TempReverseRing.Vertices);
			ReversedLateralTubeIntersectionRings.Add(TempReverseRing);
		}

		ConnectPartialRingArrayPaired(ReversedLateralTubeIntersectionRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);

		FullLateralTubeRings.Add(TubeIntersectionData.LateralTubeFirstFullRing);

		for (const FRingData Ring : TubeIntersectionData.LateralTubeRemovedVertices.Rings) {
			if (Ring.bIsComplete) {
				FullLateralTubeRings.Add(Ring);
			}
		}

		// this including the part in the struct was for testing only, can remove all of this
		TubeIntersectionData.MainTubePartialRings = FullLateralTubeRings;

		ConnectRingArray(FullLateralTubeRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
	}

	void OrderSquareIntersectionConnectionsOneCorner(
		const FTubeData& MainTubeIntersectionRings,
		const FTubeData& LateralTube,
		const TArray<FVector> RingVertices,
		const TArray<FVector> SquareVertices,
		const int32 StartIndex,
		TArray<TArray<FVector>>& OutVertexArrays,
		const bool Reversed
	) {
		if (MainTubeIntersectionRings.Rings.Num() == 0)
		{
			UE_LOG(LogTemp, Warning, TEXT("Array is empty."));
			return;
		}

		int32 CurrentIndex = StartIndex;

		TArray<TArray<FVector>> TempVertexArrays;
		TArray<FVector> PreviousVertexArray = RingVertices;
		TempVertexArrays.Add(PreviousVertexArray);

		while (true) {

			FRingData TempNextStartRing;
			FRingData TempNextEndRing;
			FVector TempNextVertex;

			TArray<FVector> CurrentVertexArray;

			for (int32 i = 0; i < PreviousVertexArray.Num(); ++i) {
				FRingData TempCurrentStartRing;
				FRingData TempCurrentEndRing;
				FindSegmentForPoint(PreviousVertexArray[i], MainTubeIntersectionRings, TempCurrentStartRing, TempCurrentEndRing);
				FVector TempCurrentCenterlinePoint = CalculateCenterLinePoint(PreviousVertexArray[i], TempCurrentStartRing, TempCurrentEndRing);
				FVector TempCurrentDirection = (PreviousVertexArray[i] - TempCurrentCenterlinePoint).GetSafeNormal();

				FVector TempCurrentVertex = FindIntersectionOnNewRing(MainTubeIntersectionRings.Rings[CurrentIndex].Vertices, TempCurrentDirection, MainTubeIntersectionRings.Rings[CurrentIndex].Center, PreviousVertexArray[i]);

				FRingData CurrentLateralTubeStartRing;
				FRingData CurrentLateralTubeEndRing;

				FindSegmentForPoint(TempCurrentVertex, LateralTube, CurrentLateralTubeStartRing, CurrentLateralTubeEndRing);

				int32 TempNextIndex = (i + 1) % PreviousVertexArray.Num();
				int32 TempPrevIndex = (i - 1 + PreviousVertexArray.Num()) % PreviousVertexArray.Num();

				bool IsInsideFrustum = IsPointInsideFrustum(CurrentLateralTubeStartRing, CurrentLateralTubeEndRing, TempCurrentVertex);
				bool IsAlongLine = IsPointOnLine(PreviousVertexArray[i], PreviousVertexArray[TempNextIndex], TempCurrentVertex) ||
					IsPointOnLine(PreviousVertexArray[TempPrevIndex], PreviousVertexArray[i], TempCurrentVertex);

				if (IsInsideFrustum && !IsAlongLine) {
					CurrentVertexArray.Add(PreviousVertexArray[i]);
				}
				else {
					CurrentVertexArray.Add(TempCurrentVertex);
				}
			}

			TempVertexArrays.Add(CurrentVertexArray);
			PreviousVertexArray = CurrentVertexArray;

			// Check the stopping condition
			if (Reversed)
			{
				if (CurrentIndex == MainTubeIntersectionRings.Rings.Num() - 1) {
					break;
				}
				++CurrentIndex; // Move forward
			}
			else
			{
				if (CurrentIndex == 0) {
					break;
				}
				--CurrentIndex; // Move backward
			}
		}

		TempVertexArrays.Add(SquareVertices);

		OutVertexArrays = TempVertexArrays;
	}

	void OrderSquareIntersectionConnections(FTwoTubeIntersectionData& TubeIntersectionData) {
		TArray<FVector> TempBottomLeftRingVertices;
		TArray<FVector> TempBottomRightRingVertices;
		TArray<FVector> TempTopLeftRingVertices;
		TArray<FVector> TempTopRightRingVertices;

		TArray<FVector> TempBottomLeftSquareVertices;
		TArray<FVector> TempBottomRightSquareVertices;
		TArray<FVector> TempTopLeftSquareVertices;
		TArray<FVector> TempTopRightSquareVertices;

		TArray<FVector> ReorderedIntersection = DoobContainerUtils::ReorderedArray(
			TubeIntersectionData.IntersectionRing.CombinedVertices,
			TubeIntersectionData.IntersectionRing.CardinalIndices[0]
		);

		int32 index = 0;
		int32 TotalIndices = ReorderedIntersection.Num();
		int32 HalfIndices = ReorderedIntersection.Num() / 2;
		int32 QuartIndices = ReorderedIntersection.Num() / 4;

		int32 RightIndex = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[1]);
		int32 BottomIndex = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[2]);
		int32 LeftIndex = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[3]);

		for (const FVector& Vertex : ReorderedIntersection) {
			FVector CenterlinePoint = CalculateCenterLinePoint(Vertex, TubeIntersectionData.MTBelowIntersectionRing, TubeIntersectionData.MTAboveIntersectionRing);
			FVector TempDirection = (Vertex - CenterlinePoint).GetSafeNormal();
			FVector TempSquarePoint;

			if (index < RightIndex) {
				TempTopRightRingVertices.Add(Vertex);
				TempSquarePoint = FindIntersectionOnNewRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TempDirection, TubeIntersectionData.MTAboveIntersectionRing.Center, Vertex);
				TempTopRightSquareVertices.Add(TempSquarePoint);
			}
			else if (index >= RightIndex && index < BottomIndex) {
				TempBottomRightRingVertices.Add(Vertex);
				TempSquarePoint = FindIntersectionOnNewRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, TempDirection, TubeIntersectionData.MTBelowIntersectionRing.Center, Vertex);
				TempBottomRightSquareVertices.Add(TempSquarePoint);
			}
			else if (index >= BottomIndex && index < LeftIndex) {
				TempBottomLeftRingVertices.Add(Vertex);
				TempSquarePoint = FindIntersectionOnNewRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, TempDirection, TubeIntersectionData.MTBelowIntersectionRing.Center, Vertex);
				TempBottomLeftSquareVertices.Add(TempSquarePoint);
			}
			else if (index >= LeftIndex) {
				TempTopLeftRingVertices.Add(Vertex);
				TempSquarePoint = FindIntersectionOnNewRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TempDirection, TubeIntersectionData.MTAboveIntersectionRing.Center, Vertex);
				TempTopLeftSquareVertices.Add(TempSquarePoint);
			}

			index++;
		}

		TempTopRightRingVertices.Add(ReorderedIntersection[RightIndex]);
		TempBottomRightRingVertices.Add(ReorderedIntersection[BottomIndex]);
		TempBottomLeftRingVertices.Add(ReorderedIntersection[LeftIndex]);
		TempTopLeftRingVertices.Add(ReorderedIntersection[0]);

		FVector TopRightCenterline = CalculateCenterLinePoint(
			ReorderedIntersection[RightIndex],
			TubeIntersectionData.MTBelowIntersectionRing,
			TubeIntersectionData.MTAboveIntersectionRing
		);

		FVector BottomRightCenterline = CalculateCenterLinePoint(
			ReorderedIntersection[BottomIndex],
			TubeIntersectionData.MTBelowIntersectionRing,
			TubeIntersectionData.MTAboveIntersectionRing
		);

		FVector BottomLeftCenterline = CalculateCenterLinePoint(
			ReorderedIntersection[LeftIndex],
			TubeIntersectionData.MTBelowIntersectionRing,
			TubeIntersectionData.MTAboveIntersectionRing
		);

		FVector TopLeftCenterline = CalculateCenterLinePoint(
			ReorderedIntersection[0],
			TubeIntersectionData.MTBelowIntersectionRing,
			TubeIntersectionData.MTAboveIntersectionRing
		);

		FVector TopRightDirection = (ReorderedIntersection[RightIndex] - TopRightCenterline).GetSafeNormal();
		FVector BottomRightDirection = (ReorderedIntersection[BottomIndex] - BottomRightCenterline).GetSafeNormal();
		FVector BottomLeftDirection = (ReorderedIntersection[LeftIndex] - BottomLeftCenterline).GetSafeNormal();
		FVector TopLeftDirection = (ReorderedIntersection[0] - TopLeftCenterline).GetSafeNormal();

		FVector TempTopRightPoint = FindIntersectionOnRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TopRightDirection, TubeIntersectionData.MTAboveIntersectionRing.Center);
		FVector TempBottomRightPoint = FindIntersectionOnRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, BottomRightDirection, TubeIntersectionData.MTBelowIntersectionRing.Center);
		FVector TempBottomLeftPoint = FindIntersectionOnRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, BottomLeftDirection, TubeIntersectionData.MTBelowIntersectionRing.Center);
		FVector TempTopLeftPoint = FindIntersectionOnRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TopLeftDirection, TubeIntersectionData.MTAboveIntersectionRing.Center);

		TempTopRightSquareVertices.Add(TempTopRightPoint);
		TempBottomRightSquareVertices.Add(TempBottomRightPoint);
		TempBottomLeftSquareVertices.Add(TempBottomLeftPoint);
		TempTopLeftSquareVertices.Add(TempTopLeftPoint);

		FVector AboveIntersectionCenter = TubeIntersectionData.MTAboveIntersectionRing.Center;
		FVector BelowIntersectionCenter = TubeIntersectionData.MTBelowIntersectionRing.Center;
		int32 MainTubeRingsNum = TubeIntersectionData.MainTube.Rings.Num();

		FTubeData TempMainTubeIntersectionRings = TubeIntersectionData.MainTube;
		TempMainTubeIntersectionRings.Rings.Empty();
		TempMainTubeIntersectionRings.Rings.Add(TubeIntersectionData.MTBelowIntersectionRing);

		for (int32 i = 0; i < MainTubeRingsNum; ++i) {
			FVector CurrentRingCenter = TubeIntersectionData.MainTube.Rings[i].Center;
			if (IsPointOnLine(BelowIntersectionCenter, AboveIntersectionCenter, CurrentRingCenter)) {
				TempMainTubeIntersectionRings.Rings.Add(TubeIntersectionData.MainTube.Rings[i]);
			}
		}

		TempMainTubeIntersectionRings.Rings.Add(TubeIntersectionData.MTAboveIntersectionRing);

		FRingData StartCenterRing;
		FRingData EndCenterRing;

		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.CardinalVertices[1], TempMainTubeIntersectionRings, StartCenterRing, EndCenterRing);

		int32 StartIntersectionRingsStartIndex = FindRingIndexByCenter(TempMainTubeIntersectionRings.Rings, StartCenterRing.Center);
		int32 EndIntersectionRingsStartIndex = FindRingIndexByCenter(TempMainTubeIntersectionRings.Rings, EndCenterRing.Center);

		OrderSquareIntersectionConnectionsOneCorner(
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempTopRightRingVertices,
			TempTopRightSquareVertices,
			EndIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.TopRightPartialRings,
			true
		);

		OrderSquareIntersectionConnectionsOneCorner(
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempBottomRightRingVertices,
			TempBottomRightSquareVertices,
			StartIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.BottomRightPartialRings,
			false
		);

		OrderSquareIntersectionConnectionsOneCorner(
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempBottomLeftRingVertices,
			TempBottomLeftSquareVertices,
			StartIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.BottomLeftPartialRings,
			false
		);

		OrderSquareIntersectionConnectionsOneCorner(
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempTopLeftRingVertices,
			TempTopLeftSquareVertices,
			EndIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.TopLeftPartialRings,
			true
		);
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                            7. Vertex and Index Utilities                                                             //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

	void FindClosestVertex(const FVector& InputVertex, const TArray<FVector>& Vertices, FVector& OutVertex, int32& OutIndex) {
		if (Vertices.Num() == 0) {
			UE_LOG(LogTemp, Error, TEXT("Vertices array is empty"));
			OutVertex = FVector::ZeroVector; // Return a default vertex
			OutIndex = -1; // Set invalid index
			return;
		}

		// init closest vertex and minimum distance
		OutIndex = 0; // Initialize to the first index
		FVector ClosestVertex = Vertices[0];
		float MinDistance = FVector::DistSquared(InputVertex, ClosestVertex);

		// loop through vertices to find the closest
		for (int32 i = 0; i < Vertices.Num(); ++i) {
			float Distance = FVector::DistSquared(InputVertex, Vertices[i]);
			if (Distance < MinDistance) {
				MinDistance = Distance;
				ClosestVertex = Vertices[i];
				OutIndex = i;
			}
		}

		OutVertex = ClosestVertex;
	}

	int32 FindVertexIndex(const TArray<FVector>& Vertices, const FVector& TargetVertex) {
		for (int32 i = 0; i < Vertices.Num(); i++) {
			if (Vertices[i] == TargetVertex) {
				return i;
			}
		}

		return -1;
	}

	void RemoveDuplicateVertices(TArray<FVector>& Vertices, float Tolerance) {
		TArray<FVector> UniqueVertices;

		for (const FVector& Vertex : Vertices) {
			bool bIsDuplicate = false;

			for (const FVector& UniqueVertex : UniqueVertices) {
				if (FMath::Abs(Vertex.X - UniqueVertex.X) < Tolerance &&
					FMath::Abs(Vertex.Y - UniqueVertex.Y) < Tolerance &&
					FMath::Abs(Vertex.Z - UniqueVertex.Z) < Tolerance) {
					bIsDuplicate = true;
					break;
				}
			}

			if (!bIsDuplicate) {
				UniqueVertices.Add(Vertex);
			}
		}

		Vertices = UniqueVertices;
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                        8. Transformations and Interpolations                                                         //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

	FVector CalculateCenterLinePoint(
		const FVector& Point,
		const FRingData& StartRing,
		const FRingData& EndRing
	) {
		// Calculate the direction vector of the centerline
		FVector CenterlineDirection = (EndRing.Center - StartRing.Center).GetSafeNormal();

		// Project the vector from StartRing.Center to the Point onto the centerline
		FVector VectorToPoint = Point - StartRing.Center;
		float ProjectedLength = FVector::DotProduct(VectorToPoint, CenterlineDirection);

		// Clamp the projection to lie between the start and end rings
		float FrustumLength = FVector::Dist(StartRing.Center, EndRing.Center);
		ProjectedLength = FMath::Clamp(ProjectedLength, 0.0f, FrustumLength);

		// Compute the point on the centerline corresponding to the projection
		FVector CenterlinePoint = StartRing.Center + (CenterlineDirection * ProjectedLength);

		// Return the centerline point as the result
		return CenterlinePoint;
	}

	float InterpolatedRingRadius(
		const FVector& CenterlinePoint,
		const FRingData& StartRing,
		const FRingData& EndRing
	) {
		// calc the distance between the 2 rings
		float TotalDistance = FVector::Dist(StartRing.Center, EndRing.Center);

		// calc the distance from start center to centerline point
		float DistanceToPoint = FVector::Dist(StartRing.Center, CenterlinePoint);

		// ensure not dividing by zero
		if (TotalDistance == 0.0f) {
			return (StartRing.Radius + EndRing.Radius) * 0.5f; // return avg if rings overlap
		}

		// calc the interpolation factor (0.0 at StartRingCenter, 1.0 at EndRingCenter)
		float InterpolationFactor = DistanceToPoint / TotalDistance;

		// linearly interpolate the radius based on the interpolation factor
		return FMath::Lerp(StartRing.Radius, EndRing.Radius, InterpolationFactor);
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                                9. Utility Functions                                                                  //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                            10. Geometry Query Functions                                                              //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //

	bool IsPointInsideFrustum(const FRingData& StartRing, const FRingData& EndRing, const FVector& Point) {
		// compute the height vector and its squared length
		FVector LengthVector = EndRing.Center - StartRing.Center;
		float LengthSquared = LengthVector.SizeSquared();

		// project the point onto the height vector
		FVector PointToStart = Point - StartRing.Center;
		float DotProduct = FVector::DotProduct(PointToStart, LengthVector);
		float t = DotProduct / LengthSquared;

		// check if the projection is within the frustum's height
		if (t < 0.0f || t > 1.0f) {
			return false;
		}

		// interpolate the radius at the projected height
		float InterpolatedRadius = FMath::Lerp(StartRing.Radius, EndRing.Radius, t);

		// compute the projected point on the frustum axis
		FVector ProjectedPoint = StartRing.Center + t * LengthVector;

		// Check the distance from the axis
		float DistanceToAxis = FVector::Dist(Point, ProjectedPoint);

		if (DistanceToAxis > InterpolatedRadius) {
			return false; // Preliminary check fails
		}

		// Final polygonal cross-section check
		const TArray<FVector>& CurrentRingVertices = (t < 0.5f) ? StartRing.Vertices : EndRing.Vertices;
		return IsPointInsidePolygon(CurrentRingVertices, Point, ProjectedPoint);
	}

	bool IsPointInsidePolygon(const TArray<FVector>& RingVertices, const FVector& Point, const FVector& RingCenter) {
		int NumVertices = RingVertices.Num();
		if (NumVertices < 3) return false; // Not a valid polygon

		// Compute the polygon's normal using the first three vertices
		FVector Edge1 = RingVertices[1] - RingVertices[0];
		FVector Edge2 = RingVertices[2] - RingVertices[0];
		FVector Normal = FVector::CrossProduct(Edge1, Edge2).GetSafeNormal();

		// Compute a local 2D basis (U, V) within the plane
		FVector U = Edge1.GetSafeNormal();
		FVector V = FVector::CrossProduct(Normal, U).GetSafeNormal();

		// Project the polygon vertices and the point into the local 2D plane
		TArray<FVector2D> ProjectedVertices;
		for (const FVector& Vertex : RingVertices) {
			FVector LocalVertex = Vertex - RingCenter;
			ProjectedVertices.Add(FVector2D(FVector::DotProduct(LocalVertex, U), FVector::DotProduct(LocalVertex, V)));
		}

		FVector LocalPoint = Point - RingCenter;
		FVector2D ProjectedPoint(FVector::DotProduct(LocalPoint, U), FVector::DotProduct(LocalPoint, V));

		// Perform the 2D point-in-polygon test
		bool bInside = false;
		for (int i = 0, j = NumVertices - 1; i < NumVertices; j = i++) {
			const FVector2D& Vertex1 = ProjectedVertices[i];
			const FVector2D& Vertex2 = ProjectedVertices[j];

			if (((Vertex1.Y > ProjectedPoint.Y) != (Vertex2.Y > ProjectedPoint.Y)) &&
				(ProjectedPoint.X < (Vertex2.X - Vertex1.X) * (ProjectedPoint.Y - Vertex1.Y) / (Vertex2.Y - Vertex1.Y) + Vertex1.X)) {
				bInside = !bInside;
			}
		}

		return bInside;
	}

	bool IsVertexInsideSquareAngle(const FVector& Vertex, const FIntersectionSquareData& IntersectionSquare) {
		// convert vertex to polar coords
		FVector2D Vertex2D(Vertex.X - IntersectionSquare.Center.X, Vertex.Y - IntersectionSquare.Center.Y);
		float VertexAngle = FMath::Atan2(Vertex2D.Y, Vertex2D.X);

		// sort square angles for easier range checking
		TArray<float> SortedAngles = IntersectionSquare.Angles;
		SortedAngles.Sort();

		// ensure the angles are in the range [0, 2 * PI]
		if (SortedAngles[0] < 0) {
			for (float& Angle : SortedAngles) {
				Angle += 2 * PI;
			}
		}

		// check if the vertex angle is within the square's angle range
		return VertexAngle >= SortedAngles[0] && VertexAngle <= SortedAngles[1];
	}

	bool IsPointOnLine(const FVector& LineStart, const FVector& LineEnd, const FVector& Point, float Tolerance) {
		// Check if the point is the same as the start or end points
		if (FVector::DistSquared(Point, LineStart) < FMath::Square(Tolerance) ||
			FVector::DistSquared(Point, LineEnd) < FMath::Square(Tolerance)) {
			return true;
		}

		// Calculate direction vectors
		FVector LineVector = LineEnd - LineStart;
		FVector PointVector = Point - LineStart;

		// Check if the point is collinear with the line
		FVector CrossProduct = FVector::CrossProduct(LineVector, PointVector);
		if (!CrossProduct.IsNearlyZero(Tolerance)) {
			return false; // Point is not on the line
		}

		// Check if the point lies between the start and end points
		float DotProduct = FVector::DotProduct(LineVector, PointVector);
		if (DotProduct < 0 || DotProduct > LineVector.SizeSquared()) {
			return false; // Point is outside the segment
		}

		return true; // Point is on the line segment
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //
	//                                                          11. Profiling and Curve Utilities                                                           //
	// ---------------------------------------------------------------------------------------------------------------------------------------------------- //	
}