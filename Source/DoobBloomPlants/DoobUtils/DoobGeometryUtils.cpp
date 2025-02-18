#include "DoobGeometryUtils.h"

#include "Math/Vector.h"
#include "Math/Transform.h"
#include "Containers/Array.h"
#include "DrawDebugHelpers.h"

//#include "DoobProfileUtils.h"
#include "DoobContainerUtils.h"
#include "DoobMathUtils.h"

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

		// Compute the centroid of the input vertices.
		FVector Centroid = ComputeCentroid(InputVertices);
		TArray<TPair<float, FVector>> VerticesWithAngles;
		int32 NumVertices = InputVertices.Num();

		// Compute an average normal for the polygon by summing cross products
		FVector AverageNormal = FVector::ZeroVector;
		for (int32 i = 0; i < NumVertices; i++)
		{
			FVector NextVertex = InputVertices[(i + 1) % NumVertices];
			FVector Edge1 = InputVertices[i] - Centroid;
			FVector Edge2 = NextVertex - Centroid;
			AverageNormal += FVector::CrossProduct(Edge1, Edge2);
		}
		AverageNormal.Normalize();

		// Establish a local 2D coordinate system in the plane of the polygon.
		// Choose a reference vector from the centroid to the first vertex.
		FVector Ref = (InputVertices[0] - Centroid).GetSafeNormal();
		// Compute a perpendicular vector in the plane using the average normal.
		FVector Perp = FVector::CrossProduct(AverageNormal, Ref).GetSafeNormal();

		// For each vertex, compute its angle relative to our local axes.
		for (const FVector& Vertex : InputVertices)
		{
			FVector Offset = Vertex - Centroid;
			float x = FVector::DotProduct(Offset, Ref);
			float y = FVector::DotProduct(Offset, Perp);
			float Angle = FMath::Atan2(y, x);
			VerticesWithAngles.Add(TPair<float, FVector>(Angle, Vertex));
		}

		// Sort vertices by angle (ascending)
		VerticesWithAngles.Sort([](const TPair<float, FVector>& A, const TPair<float, FVector>& B)
			{
				return A.Key < B.Key;
			});

		// Extract the sorted vertices.
		TArray<FVector> SortedVertices;
		for (const TPair<float, FVector>& Pair : VerticesWithAngles)
		{
			SortedVertices.Add(Pair.Value);
		}

		// Compute the signed area of the polygon in our 2D coordinate system.
		// Using the shoelace formula.
		float SignedArea = 0.0f;
		for (int32 i = 0; i < SortedVertices.Num(); i++)
		{
			int32 j = (i + 1) % SortedVertices.Num();
			FVector OffsetI = SortedVertices[i] - Centroid;
			FVector OffsetJ = SortedVertices[j] - Centroid;
			float x_i = FVector::DotProduct(OffsetI, Ref);
			float y_i = FVector::DotProduct(OffsetI, Perp);
			float x_j = FVector::DotProduct(OffsetJ, Ref);
			float y_j = FVector::DotProduct(OffsetJ, Perp);
			SignedArea += (x_i * y_j - x_j * y_i);
		}
		// In our coordinate system, a positive signed area means counter-clockwise.
		// If we need clockwise order, we want a negative area.
		if (SignedArea > 0)
		{
			Algo::Reverse(SortedVertices);
		}

		return SortedVertices;

		//FVector Centroid = ComputeCentroid(InputVertices);

		//// array to store vertices with angles
		//TArray<TPair<float, FVector>> VerticesWithAngles;

		//// Compute a better plane normal based on averaging cross products of centroid and vertices
		//FVector AverageNormal = FVector::ZeroVector;
		//for (int32 i = 0; i < InputVertices.Num(); i++) {
		//	FVector NextVertex = InputVertices[(i + 1) % InputVertices.Num()];
		//	FVector Edge1 = InputVertices[i] - Centroid;
		//	FVector Edge2 = NextVertex - Centroid;
		//	AverageNormal += FVector::CrossProduct(Edge1, Edge2);
		//}
		//AverageNormal = AverageNormal.GetSafeNormal(); // Normalize the average normal

		//// project vertices onto a plane
		//for (const FVector& Vertex : InputVertices) {
		//	FVector Offset = Vertex - Centroid;
		//	FVector Projected = Offset - FVector::DotProduct(Offset, AverageNormal) * AverageNormal;

		//	// compute angle in projected space
		//	float Angle = FMath::Atan2(Projected.Y, Projected.X);
		//	VerticesWithAngles.Add(TPair<float, FVector>(Angle, Vertex));
		//}

		//// sort vertices by angle
		//VerticesWithAngles.Sort([](const TPair<float, FVector>& A, const TPair<float, FVector>& B) {
		//	return A.Key < B.Key; // sort by angle
		//	});

		//// extract sorted vertices
		//TArray<FVector> SortedVertices;
		//for (const TPair<float, FVector>& Pair : VerticesWithAngles) {
		//	SortedVertices.Add(Pair.Value);
		//}

		//return SortedVertices;
	}

	TArray<FVector> ReorderRingVerticesToDirection(const TArray<FVector>& RingVertices, const FVector& InputDirection) {
		//TArray<FVector> OutRingVertices = RingVertices;

		int32 ClosestIndex = 0;
		float MaxProjection = -FLT_MAX;

		for (int32 i = 0; i < RingVertices.Num(); i++) {
			float projection = FVector::DotProduct(RingVertices[i], InputDirection);
			if (projection > MaxProjection) {
				MaxProjection = projection;
				ClosestIndex = i;
			}
		}

		return DoobContainerUtils::ReorderedArray(RingVertices, ClosestIndex);
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
		FindClosestVertex(StartCenter, IntersectionRing.CombinedVertices, LowestVertex, LowestIndex);
		FindClosestVertex(EndCenter, IntersectionRing.CombinedVertices, HighestVertex, HighestIndex);

		// Calculate the ring's center
		FVector RingCenter = ComputeCentroid(IntersectionRing.CombinedVertices);

		// Define the ring's axis (from StartCenter to EndCenter)
		FVector RingAxis = (EndCenter - StartCenter).GetSafeNormal();

		// new if not working ----------------------------------------------------------------------
		
		// Compute the tube's midpoint from the main tube start and end positions.
		FVector TubeMidpoint = (StartCenter + EndCenter) * 0.5f;

		// Compute the "inside vector": from the ring center toward the tube's midpoint.
		// This vector points inward (i.e. toward the inside of the tube).
		FVector InsideVector = (TubeMidpoint - RingCenter).GetSafeNormal();

		// Now, compute a "right reference" vector.
		// We want to define "right" in a way that is consistent no matter the lateral tube's sign.
		// Compute it as the cross product of the ring axis and the inside vector.
		FVector RightRef = FVector::CrossProduct(RingAxis, InsideVector).GetSafeNormal();

		// For safety, if RightRef is nearly zero (degenerate case), use an alternative.
		if (RightRef.IsNearlyZero())
		{
			RightRef = FVector::CrossProduct(RingAxis, FVector::RightVector).GetSafeNormal();
		}

		// At this point, a positive projection onto RightRef means the vertex is on the "right" side
		// (when viewing from outside the tube, where the inside is toward TubeMidpoint), 
		// and a negative projection means it is on the "left" side.

		// Find left and right points based on the dot product with RightRef.
		float MaxProjection = -FLT_MAX;
		float MinProjection = FLT_MAX;
		FVector RightVertex, LeftVertex;
		int32 RightIndex = -1, LeftIndex = -1;

		for (int32 i = 0; i < ArraySize; i++)
		{
			float Projection = FVector::DotProduct(IntersectionRing.CombinedVertices[i] - RingCenter, RightRef);

			if (Projection > MaxProjection)
			{
				MaxProjection = Projection;
				RightVertex = IntersectionRing.CombinedVertices[i];
				RightIndex = i;
			}
			if (Projection < MinProjection)
			{
				MinProjection = Projection;
				LeftVertex = IntersectionRing.CombinedVertices[i];
				LeftIndex = i;
			}
		}

		// old if not working ----------------------------------------------------------------------

		// Define a perpendicular axis
		//FVector UpVector = FVector::UpVector; // Can be replaced with another consistent vector
		//FVector PerpendicularAxis = FVector::CrossProduct(UpVector, RingAxis).GetSafeNormal();

		//// Ensure the PerpendicularAxis is valid
		//if (PerpendicularAxis.IsZero()) {
		//	PerpendicularAxis = FVector::CrossProduct(FVector::RightVector, RingAxis).GetSafeNormal();
		//}

		//// Find left and right points based on the perpendicular axis
		//float MaxProjection = -FLT_MAX;
		//float MinProjection = FLT_MAX;
		//FVector LeftVertex, RightVertex;
		//int32 LeftIndex = -1, RightIndex = -1;

		//for (int32 i = 0; i < ArraySize; i++) {
		//	float Projection = FVector::DotProduct(IntersectionRing.CombinedVertices[i] - RingCenter, PerpendicularAxis);

		//	if (Projection > MaxProjection) {
		//		MaxProjection = Projection;
		//		RightVertex = IntersectionRing.CombinedVertices[i];
		//		RightIndex = i;
		//	}

		//	if (Projection < MinProjection) {
		//		MinProjection = Projection;
		//		LeftVertex = IntersectionRing.CombinedVertices[i];
		//		LeftIndex = i;
		//	}
		//}

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
		/*bool Corner3Added = false;
		bool Corner23Added = false;*/

		int32 indexCount = 0;
		int32 indexCountFinal = 0;

		for (const FVector& Vertex : TubeIntersectionData.MTAboveIntersectionRing.Vertices) {
			FVector VertexDirection = (Vertex - TubeIntersectionData.MTAboveIntersectionRing.Center).GetSafeNormal();

			// check if the vertex direction is between the two interpolated directions
			FVector Normal = FVector::CrossProduct(TopDirections[0], TopDirections[1]).GetSafeNormal(); // Normal of the plane
			bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(TopDirections[0], VertexDirection)) >= 0 &&
				FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, TopDirections[1])) >= 0;

			if (!IsBetween) {
				if (!Corner2Added) {
					Corner2Added = true;
				}
				TempTopRing.Vertices.Add(Vertex);
				indexCount++;
			}
			//else if (Corner3Added) {
			//	//TempTopRing.Vertices.Add(IntersectionSquare.Corners[3]);
			//	Corner3Added = false;
			//	Corner2Added = true;
			//	//indexCount++;
			//}
			else if (Corner2Added) {
				//TempTopRing.Vertices.Add(IntersectionSquare.Corners[2]);
				Corner2Added = false;
				//Corner23Added = true;
				indexCountFinal = indexCount;
			}
		}

		TArray<FVector> ReorderedTopVertices;

		/*if (!Corner23Added) {
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
		else {*/
			ReorderedTopVertices = DoobContainerUtils::ReorderedArray(TempTopRing.Vertices, indexCountFinal);
		//}

		TempTopRing.Vertices = ReorderedTopVertices;

		FRingData TempBottomRing = TubeIntersectionData.MTBelowIntersectionRing;
		TempBottomRing.Vertices.Empty();

		indexCount = 0;
		indexCountFinal = 0;

		// Flags to track corner insertion
		//bool Corner0Added = false;
		bool Corner1Added = false;
		//bool Corner01Added = false;
		for (const FVector& Vertex : TubeIntersectionData.MTBelowIntersectionRing.Vertices) {
			FVector VertexDirection = (Vertex - TubeIntersectionData.MTBelowIntersectionRing.Center).GetSafeNormal();

			// check if the vertex direction is between the two interpolated directions
			FVector Normal = FVector::CrossProduct(BottomDirections[0], BottomDirections[1]).GetSafeNormal(); // Normal of the plane
			bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(BottomDirections[0], VertexDirection)) >= 0 &&
				FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, BottomDirections[1])) >= 0;

			if (!IsBetween) {
				if (!Corner1Added) {
					Corner1Added = true;
				}
				TempBottomRing.Vertices.Add(Vertex);
				indexCount++;
			}
			//else if (Corner0Added) {
			//	//TempBottomRing.Vertices.Add(IntersectionSquare.Corners[0]);
			//	Corner0Added = false;
			//	Corner1Added = true;
			//	//indexCount++;
			//}
			else if (Corner1Added) {
				//TempBottomRing.Vertices.Add(IntersectionSquare.Corners[1]);
				Corner1Added = false;
				//Corner01Added = true;
				indexCountFinal = indexCount;
			}
		}

		TArray<FVector> ReorderedBottomVertices;

		/*if (!Corner01Added) {
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
		else {*/
			ReorderedBottomVertices = DoobContainerUtils::ReorderedArray(TempBottomRing.Vertices, indexCountFinal);
		//}

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

			//bool LeftRingConnection = false;
			bool RightRingConnection = false;
			//bool BothConnection = false;

			// check and remove vertices
			for (const FVector& Vertex : Ring.Vertices) {
				FVector VertexDirection = (Vertex - Ring.Center).GetSafeNormal();

				// check if the vertex direction is between the two interpolated directions
				FVector Normal = FVector::CrossProduct(InterpolatedDirection1, InterpolatedDirection2).GetSafeNormal(); // Normal of the plane
				bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(InterpolatedDirection1, VertexDirection)) >= 0 &&
					FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, InterpolatedDirection2)) >= 0;

				if (!IsBetween) {
					if (!RightRingConnection) {
						RightRingConnection = true;
					}
					ModifiedRing.Vertices.Add(Vertex);
					indexCount++;
				}
				//else if (LeftRingConnection) {
				//	//FVector LeftRingConnectionVertex = FindIntersectionOnRing(Ring.Vertices, InterpolatedDirection1, Ring.Center);
				//	//ModifiedRing.Vertices.Add(LeftRingConnectionVertex);
				//	LeftRingConnection = false;
				//	RightRingConnection = true;
				//	//indexCount++;
				//}
				else if (RightRingConnection) {
					//FVector RightRingConnectionVertex = FindIntersectionOnRing(Ring.Vertices, InterpolatedDirection2, Ring.Center);
					//ModifiedRing.Vertices.Add(RightRingConnectionVertex);
					RightRingConnection = false;
					//BothConnection = true;
					indexCountFinal = indexCount;
				}
			}

			TArray<FVector> ReorderedModifiedVertices;

			/*if (!BothConnection) {
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
			else {*/
				ReorderedModifiedVertices = DoobContainerUtils::ReorderedArray(ModifiedRing.Vertices, indexCountFinal);
			//}

			ModifiedRing.Vertices = ReorderedModifiedVertices;
			ModifiedRing.bIsClosed = false;
			ModifiedRing.bIsComplete = false;

			TempTube.Add(ModifiedRing);
		}

		TubeIntersectionData.MainTube.Rings = TempTube;
	}

	bool IsPointOnRing(const FRingData& Ring, const FVector& TestVertex) {
		int32 VertNum = Ring.Vertices.Num();
		for (int32 i = 0; i < Ring.Vertices.Num(); i++) {
			FVector LineStart = Ring.Vertices[i];
			FVector LineEnd = Ring.Vertices[(i + 1) % VertNum];

			if (IsPointOnLine(LineStart, LineEnd, TestVertex)) return true;
		}

		return false;
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

	void RemoveInternalVertices(const FTubeData& TubeA, FTubeData& TubeB, FTubeData& OutTube) {
		TArray<FRingData> TempTube;

		// Iterate through all rings in TubeB
		for (int32 TubeRingIndexB = 0; TubeRingIndexB < TubeB.Rings.Num(); ++TubeRingIndexB) {
			FRingData CurrentRingB = TubeB.Rings[TubeRingIndexB];
			/*int32 NextRingBIndex = (TubeRingIndexB + 1) % TubeB.Rings.Num();
			FRingData NextRingB = TubeB.Rings[NextRingBIndex];*/
			
			FRingData TempRing;

			FEdgeRingData TempEdgeRing;
			TempEdgeRing.RingAIndex = TubeRingIndexB;
			TempEdgeRing.RingBIndex = TubeRingIndexB;

			bool bIntersectionFound = false;

			// Check each vertex in CurrentRingB against all frustums in TubeA
			for (int32 VertexIndex = 0; VertexIndex < CurrentRingB.Vertices.Num(); ++VertexIndex) {
				FVector CurrentVertex = CurrentRingB.Vertices[VertexIndex];
				int32 NextVertexIndex = (VertexIndex + 1) % CurrentRingB.Vertices.Num();
				FVector NextVertex = CurrentRingB.Vertices[NextVertexIndex];

				bool bInsideAnyFrustum = false;
				bool bNextInsideAnyFrustum = false;
				bool bIntersection = false;

				// Test against all frustums defined by TubeA's consecutive rings
				for (int32 TubeRingIndexA = 0; TubeRingIndexA < TubeA.Rings.Num() - 1; ++TubeRingIndexA) {
					const FRingData& CurrentRingA = TubeA.Rings[TubeRingIndexA];
					const FRingData& NextRingA = TubeA.Rings[TubeRingIndexA + 1];

					if (!bInsideAnyFrustum && IsPointInsideFrustum(CurrentRingA, NextRingA, CurrentVertex)) {
						bInsideAnyFrustum = true;
						bIntersection = true;
						//break; // No need to check further if inside any frustum
					}

					if (!bNextInsideAnyFrustum && IsPointInsideFrustum(CurrentRingA, NextRingA, NextVertex)) {
						bNextInsideAnyFrustum = true;
						bIntersection = true;
					}

					if (bInsideAnyFrustum && bNextInsideAnyFrustum) {
						bIntersection = false;
						break; // No need to check further if inside any frustum
					}
				}

				// Only keep vertices that are outside all frustums
				if (!bInsideAnyFrustum) {
					TempRing.Vertices.Add(CurrentVertex);
				}

				FEdgeData TempEdge;

				TempEdge.VertexAIndex = VertexIndex;
				TempEdge.VertexBIndex = NextVertexIndex;
				TempEdge.Length = FVector::Dist(CurrentVertex, NextVertex);
				TempEdge.VertexARemoved = bInsideAnyFrustum;
				TempEdge.VertexBRemoved = bNextInsideAnyFrustum;

				if (bIntersection && !bIntersectionFound) {
					bIntersectionFound = true;
				}

				TempEdgeRing.Edges.Add(TempEdge);
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

			TempEdgeRing.IntersectingEdge = bIntersectionFound;

			TubeB.RingEdges.Add(TempEdgeRing);
		}

		OutTube.Rings = TempTube;
	}

	void FindSegmentForPoint(
		const FVector& Point,
		const FTubeData& Tube,
		FRingData& StartRing,
		FRingData& EndRing
	) {
		// Initialize with fallback values
		/*StartRing = FRingData();
		EndRing = FRingData();*/

		TArray<TArray<FRingData>> RingPairs;

		// iterate through the rings to find the segments the point is between
		for (int32 i = 0; i < Tube.Rings.Num() - 1; ++i) {
			// get the current ring and the next ring
			FRingData CurrentRing = Tube.Rings[i];
			FRingData NextRing = Tube.Rings[i + 1];

			// Skip pairs with empty vertices
			if (CurrentRing.Vertices.IsEmpty() || NextRing.Vertices.IsEmpty())
			{
				UE_LOG(LogTemp, Warning, TEXT("Skipping ring pair %d-%d (empty vertices)"), i, i + 1);
				continue;
			}

			// calc the distance between the 2 rings
			float SegmentLength = FVector::Dist(CurrentRing.Center, NextRing.Center);

			// calc the length of the point relative to current ring - for interpolation
			FVector PointToCurrentRing = Point - CurrentRing.Center;
			float PointLength = FVector::DotProduct(PointToCurrentRing, (NextRing.Center - CurrentRing.Center).GetSafeNormal());

			// check if point is within segments
			if (PointLength >= 0 && PointLength <= SegmentLength) {
				TArray<FRingData> RingPair;
				RingPair.Add(CurrentRing);
				RingPair.Add(NextRing);
				RingPairs.Add(RingPair);
				StartRing = CurrentRing;
				EndRing = NextRing;
			}
		}

		float MinDistance = FLT_MAX;

		for (int32 i = 0; i < RingPairs.Num(); i++) {
			float DistanceToStart = FVector::DistSquared(Point, RingPairs[i][0].Center);
			float DistanceToEnd = FVector::DistSquared(Point, RingPairs[i][1].Center);

			float Distance;
			if (DistanceToStart > DistanceToEnd) {
				Distance = DistanceToStart;
			}
			else {
				Distance = DistanceToEnd;
			}

			if (Distance < MinDistance) {
				MinDistance = Distance;
				StartRing = RingPairs[i][0];
				EndRing = RingPairs[i][1];
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

	FTubeData ReorderTubeVerticesToDirection(const FTubeData& InputTube, const FVector& InputDirection) {
		FTubeData OutTube = InputTube;
		OutTube.Rings.Empty();

		for (int32 i = 0; i < InputTube.Rings.Num(); ++i) {
			FRingData TempRing = InputTube.Rings[i];
			TempRing.Vertices = ReorderRingVerticesToDirection(InputTube.Rings[i].Vertices, InputDirection);
			OutTube.Rings.Add(TempRing);
		}

		return OutTube;
	}

	void GenerateTubeEdgesQuadrilaterals(FTubeData& Tube) {
		UE_LOG(LogTemp, Log, TEXT("Generate Quads and Edges"))
		int32 NumRings = Tube.RingEdges.Num();
		for (int32 i = 0; i < NumRings - 1; i++) {
			UE_LOG(LogTemp, Log, TEXT("Ring num: %d"), i);
			FEdgeRingData CurrentEdgeRing = Tube.RingEdges[i];
			int32 NextIndex = (i + 1) % NumRings;
			FEdgeRingData NextEdgeRing = Tube.RingEdges[NextIndex];

			FEdgeRingData TempEdgeRing;
			TempEdgeRing.RingAIndex = i;
			TempEdgeRing.RingBIndex = NextIndex;

			FQuadrilateralRingData TempQuadRing;
			TempQuadRing.RingAIndex = i;
			TempQuadRing.RingBIndex = NextIndex;
			TempQuadRing.Height = FVector::Dist(
				Tube.Rings[i].Center, 
				Tube.Rings[NextIndex].Center
			);

			bool bIntersection = false;

			for (int32 j = 0; j < CurrentEdgeRing.Edges.Num(); j++) {
				FEdgeData EdgeA;
				FEdgeData EdgeB;
				FQuadrilateralData Quadrilateral;

				int32 NextVertexIndex = (j + 1) % CurrentEdgeRing.Edges.Num();

				EdgeA.VertexAIndex = CurrentEdgeRing.Edges[j].VertexAIndex;
				EdgeA.VertexBIndex = NextEdgeRing.Edges[j].VertexAIndex;
				EdgeA.VertexARemoved = CurrentEdgeRing.Edges[j].VertexARemoved;
				EdgeA.VertexBRemoved = NextEdgeRing.Edges[j].VertexARemoved;
				EdgeA.Length = FVector::Dist(
					Tube.Rings[i].Vertices[j],
					Tube.Rings[NextIndex].Vertices[j]
				);

				bool bForwardIntersection = EdgeA.VertexARemoved && !EdgeA.VertexBRemoved;
				bool bBackwardIntersection = !EdgeA.VertexARemoved && EdgeA.VertexBRemoved;
				if (bBackwardIntersection || bForwardIntersection) {
					bIntersection = true;
				}

				EdgeB.VertexAIndex = CurrentEdgeRing.Edges[j].VertexBIndex;
				EdgeB.VertexBIndex = NextEdgeRing.Edges[j].VertexBIndex;
				EdgeB.VertexARemoved = CurrentEdgeRing.Edges[j].VertexBRemoved;
				EdgeB.VertexBRemoved = NextEdgeRing.Edges[j].VertexBRemoved;
				EdgeB.Length = FVector::Dist(
					Tube.Rings[i].Vertices[NextVertexIndex],
					Tube.Rings[NextIndex].Vertices[NextVertexIndex]
				);

				Quadrilateral.EdgeA = EdgeA;
				Quadrilateral.EdgeB = EdgeB;
				Quadrilateral.BottomWidth = FVector::Dist(
					Tube.Rings[i].Vertices[j],
					Tube.Rings[i].Vertices[NextVertexIndex]
				);
				Quadrilateral.TopWidth = FVector::Dist(
					Tube.Rings[NextIndex].Vertices[j],
					Tube.Rings[NextIndex].Vertices[NextVertexIndex]
				);

				TempEdgeRing.Edges.Add(EdgeA);
				TempQuadRing.Quadrilaterals.Add(Quadrilateral);
			}

			TempEdgeRing.IntersectingEdge = bIntersection; // added checks per ring need to utileze in genterating intersection ring next <-------------------------------------

			Tube.SideEdges.Add(TempEdgeRing);
			Tube.Quadrilaterals.Add(TempQuadRing);
		}
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
		TArray<FVector> MainTubeRingIntersectionVerts;
		TArray<FVector> LateralTubeRingIntersectionVerts;

		GenerateHalfIntersectionRing(MainTube, LateralTube, OutRingData.MainTubeVertices);
		GenerateHalfIntersectionRing(LateralTube, MainTube, OutRingData.LateralTubeVertices);

		FVector MainTubeReorderDirection = -DoobMathUtils::GetPerpendicularDirection(MainTube.Direction, LateralTube.Direction);
		FTubeData MainTubeReordered = ReorderTubeVerticesToDirection(MainTube, MainTubeReorderDirection);

		FVector LateralTubeReorderDirection = -DoobMathUtils::GetPerpendicularDirection(LateralTube.Direction, MainTube.Direction);
		FTubeData LateralTubeReordered = ReorderTubeVerticesToDirection(LateralTube, LateralTubeReorderDirection);

		GenerateHalfIntersectionRing(MainTube, LateralTube, MainTubeRingIntersectionVerts, true);
		OutRingData.MainTubeVertices.Append(MainTubeRingIntersectionVerts);

		GenerateHalfIntersectionRing(LateralTube, MainTube, LateralTubeRingIntersectionVerts, true);
		OutRingData.LateralTubeVertices.Append(LateralTubeRingIntersectionVerts);

		CombinedVertices = OutRingData.MainTubeVertices;
		CombinedVertices.Append(OutRingData.LateralTubeVertices);

		RemoveDuplicateVertices(CombinedVertices);

		MergeClosePoints(CombinedVertices, OutRingData.MinDistanceBetweenVertices);

		OutRingData.CombinedVertices = OrderRingVertices(CombinedVertices);

		OutRingData.Centroid = ComputeCentroid(OutRingData.CombinedVertices);
	}

	void GenerateHalfIntersectionRing(
		const FTubeData& TubeA,
		const FTubeData& TubeB,
		TArray<FVector>& OutRingVertices,
		bool bRingEdges,
		/*bool bForwardAndReverse,
		bool bOnePerRing,*/
		float Precision
	) {
		int32 TubeANumRings/* = TubeA.SideEdges.Num()*/;

		if (bRingEdges) {
			TubeANumRings = TubeA.RingEdges.Num();
		}
		else {
			TubeANumRings = TubeA.SideEdges.Num();
		}

		UE_LOG(LogTemp, Log, TEXT("Test Intersection Ring Logic"))

		for (int32 EdgeRingIndexA = 0; EdgeRingIndexA < TubeANumRings; EdgeRingIndexA++) {
			UE_LOG(LogTemp, Log, TEXT("Edge Ring num: %d"), EdgeRingIndexA);

			FEdgeRingData CurrentEdgeRing/* = TubeA.SideEdges[EdgeRingIndexA]*/;

			if (bRingEdges) {
				CurrentEdgeRing = TubeA.RingEdges[EdgeRingIndexA];
			}
			else {
				CurrentEdgeRing = TubeA.SideEdges[EdgeRingIndexA];
			}

			if (!CurrentEdgeRing.IntersectingEdge) {
				continue;
			}

			int32 NumEdges = CurrentEdgeRing.Edges.Num();

			for (int32 EdgeIndexA = 0; EdgeIndexA < NumEdges; EdgeIndexA++) {
				FEdgeData CurrentEdge = CurrentEdgeRing.Edges[EdgeIndexA];
				bool bBothRemoved = CurrentEdge.VertexARemoved && CurrentEdge.VertexBRemoved;
				bool bNeitherRemoved = !CurrentEdge.VertexARemoved && !CurrentEdge.VertexBRemoved;

				FVector LineStartA;
				FVector LineEndA;

				if (bBothRemoved) {
					UE_LOG(LogTemp, Log, TEXT("No Intersection both removed"));
					continue;
				}
				/*else if (bNeitherRemoved) {

				}*/
				else if (CurrentEdge.VertexARemoved) {
					LineStartA = TubeA.Rings[CurrentEdgeRing.RingBIndex].Vertices[CurrentEdge.VertexBIndex];
					LineEndA = TubeA.Rings[CurrentEdgeRing.RingAIndex].Vertices[CurrentEdge.VertexAIndex];
				}
				else {
					LineStartA = TubeA.Rings[CurrentEdgeRing.RingAIndex].Vertices[CurrentEdge.VertexAIndex];
					LineEndA = TubeA.Rings[CurrentEdgeRing.RingBIndex].Vertices[CurrentEdge.VertexBIndex];
				}

				UE_LOG(LogTemp, Log, TEXT("LineStartA: %s, LineEndA; %s"), *LineStartA.ToString(), *LineEndA.ToString());

				int32 TubeBNumRings = TubeB.Quadrilaterals.Num();

				for (int32 QuadRingIndexB = 0; QuadRingIndexB < TubeBNumRings; QuadRingIndexB++) {
					FQuadrilateralRingData CurrentQuadRing = TubeB.Quadrilaterals[QuadRingIndexB];

					int32 NumQuads = CurrentQuadRing.Quadrilaterals.Num();

					for (int32 QuadIndexB = 0; QuadIndexB < NumQuads; QuadIndexB++) {
						FQuadrilateralData CurrentQuad = CurrentQuadRing.Quadrilaterals[QuadIndexB];
						FVector V0 = TubeB.Rings[CurrentQuadRing.RingAIndex].Vertices[CurrentQuad.EdgeA.VertexAIndex];
						FVector V1 = TubeB.Rings[CurrentQuadRing.RingBIndex].Vertices[CurrentQuad.EdgeA.VertexBIndex];
						FVector V2 = TubeB.Rings[CurrentQuadRing.RingBIndex].Vertices[CurrentQuad.EdgeB.VertexBIndex];
						FVector V3 = TubeB.Rings[CurrentQuadRing.RingAIndex].Vertices[CurrentQuad.EdgeB.VertexAIndex];

						/*if (bNeitherRemoved) {
							FVector
						}*/

						// check for intersection with triangles
						FVector IntersectionPoint;

						bool DoesIntersect = LineSegmentIntersectsQuadrilateral(LineStartA, LineEndA, V0, V1, V2, V3, IntersectionPoint);

						if (DoesIntersect) {
							OutRingVertices.Add(IntersectionPoint);
						}
					}
				}
			}
		}

		// old -- remove when done

		//int32 TubeANumRings = TubeA.Rings.Num();

		//for (int32 LineIndexA = 0; LineIndexA < TubeANumRings - 1; ++LineIndexA) {
		//	int32 CurrentRingAIndex;
		//	int32 NextRingAIndex;

		//	if (bReverse) {
		//		CurrentRingAIndex = LineIndexA;
		//		NextRingAIndex = (LineIndexA + 1) % TubeANumRings;
		//	}
		//	else {
		//		CurrentRingAIndex = (TubeANumRings - 1 - LineIndexA + TubeANumRings) % TubeANumRings;
		//		NextRingAIndex = (TubeANumRings - 1 - (LineIndexA + 1) + TubeANumRings) % TubeANumRings;
		//	}

		//	const FRingData& CurrentRingA = TubeA.Rings[CurrentRingAIndex];
		//	const FRingData& NextRingA = TubeA.Rings[NextRingAIndex];

		//	int32 NumVerticesA = CurrentRingA.Vertices.Num();

		//	for (int32 VertexIndexA = 0; VertexIndexA < NumVerticesA; ++VertexIndexA) {
		//		FVector LineStartA = CurrentRingA.Vertices[VertexIndexA];
		//		FVector LineEndA = NextRingA.Vertices[VertexIndexA];

		//		FRingData StartRingA;
		//		FRingData EndRingA;
		//		FRingData StartRingB;
		//		FRingData EndRingB;

		//		FindSegmentForPoint(LineStartA, TubeB, StartRingA, EndRingA);
		//		FindSegmentForPoint(LineEndA, TubeB, StartRingB, EndRingB);

		//		TArray<FRingData> TempTubeB;
		//		TempTubeB.Add(StartRingA);
		//		TempTubeB.Add(EndRingA);
		//		TempTubeB.Add(StartRingB);
		//		TempTubeB.Add(EndRingB);

		//		// loop through TubeB rectangles
		//		for (int32 RingIndexB = 0; RingIndexB < TempTubeB.Num() - 1; ++RingIndexB) {
		//			const FRingData& CurrentRingB = TempTubeB[RingIndexB];
		//			const FRingData& NextRingB = TempTubeB[RingIndexB + 1];
		//			int32 NumVerticesB = CurrentRingB.Vertices.Num();

		//			for (int32 VertexIndexB = 0; VertexIndexB < NumVerticesB; ++VertexIndexB) {
		//				/*FVector V0 = CurrentRingB.Vertices[VertexIndexB];
		//				FVector V1 = CurrentRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];
		//				FVector V2 = NextRingB.Vertices[VertexIndexB];
		//				FVector V3 = NextRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];*/

		//				FVector V0 = CurrentRingB.Vertices[VertexIndexB];
		//				FVector V1 = NextRingB.Vertices[VertexIndexB];
		//				FVector V2 = NextRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];
		//				FVector V3 = CurrentRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];

		//				// check for intersection with triangles
		//				FVector IntersectionPoint;

		//				bool DoesIntersect = LineSegmentIntersectsQuadrilateral(LineStartA, LineEndA, V0, V1, V2, V3, IntersectionPoint);

		//				if (/*LineSegmentIntersectsTriangle(LineStartA, LineEndA, V0, V1, V2, IntersectionPoint) ||
		//					LineSegmentIntersectsTriangle(LineStartA, LineEndA, V1, V2, V3, IntersectionPoint)*/DoesIntersect) {
		//					OutRingVertices.Add(IntersectionPoint);
		//				}
		//			}
		//		}
		//	}
		//}
	}

	void GenerateHalfIntersectionRingUsingCircumference(
		const FTubeData& TubeA,
		const FTubeData& TubeB,
		TArray<FVector>& OutRingVertices,
		bool OnePerRing,
		float Precision
	) {
		for (int32 LineIndexA = 0; LineIndexA < TubeA.Rings.Num(); ++LineIndexA) {
			const FRingData& CurrentRingA = TubeA.Rings[LineIndexA];

			int32 NumVerticesA = CurrentRingA.Vertices.Num();

			bool ForwardIntersectionAdded = false;
			bool ReverseIntersectionAdded = false;

			//UE_LOG(LogTemp, Warning, TEXT("Processing Ring %d with %d vertices"), LineIndexA, NumVerticesA);

			for (int32 VertexIndexA = 0; VertexIndexA < NumVerticesA; ++VertexIndexA) {
				if (ForwardIntersectionAdded && ReverseIntersectionAdded && OnePerRing) continue;

				FVector LineStartA = CurrentRingA.Vertices[VertexIndexA];
				FVector LineEndA = CurrentRingA.Vertices[(VertexIndexA + 1) % NumVerticesA];

				FVector RevLineStartA = CurrentRingA.Vertices[(NumVerticesA - 1 - VertexIndexA + NumVerticesA) % NumVerticesA];
				FVector RevLineEndA = CurrentRingA.Vertices[(NumVerticesA - 1 - (VertexIndexA + 1) + NumVerticesA) % NumVerticesA];

				/*UE_LOG(LogTemp, Warning, TEXT("Forward Line: Start (%.3f, %.3f, %.3f) -> End (%.3f, %.3f, %.3f)"),
					LineStartA.X, LineStartA.Y, LineStartA.Z,
					LineEndA.X, LineEndA.Y, LineEndA.Z);

				UE_LOG(LogTemp, Warning, TEXT("Reverse Line: Start (%.3f, %.3f, %.3f) -> End (%.3f, %.3f, %.3f)"),
					RevLineStartA.X, RevLineStartA.Y, RevLineStartA.Z,
					RevLineEndA.X, RevLineEndA.Y, RevLineEndA.Z);*/

				FRingData StartRingA;
				FRingData EndRingA;
				FRingData StartRingB;
				FRingData EndRingB;

				FRingData RevStartRingA;
				FRingData RevEndRingA;
				FRingData RevStartRingB;
				FRingData RevEndRingB;

				FindSegmentForPoint(LineStartA, TubeB, StartRingA, EndRingA);
				FindSegmentForPoint(LineEndA, TubeB, StartRingB, EndRingB);

				FindSegmentForPoint(RevLineStartA, TubeB, RevStartRingA, RevEndRingA);
				FindSegmentForPoint(RevLineEndA, TubeB, RevStartRingB, RevEndRingB);

				//UE_LOG(LogTemp, Warning, TEXT("Found segments for forward and reverse line"));

				TArray<FRingData> TempTubeB;
				TempTubeB.Add(StartRingA);
				TempTubeB.Add(EndRingA);
				TempTubeB.Add(StartRingB);
				TempTubeB.Add(EndRingB);

				TempTubeB.Add(RevStartRingA);
				TempTubeB.Add(RevEndRingA);
				TempTubeB.Add(RevStartRingB);
				TempTubeB.Add(RevEndRingB);

				

				// loop through TubeB rectangles
				for (int32 RingIndexB = 0; RingIndexB < TempTubeB.Num() - 1; ++RingIndexB) {
					if (ForwardIntersectionAdded && ReverseIntersectionAdded && OnePerRing) continue;

					const FRingData& CurrentRingB = TempTubeB[RingIndexB];
					const FRingData& NextRingB = TempTubeB[RingIndexB + 1];
					int32 NumVerticesB = CurrentRingB.Vertices.Num();

					/*UE_LOG(LogTemp, Warning, TEXT("Processing TempTubeB Rings %d and %d with %d vertices"),
						RingIndexB, RingIndexB + 1, NumVerticesB);*/

					for (int32 VertexIndexB = 0; VertexIndexB < NumVerticesB; ++VertexIndexB) {
						if (ForwardIntersectionAdded && ReverseIntersectionAdded && OnePerRing) continue;

						FVector V0 = CurrentRingB.Vertices[VertexIndexB];
						FVector V1 = NextRingB.Vertices[VertexIndexB];
						FVector V2 = NextRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];
						FVector V3 = CurrentRingB.Vertices[(VertexIndexB + 1) % NumVerticesB];

						// check for intersection with triangles
						FVector IntersectionPoint;
						bool DoesForwardIntersect = LineSegmentIntersectsQuadrilateral(LineStartA, LineEndA, V0, V1, V2, V3, IntersectionPoint);
							/*LineSegmentIntersectsTriangle(LineStartA, LineEndA, V0, V1, V2, IntersectionPoint) ||
							LineSegmentIntersectsTriangle(LineStartA, LineEndA, V1, V2, V3, IntersectionPoint);*/

						FVector RevIntersectionPoint;
						bool DoesReverseIntersect = LineSegmentIntersectsQuadrilateral(RevLineStartA, RevLineEndA, V0, V1, V2, V3, RevIntersectionPoint);
							/*LineSegmentIntersectsTriangle(RevLineStartA, RevLineEndA, V0, V1, V2, RevIntersectionPoint) ||
							LineSegmentIntersectsTriangle(RevLineStartA, RevLineEndA, V1, V2, V3, RevIntersectionPoint);*/
						if (OnePerRing) {
							if (DoesForwardIntersect && !ForwardIntersectionAdded) {
								/*UE_LOG(LogTemp, Warning, TEXT("Forward intersection found at (%.3f, %.3f, %.3f)"),
									IntersectionPoint.X, IntersectionPoint.Y, IntersectionPoint.Z);*/
								OutRingVertices.Add(IntersectionPoint);
								ForwardIntersectionAdded = true;
							}

							if (DoesReverseIntersect && !ReverseIntersectionAdded) {
								/*UE_LOG(LogTemp, Warning, TEXT("Reverse intersection found at (%.3f, %.3f, %.3f)"),
									RevIntersectionPoint.X, RevIntersectionPoint.Y, RevIntersectionPoint.Z);*/
								OutRingVertices.Add(RevIntersectionPoint);
								ReverseIntersectionAdded = true;
							}
						}
						else {
							if (DoesForwardIntersect) {
								/*UE_LOG(LogTemp, Warning, TEXT("Forward intersection found at (%.3f, %.3f, %.3f)"),
									IntersectionPoint.X, IntersectionPoint.Y, IntersectionPoint.Z);*/
								OutRingVertices.Add(IntersectionPoint);
								//ForwardIntersectionAdded = true;
							}

							if (DoesReverseIntersect) {
								/*UE_LOG(LogTemp, Warning, TEXT("Reverse intersection found at (%.3f, %.3f, %.3f)"),
									RevIntersectionPoint.X, RevIntersectionPoint.Y, RevIntersectionPoint.Z);*/
								OutRingVertices.Add(RevIntersectionPoint);
								//ReverseIntersectionAdded = true;
							}
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
		FVector BottomLeft = FindIntersectionOnNewRing(LowestRing.Vertices, LeftDirection, LowestRing.Center, LeftVertex);
		FVector BottomRight = FindIntersectionOnNewRing(LowestRing.Vertices, RightDirection, LowestRing.Center, RightVertex);
		FVector TopLeft = FindIntersectionOnNewRing(HighestRing.Vertices, LeftDirection, HighestRing.Center, LeftVertex);
		FVector TopRight = FindIntersectionOnNewRing(HighestRing.Vertices, RightDirection, HighestRing.Center, RightVertex);

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

		// need to check if any vertexes are in the lateral tube, and add vertices where these ring pass into the lateral tube.


	}

	void GenerateLateralTubeIntersectionRings(FTwoTubeIntersectionData& TubeIntersectionData) {
		TubeIntersectionData.LateralTubeIntersectionRings = TubeIntersectionData.LateralTube;
		TubeIntersectionData.LateralTubeIntersectionRings.Rings.Empty();

		// create temp ring to reorder ring to north point
		TArray<FVector> ReorderedIntersection = DoobContainerUtils::ReorderedArray(
			TubeIntersectionData.IntersectionRing.CombinedVertices,
			TubeIntersectionData.IntersectionRing.CardinalIndices[0]
		);

		TArray<FVector> CenterlineVertices;

		for (int32 i = 0; i < ReorderedIntersection.Num(); i++) {
			FRingData CurrentStartRing;
			FRingData CurrentEndRing;
			FindSegmentForPoint(ReorderedIntersection[i], TubeIntersectionData.LateralTube, CurrentStartRing, CurrentEndRing);

			FVector CurrentVertexCenterPoint = CalculateCenterLinePoint(ReorderedIntersection[i], CurrentStartRing, CurrentEndRing);
			CenterlineVertices.Add(CurrentVertexCenterPoint);
		}

		// log the intersection centerpoints here

		UE_LOG(LogTemp, Log, TEXT("Centerline Vertices Created"));
		for (int32 i = 0; i < CenterlineVertices.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("i: %d, Vector: %s"), i, *CenterlineVertices[i].ToString());
		}

		// < -------------------------- ************************************************************************************************************************************************
		// < -------------------------- ************************************************************************************************************************************************
		// < -------------------------- need to use these centerline points to create the lateral tube intersection rings. they are now in a proper order, and any duplicates have been removed
		// < -------------------------- also now all vertices in the intersection ring will have an accompanying lateral intersection ring
		// < -------------------------- ************************************************************************************************************************************************
		RemoveDuplicateVertices(CenterlineVertices);

		// log the intersection centerpoints here
		UE_LOG(LogTemp, Log, TEXT("Centerline Vertices Dupes remove"));
		for (int32 i = 0; i < CenterlineVertices.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("i: %d, Vector: %s"), i, *CenterlineVertices[i].ToString());
		}

		OrderVerticesByDistanceToPoint(CenterlineVertices, TubeIntersectionData.LateralTube.StartPosition);

		UE_LOG(LogTemp, Log, TEXT("Centerline Vertices Ordered"));
		for (int32 i = 0; i < CenterlineVertices.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("i: %d, Vector: %s"), i, *CenterlineVertices[i].ToString());
		}

		// log the intersection centerpoints here

		//int32 StartIndex = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[1]) - 1;
		//int32 EndIndexAlt = FindVertexIndex(ReorderedIntersection, TubeIntersectionData.IntersectionRing.CardinalVertices[3]) + 1;

		/*TArray<FVector> NorthSouthVerts = {
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
		}*/

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

		FirstLateralIntersectionRing.Vertices = ReorderedIntersection;

		/*for (int32 i = StartIndex; i <= EndIndexAlt; ++i) {
			FVector VertexToAdd = ReorderedIntersection[i];

			FirstLateralIntersectionRing.Vertices.Add(VertexToAdd);

			FVector SuperTempCenter = CalculateCenterLinePoint(VertexToAdd, IntersectionRingStartRing, IntersectionRingEndRing);
			FVector SuperTempDirection = (VertexToAdd - SuperTempCenter).GetSafeNormal();
		}*/

		TubeIntersectionData.LateralTubeIntersectionRings.Rings.Add(FirstLateralIntersectionRing);

		int32 EndIndex = 0;

		FRingData CurrentLateralIntersectionRing;
		FRingData PreviousLateralIntersectionRing;

		bool bPointsIntersect = false;

		for (int32 i = /*StartIndex*/0; i < CenterlineVertices.Num()/*>= EndIndex*/; /*--i*/++i) {
			bPointsIntersect = false;

			FRingData CurrentStartRing;
			FRingData CurrentEndRing;

			TArray<FVector> TempCurrentVertices;
			TArray<FVector> TempPreviousVertices;

			/*FindSegmentForPoint(ReorderedIntersection[i], TubeIntersectionData.LateralTube, CurrentStartRing, CurrentEndRing);

			FVector CurrentVertexCenterPoint = CalculateCenterLinePoint(ReorderedIntersection[i], CurrentStartRing, CurrentEndRing);*/
			FindSegmentForPoint(CenterlineVertices[i], TubeIntersectionData.LateralTube, CurrentStartRing, CurrentEndRing);

			FVector CurrentVertexCenterPoint = CenterlineVertices[i];

			float CurrentVertexRingRadius = InterpolatedRingRadius(CurrentVertexCenterPoint, CurrentStartRing, CurrentEndRing);

			if (i == /*StartIndex*/0) {
				PreviousLateralIntersectionRing = FirstLateralIntersectionRing;
			}
			else {
				PreviousLateralIntersectionRing = CurrentLateralIntersectionRing;
				/*int32 ReorderedIntersectionNum = ReorderedIntersection.Num();
				TempPreviousVertices.Add(ReorderedIntersection[i]);
				TempPreviousVertices.Append(CurrentLateralIntersectionRing.Vertices);*/

				/*if (EndIndexAlt >= ReorderedIntersectionNum) {
					EndIndexAlt = 0;
				}*/

				/*TempPreviousVertices.Add(ReorderedIntersection[EndIndexAlt]);
				PreviousLateralIntersectionRing.Vertices = TempPreviousVertices;*/
			}

			/*EndIndexAlt++;*/

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
					bPointsIntersect = true;
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

			//TubeIntersectionData.LateralTubeIntersectionRings.Rings.Add(PreviousLateralIntersectionRing);
			TubeIntersectionData.LateralTubeIntersectionRings.Rings.Add(CurrentLateralIntersectionRing);
		}

		// < -------------------------- ************************************************************************************************************************************************
		//< ---------------------------------- maybe remove this bool and have this run everytime. migh connect better with the other problem we were having
		// < -------------------------- ************************************************************************************************************************************************
		//if (bPointsIntersect) {
			FRingData FirstFullLateralRing;
			for (const FRingData& Ring : TubeIntersectionData.LateralTubeRemovedVertices.Rings) {
				if (Ring.bIsComplete) {
					FirstFullLateralRing = Ring;
					break;
				}
			}

			TubeIntersectionData.LateralTubeFirstFullRing = FirstFullLateralRing;

			TArray<FVector> TempCurrentVertices;
			TArray<FVector> LastRingVertices = TubeIntersectionData.LateralTubeIntersectionRings.Rings.Last().Vertices;

			int32 NumLastRingVertices = LastRingVertices.Num();

			for (int32 i = 0; i < NumLastRingVertices; i++) {
				FRingData TempLastStartRing;
				FRingData TempLastEndRing;
				FVector TempLastVertex = LastRingVertices[i];
				FindSegmentForPoint(TempLastVertex, TubeIntersectionData.LateralTube, TempLastStartRing, TempLastEndRing);
				FVector TempLastCenterlinePoint = CalculateCenterLinePoint(TempLastVertex, TempLastStartRing, TempLastEndRing);
				FVector TempCurrentDirection = (TempLastVertex - TempLastCenterlinePoint).GetSafeNormal();
				FVector TempCurrentVertex = FindIntersectionOnNewRing(FirstFullLateralRing.Vertices, TempCurrentDirection, FirstFullLateralRing.Center, TempLastVertex);
				TempCurrentVertices.Add(TempCurrentVertex);
			}

			FirstFullLateralRing.Vertices = TempCurrentVertices;
			TubeIntersectionData.LateralTubeIntersectionRings.Rings.Add(FirstFullLateralRing);
		//}
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

	//	// --- 1. Compute the Ring’s Normal ---
	//// Assume the ring is planar. If there are at least 3 vertices, compute a normal.
	//	FVector RingNormal = FVector::UpVector; // default fallback
	//	if (RingVertices.Num() >= 3)
	//	{
	//		RingNormal = FVector::CrossProduct(RingVertices[1] - RingVertices[0], RingVertices[2] - RingVertices[0]).GetSafeNormal();
	//	}

	//	// --- 2. Build a 2D Basis for the Ring’s Plane ---
	//	// Choose X-axis as the direction from Center to the first ring vertex.
	//	FVector XAxis = (RingVertices[0] - Center).GetSafeNormal();
	//	// Y-axis is perpendicular to XAxis within the plane.
	//	FVector YAxis = FVector::CrossProduct(RingNormal, XAxis).GetSafeNormal();

	//	// --- 3. Helper: Project a 3D point to 2D coordinates in the ring’s plane ---
	//	auto ProjectTo2D = [Center, XAxis, YAxis](const FVector& V) -> FVector2D
	//		{
	//			FVector Vec = V - Center;
	//			return FVector2D(FVector::DotProduct(Vec, XAxis), FVector::DotProduct(Vec, YAxis));
	//		};

	//	// --- 4. Project the ring vertices and the ray direction into 2D ---
	//	TArray<FVector2D> ProjectedVertices;
	//	for (const FVector& V : RingVertices)
	//	{
	//		ProjectedVertices.Add(ProjectTo2D(V));
	//	}
	//	// Project the direction. (Direction is a vector so we treat it like a point from the origin.)
	//	FVector2D Dir2D = FVector2D(FVector::DotProduct(Direction, XAxis), FVector::DotProduct(Direction, YAxis));

	//	// --- 5. Set up the ray ---
	//	// Here we interpret the ray as originating from 'Center' and going in the 2D direction Dir2D.
	//	// (This matches your original code that uses Center in the computation.)
	//	FVector2D RayOrigin2D = ProjectTo2D(Center);
	//	FVector2D RayDir2D = Dir2D.GetSafeNormal();

	//	// --- 6. Find intersections between the ray and each segment of the ring ---
	//	// We use standard 2D line intersection math.
	//	// Given:
	//	//   Segment: P + t*(R), with t in [0, 1]
	//	//   Ray: Q + u*(S), with u >= 0
	//	// Then the intersection satisfies:
	//	//   P + t*R = Q + u*S
	//	// and we can solve:
	//	//   t = Cross2D(Q - P, S) / Cross2D(R, S)
	//	//   u = Cross2D(Q - P, R) / Cross2D(R, S)
	//	//
	//	// We'll define a helper lambda for the signed 2D cross product:
	//	auto Cross2D = [](const FVector2D& A, const FVector2D& B) -> float
	//		{
	//			return A.X * B.Y - A.Y * B.X;
	//		};

	//	TArray<FVector2D> IntersectionPoints2D;
	//	int32 VertexCount = ProjectedVertices.Num();
	//	for (int32 i = 0; i < VertexCount; ++i)
	//	{
	//		FVector2D P = ProjectedVertices[i];
	//		FVector2D R = ProjectedVertices[(i + 1) % VertexCount] - P; // segment vector

	//		// Ray defined as Q + u * S, with Q = RayOrigin2D and S = RayDir2D.
	//		FVector2D Q = RayOrigin2D;
	//		FVector2D S = RayDir2D;

	//		// Compute denominator (the cross product of R and S).
	//		float denom = Cross2D(R, S);
	//		if (FMath::IsNearlyZero(denom))
	//		{
	//			continue; // They are parallel; no unique intersection.
	//		}

	//		FVector2D QminusP = Q - P;
	//		float t = Cross2D(QminusP, S) / denom; // parameter along the segment
	//		float u = Cross2D(QminusP, R) / denom; // parameter along the ray

	//		// If t is within the segment and u is along the ray, we have an intersection.
	//		if (t >= 0.0f && t <= 1.0f && u >= 0.0f)
	//		{
	//			FVector2D Intersection2D = P + t * R;
	//			IntersectionPoints2D.Add(Intersection2D);
	//		}
	//	}

	//	// --- 7. Handle No Intersection ---
	//	if (IntersectionPoints2D.Num() == 0)
	//	{
	//		UE_LOG(LogTemp, Log, TEXT("No intersection found on the ring."));
	//		return FVector::ZeroVector;
	//	}

	//	// --- 8. Choose the Intersection Closest to the TargetVertex ---
	//	// Project the target vertex into our 2D space.
	//	FVector2D Target2D = ProjectTo2D(TargetVertex);
	//	FVector2D ClosestPoint2D = IntersectionPoints2D[0];
	//	float MinDistance2D = FVector2D::DistSquared(ClosestPoint2D, Target2D);

	//	for (const FVector2D& P2D : IntersectionPoints2D)
	//	{
	//		float dist2D = FVector2D::DistSquared(P2D, Target2D);
	//		if (dist2D < MinDistance2D)
	//		{
	//			ClosestPoint2D = P2D;
	//			MinDistance2D = dist2D;
	//		}
	//	}

	//	// --- 9. Convert the 2D intersection back to 3D ---
	//	// In our coordinate system, any 2D point (x, y) corresponds to:
	//	//   3D point = Center + x*XAxis + y*YAxis.
	//	FVector ClosestPoint3D = Center + ClosestPoint2D.X * XAxis + ClosestPoint2D.Y * YAxis;
	//	return ClosestPoint3D;
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

	bool LineSegmentIntersectsQuadrilateral(
		const FVector& LineStart,
		const FVector& LineEnd,
		const FVector& V0,
		const FVector& V1,
		const FVector& V2,
		const FVector& V3,
		FVector& OutIntersectionPoint
	) {
		// compute plane normal using cross product of two edges
		FVector Edge1 = V1 - V0;
		FVector Edge2 = V3 - V0;
		FVector PlaneNormal = FVector::CrossProduct(Edge1, Edge2).GetSafeNormal();

		// compute plane equation: N . (P - V0) = 0
		float PlaneD = -FVector::DotProduct(PlaneNormal, V0);

		// compute line direction
		FVector LineDir = LineEnd - LineStart;
		float Denom = FVector::DotProduct(PlaneNormal, LineDir);

		// check if the line is parallel to the plane
		if (FMath::Abs(Denom) < KINDA_SMALL_NUMBER) {
			//UE_LOG(LogTemp, Log, TEXT("Line is parallel to the plane V0: %s, V1: %s, V2: %s, V3: %s, LineStart: %s, LineEnd: %s"), *V0.ToString(), *V1.ToString(), *V2.ToString(), *V3.ToString(), *LineStart.ToString(), *LineEnd.ToString());
			return false;
		}

		// compute the intersection point with the plane
		float T = -(FVector::DotProduct(PlaneNormal, LineStart) + PlaneD) / Denom;
		if (T < 0.0f || T > 1.0f) {
			//UE_LOG(LogTemp, Log, TEXT("Intersection is outside the line segment V0: %s, V1: %s, V2: %s, V3: %s, LineStart: %s, LineEnd: %s"), *V0.ToString(), *V1.ToString(), *V2.ToString(), *V3.ToString(), *LineStart.ToString(), *LineEnd.ToString());
			return false; // intersection outside the line segment
		}

		OutIntersectionPoint = LineStart + T * LineDir;

		//// check if the intersection point is inside the quadrilateral
		//FVector QuadEdges[4] = { V1 - V0, V2 - V1, V3 - V2, V0 - V3 };
		//FVector ToPoint[4] = { 
		//	OutIntersectionPoint - V0, 
		//	OutIntersectionPoint - V1, 
		//	OutIntersectionPoint - V2, 
		//	OutIntersectionPoint - V3 
		//};

		//for (int32 i = 0; i < 4; i++) {
		//	FVector Cross = FVector::CrossProduct(QuadEdges[i], ToPoint[i]);
		//	if (FVector::DotProduct(Cross, PlaneNormal) < 0) {
		//		UE_LOG(LogTemp, Log, TEXT("Point is outside the quadrilateral V0: %s, V1: %s, V2: %s, V3: %s, LineStart: %s, LineEnd: %s"), *V0.ToString(), *V1.ToString(), *V2.ToString(), *V3.ToString(), *LineStart.ToString(), *LineEnd.ToString());
		//		return false; // point is outside the quadrilateral
		//	}
		//}

		//UE_LOG(LogTemp, Log, TEXT("Intersection is inside the quadrilateral V0: %s, V1: %s, V2: %s, V3: %s, LineStart: %s, LineEnd: %s"), *V0.ToString(), *V1.ToString(), *V2.ToString(), *V3.ToString(), *LineStart.ToString(), *LineEnd.ToString());

		//return true; // intersection confirmed inside the quadrilateral

		// Check if the intersection point is inside the quadrilateral using the winding number method
		FVector QuadVertices[4] = { V0, V1, V2, V3 };
		float TotalAngle = 0.0f;

		for (int i = 0; i < 4; i++) {
			FVector A = (QuadVertices[i] - OutIntersectionPoint).GetSafeNormal();
			FVector B = (QuadVertices[(i + 1) % 4] - OutIntersectionPoint).GetSafeNormal();
			float Angle = FMath::Acos(FVector::DotProduct(A, B));
			TotalAngle += Angle;
		}

		// If the total angle is close to 2pi (360 degrees), the point is inside
		if (FMath::Abs(TotalAngle - 2.0f * PI) < 0.01f) {
			//UE_LOG(LogTemp, Warning, TEXT("Intersection is inside the quadrilateral V0: %s, V1: %s, V2: %s, V3: %s, LineStart: %s, LineEnd: %s"), *V0.ToString(), *V1.ToString(), *V2.ToString(), *V3.ToString(), *LineStart.ToString(), *LineEnd.ToString());
			return true;
		}

		//UE_LOG(LogTemp, Log, TEXT("Point is outside the quadrilateral V0: %s, V1: %s, V2: %s, V3: %s, LineStart: %s, LineEnd: %s"), *V0.ToString(), *V1.ToString(), *V2.ToString(), *V3.ToString(), *LineStart.ToString(), *LineEnd.ToString());

		return false;
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

		FRingData IntersectionHalfWayStartRing;
		FRingData IntersectionHalfWayEndRing;

		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.CardinalVertices[1], TubeIntersectionData.MainTube, IntersectionHalfWayStartRing, IntersectionHalfWayEndRing);

		//FVector IntersectionHalfWay = CalculateCenterLinePoint(TubeIntersectionData.IntersectionRing.CardinalVertices[1], IntersectionHalfWayStartRing, IntersectionHalfWayEndRing);

		int32 StartHalfwayIndex = FindRingIndexByCenter(TubeIntersectionData.MainTube.Rings, IntersectionHalfWayStartRing.Center);
		int32 EndHalfwayIndex = FindRingIndexByCenter(TubeIntersectionData.MainTube.Rings, IntersectionHalfWayEndRing.Center);

		PartialIntersectionRings.Add(TubeIntersectionData.MTBelowIntersectionRingPartial);

		int32 MainTubeRingIndex = 0;
		int32 BaseIndex = 0;

		bool AboveBelowAdded = false;

		for (const FRingData Ring : TubeIntersectionData.MainTube.Rings) {

			if (Ring.bIsComplete && MainTubeRingIndex <= /*(TubeIntersectionData.MainTube.Rings.Num() / 2) - 1*/StartHalfwayIndex) {
				BelowIntersectionRings.Add(Ring);
			}
			else if (!Ring.bIsComplete) {
				PartialIntersectionRings.Add(Ring);
			}
			else if (Ring.bIsComplete && MainTubeRingIndex > /*(TubeIntersectionData.MainTube.Rings.Num() / 2) - 1*/EndHalfwayIndex) {
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

		//ConnectPartialRingArrayPaired(ReversedLateralTubeIntersectionRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
		ConnectRingArray(ReversedLateralTubeIntersectionRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);

		FullLateralTubeRings.Add(TubeIntersectionData.LateralTubeFirstFullRing);

		for (const FRingData Ring : TubeIntersectionData.LateralTubeRemovedVertices.Rings) {
			if (Ring.bIsComplete) {
				FullLateralTubeRings.Add(Ring);
			}
		}

		// this including the part in the struct was for testing only, can remove all of this
		TubeIntersectionData.MainTubePartialRings = BelowIntersectionRings;

		ConnectRingArray(FullLateralTubeRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
	}

	void OrderSquareIntersectionConnectionsOneCorner(
		const FTubeData& MainTube,
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

		// < ------------------ rework all of this so it works more like the lateral tube intersection rings. we find a centerline point for each of the points on the corner of the ring
		// < ------------------ put these into an array and order them properly, then use these to construct the rings, this should avoid alot of the problems that made this crazy
		// should make this its own function when we are doing that stuff
		TArray<FVector> CenterlineVertices;
		for (int32 i = 0; i < RingVertices.Num(); i++) {
			FRingData CurrentStartRing;
			FRingData CurrentEndRing;
			FindSegmentForPoint(RingVertices[i], MainTube, CurrentStartRing, CurrentEndRing);

			FVector CurrentVertexCenterPoint = CalculateCenterLinePoint(RingVertices[i], CurrentStartRing, CurrentEndRing);
			CenterlineVertices.Add(CurrentVertexCenterPoint);
		}

		UE_LOG(LogTemp, Log, TEXT("Centerline Vertices Created"));
		for (int32 i = 0; i < CenterlineVertices.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("i: %d, Vector: %s"), i, *CenterlineVertices[i].ToString());
		}

		RemoveDuplicateVertices(CenterlineVertices);

		// log the intersection centerpoints here
		UE_LOG(LogTemp, Log, TEXT("Centerline Vertices Dupes remove"));
		for (int32 i = 0; i < CenterlineVertices.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("i: %d, Vector: %s"), i, *CenterlineVertices[i].ToString());
		}

		if (!Reversed) {
			OrderVerticesByDistanceToPoint(CenterlineVertices, MainTube.EndPosition);
		}
		else {
			OrderVerticesByDistanceToPoint(CenterlineVertices, MainTube.StartPosition);
		}

		UE_LOG(LogTemp, Log, TEXT("Centerline Vertices Ordered"));
		for (int32 i = 0; i < CenterlineVertices.Num(); ++i) {
			UE_LOG(LogTemp, Log, TEXT("i: %d, Vector: %s"), i, *CenterlineVertices[i].ToString());
		}

		TArray<TArray<FVector>> TempVertexArrays;
		TArray<FVector> PreviousVertexArray = RingVertices;
		TempVertexArrays.Add(PreviousVertexArray);

		for (int32 i = 0; i < CenterlineVertices.Num(); i++) {
			FRingData CurrentStartRing;
			FRingData CurrentEndRing;
			FVector TempNextVertex;

			TArray<FVector> CurrentVertexArray;

			FindSegmentForPoint(CenterlineVertices[i], MainTube, CurrentStartRing, CurrentEndRing);

			float TempCurrentRingRadius = InterpolatedRingRadius(CenterlineVertices[i], CurrentStartRing, CurrentEndRing);

			// *** NEW: Compute a local up vector for the generated ring ***
			// Instead of using the global MainTube.UpVector, derive a local up from the actual geometry.
			// We choose a "representative" previous vertex from the array. In our case, we simply take the first vertex
			// from the PreviousVertexArray. (You may choose another strategy if needed.)
			//FVector RepresentativePrevVertex;
			//if (PreviousVertexArray.Num() > 0) {
			//	RepresentativePrevVertex = PreviousVertexArray[0];  // Use first valid vertex if available
			//}
			//else {
			//	RepresentativePrevVertex = CenterlineVertices[i] + (MainTube.UpVector * TempCurrentRingRadius);
			//	// Fallback: Place the representative vertex slightly above the centerline vertex
			//}

			// The local up vector points from the centerline to the representative previous vertex.
			// This ensures that the generated ring is oriented to match the actual intersection.
			FVector LocalUp = CurrentStartRing.UpVector;

			FRingData TempCurrentRing;
			GenerateRing(
				CenterlineVertices[i],
				MainTube.Direction,
				MainTube.UpVector,
				TempCurrentRingRadius,
				MainTube.NumSides,
				TempCurrentRing
			);

			for (int32 j = 0; j < PreviousVertexArray.Num(); j++) {
				FVector TempPreviousVertex = PreviousVertexArray[j];

				FRingData TempCurrentStartRing;
				FRingData TempCurrentEndRing;
				FindSegmentForPoint(TempPreviousVertex, MainTubeIntersectionRings, TempCurrentStartRing, TempCurrentEndRing);
				FVector TempCurrentCenterlinePoint = CalculateCenterLinePoint(TempPreviousVertex, TempCurrentStartRing, TempCurrentEndRing);
				FVector TempCurrentDirection = (TempPreviousVertex - TempCurrentCenterlinePoint).GetSafeNormal();

				FVector TempCurrentVertex = FindIntersectionOnNewRing(TempCurrentRing.Vertices, TempCurrentDirection, TempCurrentRing.Center, TempPreviousVertex);

				FRingData CurrentLateralTubeStartRing;
				FRingData CurrentLateralTubeEndRing;

				FindSegmentForPoint(TempCurrentVertex, LateralTube, CurrentLateralTubeStartRing, CurrentLateralTubeEndRing);

				int32 TempNextIndex = (j + 1) % PreviousVertexArray.Num();
				int32 TempPrevIndex = (j - 1 + PreviousVertexArray.Num()) % PreviousVertexArray.Num();

				int32 CurrentLateralTubeStartRingIndex = FindRingIndexByCenter(LateralTube.Rings, CurrentLateralTubeStartRing.Center);

				int32 LateralTubeNumRings = LateralTube.Rings.Num();
				FRingData PreviousLateralTubeRing = LateralTube.Rings[(CurrentLateralTubeStartRingIndex - 1 + LateralTubeNumRings) % LateralTubeNumRings];
				FRingData NextLateralTubeRing = LateralTube.Rings[(CurrentLateralTubeStartRingIndex + 2) % LateralTubeNumRings];

				FVector StartOrEnd;

				if (Reversed) {
					StartOrEnd = MainTube.EndPosition;
				}
				else {
					StartOrEnd = MainTube.StartPosition;
				}

				float DistanceToEndCurrentCenter = FVector::DistSquared(StartOrEnd, TempCurrentRing.Center);
				float DistanceToEndPrevVertexCenter = FVector::DistSquared(StartOrEnd, TempCurrentCenterlinePoint);

				bool IsInsideFrustum = IsPointInsideFrustum(CurrentLateralTubeStartRing, CurrentLateralTubeEndRing, TempCurrentVertex) ||
					IsPointInsideFrustum(PreviousLateralTubeRing, CurrentLateralTubeStartRing, TempCurrentVertex) ||
					IsPointInsideFrustum(CurrentLateralTubeEndRing, NextLateralTubeRing, TempCurrentVertex);
				bool IsAlongLine = IsPointOnLine(TempPreviousVertex, PreviousVertexArray[TempNextIndex], TempCurrentVertex) ||
					IsPointOnLine(PreviousVertexArray[TempPrevIndex], TempPreviousVertex, TempCurrentVertex);
				bool IsInsideFrustumNotOnLine = IsInsideFrustum/* && !IsAlongLine*/;
				bool IsCurrentCloserToEnd = DistanceToEndCurrentCenter < DistanceToEndPrevVertexCenter;

				if (IsInsideFrustumNotOnLine || !IsCurrentCloserToEnd) {
					UE_LOG(LogTemp, Log, TEXT("Intersection Corner Generation inside tube j: %d, TempCurrentVertex: %s, PreviousVertex: %s, CenterlinePoint: %s, StartRingCenter: %s, EndRingCenter: %s"), j, *TempCurrentVertex.ToString(), *PreviousVertexArray[j].ToString(), *TempCurrentCenterlinePoint.ToString(), *TempCurrentStartRing.Center.ToString(), *TempCurrentEndRing.Center.ToString());
					FVector ClosestExteriorVertex = TempPreviousVertex;
					//int32 ClosestExteriorIndex;

					// < --------------------------------need to check if the point is below the current ring using the ring indexing, and find ring index by center.
					// i think we need to check each vertex as well. problem is that it correctly picks new ones on line, then one is inside because of weird shape and it puts it back on the intersection ring so it jumps up
					if (IsInsideFrustumNotOnLine && IsCurrentCloserToEnd) {
						//UE_LOG(LogTemp, Log, TEXT("--------------- > it hit here %d times"), testCount);
						//testCount++;
						bool bFoundOne = false;
						for (const FVector& Vertex : PreviousVertexArray) {
							if (IsPointOnRing(TempCurrentRing, Vertex)) {
								UE_LOG(LogTemp, Log, TEXT("Point is bellow current ring: %d"), i);
								CurrentVertexArray.Add(Vertex);
								bFoundOne = true;
								break;
							}
						}

						if (!bFoundOne && j > 0) {
							FVector PreviousVertex = CurrentVertexArray[j - 1];
							CurrentVertexArray.Add(PreviousVertex);
						}
						/*else {
							CurrentVertexArray.Add(ClosestExteriorVertex);
						}*/
					}
					else {
						CurrentVertexArray.Add(ClosestExteriorVertex);
					}

					//CurrentVertexArray.Add(ClosestExteriorVertex);

					//FindClosestVertex(TempCurrentVertex, PreviousVertexArray, ClosestExteriorVertex, ClosestExteriorIndex);
					//CurrentVertexArray.Add(ClosestExteriorVertex);
					//UE_LOG(LogTemp, Log, TEXT("Point is inside tube, using vertex from previous array at i: %d"), ClosestExteriorIndex);
				}
				else {
					UE_LOG(LogTemp, Log, TEXT("Intersection Corner Generation outside tube j: %d, TempCurrentVertex: %s, PreviousVertex: %s, CenterlinePoint: %s, StartRingCenter: %s, EndRingCenter: %s"), j, *TempCurrentVertex.ToString(), *PreviousVertexArray[j].ToString(), *TempCurrentCenterlinePoint.ToString(), *TempCurrentStartRing.Center.ToString(), *TempCurrentEndRing.Center.ToString());
					CurrentVertexArray.Add(TempCurrentVertex);
				}
			}

			TempVertexArrays.Add(CurrentVertexArray);
			PreviousVertexArray = CurrentVertexArray;
		}

		//int32 CurrentIndex = StartIndex;
		////int32 testCount = 0;

		//TArray<TArray<FVector>> TempVertexArrays;
		//TArray<FVector> PreviousVertexArray = RingVertices;
		//TempVertexArrays.Add(PreviousVertexArray);

		///*UE_LOG(LogTemp, Log, TEXT("PreviousVertexArray contains %d vertices:"), PreviousVertexArray.Num());
		//for (int32 i = 0; i < PreviousVertexArray.Num(); ++i)
		//{
		//	UE_LOG(LogTemp, Log, TEXT("[%d] %s"), i, *PreviousVertexArray[i].ToString());
		//}*/

		//while (true) {

		//	FRingData TempNextStartRing;
		//	FRingData TempNextEndRing;
		//	FVector TempNextVertex;

		//	TArray<FVector> CurrentVertexArray;

		//	for (int32 i = 0; i < PreviousVertexArray.Num(); ++i) {
		//		FRingData TempCurrentStartRing;
		//		FRingData TempCurrentEndRing;
		//		FindSegmentForPoint(PreviousVertexArray[i], MainTubeIntersectionRings, TempCurrentStartRing, TempCurrentEndRing);
		//		FVector TempCurrentCenterlinePoint = CalculateCenterLinePoint(PreviousVertexArray[i], TempCurrentStartRing, TempCurrentEndRing);
		//		FVector TempCurrentDirection = (PreviousVertexArray[i] - TempCurrentCenterlinePoint).GetSafeNormal();

		//		//int32 CurrentEndRingIndex = FindRingIndexByCenter(MainTubeIntersectionRings.Rings, TempCurrentEndRing.Center);

		//		FVector TempCurrentVertex = FindIntersectionOnNewRing(MainTubeIntersectionRings.Rings[CurrentIndex].Vertices, TempCurrentDirection, MainTubeIntersectionRings.Rings[CurrentIndex].Center, PreviousVertexArray[i]);

		//		FRingData CurrentLateralTubeStartRing;
		//		FRingData CurrentLateralTubeEndRing;

		//		FindSegmentForPoint(TempCurrentVertex, LateralTube, CurrentLateralTubeStartRing, CurrentLateralTubeEndRing);

		//		int32 CurrentLateralTubeStartRingIndex = FindRingIndexByCenter(LateralTube.Rings, CurrentLateralTubeStartRing.Center);

		//		/*FRingData PreviousLateralTubeRing = LateralTube.Rings[CurrentLateralTubeStartRingIndex - 1];
		//		FRingData NextLateralTubeRing = LateralTube.Rings[CurrentLateralTubeStartRingIndex + 2];*/

		//		int32 TempNextIndex = (i + 1) % PreviousVertexArray.Num();
		//		int32 TempPrevIndex = (i - 1 + PreviousVertexArray.Num()) % PreviousVertexArray.Num();

		//		FVector StartOrEnd;

		//		if (Reversed) {
		//			StartOrEnd = MainTubeIntersectionRings.EndPosition;
		//		}
		//		else {
		//			StartOrEnd = MainTubeIntersectionRings.StartPosition;
		//		}

		//		float DistanceToEndCurrentCenter = FVector::DistSquared(StartOrEnd, MainTubeIntersectionRings.Rings[CurrentIndex].Center);
		//		float DistanceToEndPrevVertexCenter = FVector::DistSquared(StartOrEnd, TempCurrentCenterlinePoint);

		//		bool IsInsideFrustum = IsPointInsideFrustum(CurrentLateralTubeStartRing, CurrentLateralTubeEndRing, TempCurrentVertex) /*||
		//			IsPointInsideFrustum(PreviousLateralTubeRing, CurrentLateralTubeStartRing, TempCurrentVertex) ||
		//			IsPointInsideFrustum(CurrentLateralTubeEndRing, NextLateralTubeRing, TempCurrentVertex)*/;
		//		bool IsAlongLine = IsPointOnLine(PreviousVertexArray[i], PreviousVertexArray[TempNextIndex], TempCurrentVertex) ||
		//			IsPointOnLine(PreviousVertexArray[TempPrevIndex], PreviousVertexArray[i], TempCurrentVertex);
		//		bool IsInsideFrustumNotOnLine = IsInsideFrustum && !IsAlongLine;
		//		bool IsCurrentCloserToEnd = DistanceToEndCurrentCenter < DistanceToEndPrevVertexCenter;

		//		if (IsInsideFrustumNotOnLine || !IsCurrentCloserToEnd) {
		//			UE_LOG(LogTemp, Log, TEXT("Intersection Corner Generation inside tube i: %d, TempCurrentVertex: %s, PreviousVertex: %s, CenterlinePoint: %s, StartRingCenter: %s, EndRingCenter: %s"), i, *TempCurrentVertex.ToString(), *PreviousVertexArray[i].ToString(), *TempCurrentCenterlinePoint.ToString(), *TempCurrentStartRing.Center.ToString(), *TempCurrentEndRing.Center.ToString());
		//			FVector ClosestExteriorVertex = PreviousVertexArray[i];
		//			//int32 ClosestExteriorIndex;

		//			// < --------------------------------need to check if the point is below the current ring using the ring indexing, and find ring index by center.
		//			// i think we need to check each vertex as well. problem is that it correctly picks new ones on line, then one is inside because of weird shape and it puts it back on the intersection ring so it jumps up
		//			if (IsInsideFrustumNotOnLine && IsCurrentCloserToEnd) {
		//				//UE_LOG(LogTemp, Log, TEXT("--------------- > it hit here %d times"), testCount);
		//				//testCount++;
		//				for (const FVector& Vertex : PreviousVertexArray) {
		//					if (IsPointOnRing(MainTubeIntersectionRings.Rings[CurrentIndex], Vertex)) {
		//						//UE_LOG(LogTemp, Log, TEXT("Point is bellow current ring: %d"), CurrentIndex);
		//						CurrentVertexArray.Add(Vertex);
		//						break;
		//					}
		//				}
		//			}
		//			else {
		//				CurrentVertexArray.Add(ClosestExteriorVertex);
		//			}

		//			//FindClosestVertex(TempCurrentVertex, PreviousVertexArray, ClosestExteriorVertex, ClosestExteriorIndex);
		//			//CurrentVertexArray.Add(ClosestExteriorVertex);
		//			//UE_LOG(LogTemp, Log, TEXT("Point is inside tube, using vertex from previous array at i: %d"), ClosestExteriorIndex);
		//		}
		//		else {
		//			UE_LOG(LogTemp, Log, TEXT("Intersection Corner Generation outside tube i: %d, TempCurrentVertex: %s, PreviousVertex: %s, CenterlinePoint: %s, StartRingCenter: %s, EndRingCenter: %s"), i, *TempCurrentVertex.ToString(), *PreviousVertexArray[i].ToString(), *TempCurrentCenterlinePoint.ToString(), *TempCurrentStartRing.Center.ToString(), *TempCurrentEndRing.Center.ToString());
		//			CurrentVertexArray.Add(TempCurrentVertex);
		//		}

		//		
		//	}

		//	TempVertexArrays.Add(CurrentVertexArray);
		//	PreviousVertexArray = CurrentVertexArray;

		//	// Check the stopping condition
		//	if (Reversed)
		//	{
		//		if (CurrentIndex >= MainTubeIntersectionRings.Rings.Num() - 1) {
		//			break;
		//		}
		//		++CurrentIndex; // Move forward
		//	}
		//	else
		//	{
		//		if (CurrentIndex <= 0) {
		//			break;
		//		}
		//		--CurrentIndex; // Move backward
		//	}
		//}

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

		FVector TempTopRightPoint = FindIntersectionOnNewRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TopRightDirection, TubeIntersectionData.MTAboveIntersectionRing.Center, ReorderedIntersection[RightIndex]);
		FVector TempBottomRightPoint = FindIntersectionOnNewRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, BottomRightDirection, TubeIntersectionData.MTBelowIntersectionRing.Center, ReorderedIntersection[BottomIndex]);
		FVector TempBottomLeftPoint = FindIntersectionOnNewRing(TubeIntersectionData.MTBelowIntersectionRing.Vertices, BottomLeftDirection, TubeIntersectionData.MTBelowIntersectionRing.Center, ReorderedIntersection[LeftIndex]);
		FVector TempTopLeftPoint = FindIntersectionOnNewRing(TubeIntersectionData.MTAboveIntersectionRing.Vertices, TopLeftDirection, TubeIntersectionData.MTAboveIntersectionRing.Center, ReorderedIntersection[0]);

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

		// < ---------------------------------------------------------- need to add a connection vertex to both the first and last vertex of partial rings
		// We need to remove the square around the intersection as it is now. Need to go back to the last corner instead. 
		// get perpendicular directions for top and bottom rings
		TArray<FVector> TopDirections, BottomDirections;
		TopDirections.Add((TubeIntersectionData.IntersectionSquare.Corners[3] - TubeIntersectionData.MTAboveIntersectionRing.Center).GetSafeNormal());
		TopDirections.Add((TubeIntersectionData.IntersectionSquare.Corners[2] - TubeIntersectionData.MTAboveIntersectionRing.Center).GetSafeNormal());

		BottomDirections.Add((TubeIntersectionData.IntersectionSquare.Corners[0] - TubeIntersectionData.MTBelowIntersectionRing.Center).GetSafeNormal());
		BottomDirections.Add((TubeIntersectionData.IntersectionSquare.Corners[1] - TubeIntersectionData.MTBelowIntersectionRing.Center).GetSafeNormal());

		FVector RightCardVertCenterPoint = CalculateCenterLinePoint(TubeIntersectionData.IntersectionRing.CardinalVertices[1], StartCenterRing, EndCenterRing);
		float RightCardVertRingRadius = InterpolatedRingRadius(RightCardVertCenterPoint, StartCenterRing, EndCenterRing);

		FRingData RightCardVertRing;

		GenerateRing(
			RightCardVertCenterPoint,
			StartCenterRing.Direction,
			StartCenterRing.UpVector,
			RightCardVertRingRadius,
			TubeIntersectionData.MainTube.NumSides,
			RightCardVertRing
		);

		TArray<FVector> ModifiedRingVerts;

		int32 indexCount = 0;
		int32 indexCountFinal = 0;

		bool Corner1Added = false;
		//bool Corner01Added = false;
		for (const FVector& Vertex : RightCardVertRing.Vertices) {
			FVector VertexDirection = (Vertex - RightCardVertRing.Center).GetSafeNormal();

			// check if the vertex direction is between the two interpolated directions
			FVector Normal = FVector::CrossProduct(BottomDirections[0], BottomDirections[1]).GetSafeNormal(); // Normal of the plane
			bool IsBetween = FVector::DotProduct(Normal, FVector::CrossProduct(BottomDirections[0], VertexDirection)) >= 0 &&
				FVector::DotProduct(Normal, FVector::CrossProduct(VertexDirection, BottomDirections[1])) >= 0;

			if (!IsBetween) {
				if (!Corner1Added) {
					Corner1Added = true;
				}
				ModifiedRingVerts.Add(Vertex);
				indexCount++;
			}
			else if (Corner1Added) {
				Corner1Added = false;
				indexCountFinal = indexCount;
			}
		}

		TArray<FVector> ReorderedModifiedRingVerts = DoobContainerUtils::ReorderedArray(ModifiedRingVerts, indexCountFinal);

		TempTopRightRingVertices.Add(ReorderedModifiedRingVerts[0]);
		FVector RightDirection = (ReorderedModifiedRingVerts[0] - RightCardVertRing.Center).GetSafeNormal();
		FVector TopRightLastVertex = FindIntersectionOnNewRing(
			TubeIntersectionData.MTAboveIntersectionRing.Vertices, 
			RightDirection, 
			TubeIntersectionData.MTAboveIntersectionRing.Center, 
			ReorderedModifiedRingVerts[0]
		);
		TempTopRightSquareVertices.Add(TopRightLastVertex);

		//TubeIntersectionData.TestVerts = ReorderedModifiedRingVerts;

		TempBottomRightRingVertices.Insert(ReorderedModifiedRingVerts[0], 0);
		FVector BottomRightFirstVertex = FindIntersectionOnNewRing(
			TubeIntersectionData.MTBelowIntersectionRing.Vertices,
			RightDirection,
			TubeIntersectionData.MTBelowIntersectionRing.Center,
			ReorderedModifiedRingVerts[0]
		);
		TempBottomRightSquareVertices.Insert(BottomRightFirstVertex, 0);

		//int32 ReorderedModifiedRingLastIndex = ReorderedModifiedRingVerts.Num() - 1;

		TempBottomLeftRingVertices.Add(ReorderedModifiedRingVerts.Last());
		FVector LeftDirection = (ReorderedModifiedRingVerts.Last() - RightCardVertRing.Center).GetSafeNormal();
		FVector BottomLeftLastVertex = FindIntersectionOnNewRing(
			TubeIntersectionData.MTBelowIntersectionRing.Vertices,
			LeftDirection,
			TubeIntersectionData.MTBelowIntersectionRing.Center,
			ReorderedModifiedRingVerts.Last()
		);
		TempBottomLeftSquareVertices.Add(BottomLeftLastVertex);

		TempTopLeftRingVertices.Insert(ReorderedModifiedRingVerts.Last(), 0);
		FVector TopLeftLastVertex = FindIntersectionOnNewRing(
			TubeIntersectionData.MTAboveIntersectionRing.Vertices,
			LeftDirection,
			TubeIntersectionData.MTAboveIntersectionRing.Center,
			ReorderedModifiedRingVerts.Last()
		);
		TempTopLeftSquareVertices.Insert(ReorderedModifiedRingVerts.Last(), 0);

		OrderSquareIntersectionConnectionsOneCorner(
			TubeIntersectionData.MainTube,
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempTopRightRingVertices,
			TempTopRightSquareVertices,
			EndIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.TopRightPartialRings,
			true
		);

		UE_LOG(LogTemp, Log, TEXT("Top Right Intersection connections"));

		OrderSquareIntersectionConnectionsOneCorner(
			TubeIntersectionData.MainTube,
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempBottomRightRingVertices,
			TempBottomRightSquareVertices,
			StartIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.BottomRightPartialRings,
			false
		);

		UE_LOG(LogTemp, Log, TEXT("Bottom Right Intersection connections"));

		// < ----------------------- this reversal and the one after conrner is connected are janky as hell, need to fix
		// NOTE: only this corner is having this issue, index 0 for the second and all after iterations are inaccurate as all get up cannot fathom why

		Algo::Reverse(TempBottomLeftRingVertices);
		Algo::Reverse(TempBottomLeftSquareVertices);

		OrderSquareIntersectionConnectionsOneCorner(
			TubeIntersectionData.MainTube,
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempBottomLeftRingVertices,
			TempBottomLeftSquareVertices,
			StartIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.BottomLeftPartialRings,
			false
		);

		TArray<TArray<FVector>> TempBottomLeftPartialRings = TubeIntersectionData.IntersectionSquare.BottomLeftPartialRings;

		for (int32 i = 0; i < TempBottomLeftPartialRings.Num(); i++) {
			Algo::Reverse(TempBottomLeftPartialRings[i]);
		}

		TubeIntersectionData.IntersectionSquare.BottomLeftPartialRings = TempBottomLeftPartialRings;

		// < ----------------------------------------------------------------------------------------------------------------

		UE_LOG(LogTemp, Log, TEXT("Bottom Left Intersection connections"));

		OrderSquareIntersectionConnectionsOneCorner(
			TubeIntersectionData.MainTube,
			TempMainTubeIntersectionRings,
			TubeIntersectionData.LateralTube,
			TempTopLeftRingVertices,
			TempTopLeftSquareVertices,
			EndIntersectionRingsStartIndex,
			TubeIntersectionData.IntersectionSquare.TopLeftPartialRings,
			true
		);

		UE_LOG(LogTemp, Log, TEXT("Top Left Intersection connections"));
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

	void MergeClosePoints(TArray<FVector>& Points, float MinDistance) {
		bool bMerged = false;

		for (int32 i = 0; i < Points.Num(); ++i) {
			for (int32 j = i + 1; j < Points.Num(); ++j) {
				if (FVector::Dist(Points[i], Points[j]) < MinDistance) {
					FVector MidPoint = (Points[i] + Points[j]) * 0.5f;

					// remove original points
					Points.RemoveAt(j); // removes last first to avoid issues
					Points.RemoveAt(i);

					// add midpoint
					Points.Add(MidPoint);

					bMerged = true;

					// restart to check no new pairs exist
					MergeClosePoints(Points, MinDistance);
					return;
				}
			}
		}
	}

	void OrderVerticesByDistanceToPoint(TArray<FVector>& Vertices, const FVector& ReferencePoint) {
		// Sort the vertices based on their distance to the ReferencePoint
		Vertices.Sort([&ReferencePoint](const FVector& A, const FVector& B) {
			float DistanceA = FVector::DistSquared(A, ReferencePoint); // Squared distance for performance
			float DistanceB = FVector::DistSquared(B, ReferencePoint);
			return DistanceA < DistanceB; // Sort in ascending order of distance
		});
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

		const float Epsilon = 0.01f; // Adjust epsilon as needed for your application.

		if (DistanceToAxis > InterpolatedRadius + Epsilon) {
			return false; // Preliminary check fails
		}

		// Final polygonal cross-section check
		const TArray<FVector>& CurrentRingVertices = (t < 0.5f) ? StartRing.Vertices : EndRing.Vertices;
		return IsPointInsidePolygonWinding(CurrentRingVertices, Point, ProjectedPoint);
		
		// 5. Blend the two cross-sectional polygons.
		//TArray<FVector> BlendedPolygon;
		//if (StartRing.Vertices.Num() == EndRing.Vertices.Num() && StartRing.Vertices.Num() > 0)
		//{
		//	int32 NumVertices = StartRing.Vertices.Num();
		//	for (int32 i = 0; i < NumVertices; i++)
		//	{
		//		// Interpolate corresponding vertices based on t.
		//		FVector BlendedVertex = FMath::Lerp(StartRing.Vertices[i], EndRing.Vertices[i], t);
		//		BlendedPolygon.Add(BlendedVertex);
		//	}
		//}
		//else
		//{
		//	// Fallback: if the vertex counts differ, choose one based on t.
		//	BlendedPolygon = (t < 0.5f ? StartRing.Vertices : EndRing.Vertices);
		//}

		//// 6. Compute the centroid of the blended polygon.
		//FVector BlendedCentroid = FVector::ZeroVector;
		//for (const FVector& V : BlendedPolygon)
		//{
		//	BlendedCentroid += V;
		//}
		//BlendedCentroid /= BlendedPolygon.Num();

		//// 7. Use the winding number method on the blended polygon.
		//return IsPointInsidePolygonWinding(BlendedPolygon, Point, BlendedCentroid);

		// 4. Interpolate (blend) between the two cross-sectional polygons.
		//TArray<FVector> BlendedPolygon;
		//if (StartRing.Vertices.Num() == EndRing.Vertices.Num() && StartRing.Vertices.Num() > 0)
		//{
		//	int32 NumVertices = StartRing.Vertices.Num();
		//	for (int32 i = 0; i < NumVertices; i++)
		//	{
		//		// Blend corresponding vertices based on t.
		//		FVector BlendedVertex = FMath::Lerp(StartRing.Vertices[i], EndRing.Vertices[i], t);
		//		BlendedPolygon.Add(BlendedVertex);
		//	}
		//}
		//else
		//{
		//	// Fallback: if vertex counts differ, choose one of the polygons based on t.
		//	BlendedPolygon = (t < 0.5f ? StartRing.Vertices : EndRing.Vertices);
		//}

		//// 5. Compute the centroid of the blended polygon.
		//FVector BlendedCentroid = FVector::ZeroVector;
		//for (const FVector& Vertex : BlendedPolygon)
		//{
		//	BlendedCentroid += Vertex;
		//}
		//BlendedCentroid /= BlendedPolygon.Num();

		//// 6. Use the blended polygon and its centroid for the 2D point-in-polygon test.
		//return IsPointInsidePolygon(BlendedPolygon, Point, BlendedCentroid);

		// Choose the appropriate polygon based on t.
		//const TArray<FVector>& CurrentRingVertices = (t < 0.5f) ? StartRing.Vertices : EndRing.Vertices;

		//// Compute the centroid of the chosen polygon.
		//FVector PolygonCentroid = FVector::ZeroVector;
		//for (const FVector& Vertex : CurrentRingVertices)
		//{
		//	PolygonCentroid += Vertex;
		//}
		//PolygonCentroid /= CurrentRingVertices.Num();

		//// Use the polygon's centroid for the final 2D point-in-polygon check.
		//return IsPointInsidePolygon(CurrentRingVertices, Point, PolygonCentroid);



		// --- Begin Interpolated Cross-Sectional Check ---
		// Instead of choosing one polygon, blend the two based on t.
		//TArray<FVector> BlendedPolygon;
		//if (StartRing.Vertices.Num() == EndRing.Vertices.Num() && StartRing.Vertices.Num() > 0)
		//{
		//	int32 NumPoly = StartRing.Vertices.Num();
		//	for (int32 i = 0; i < NumPoly; i++)
		//	{
		//		// Interpolate each corresponding vertex.
		//		FVector BlendedVertex = FMath::Lerp(StartRing.Vertices[i], EndRing.Vertices[i], t);
		//		BlendedPolygon.Add(BlendedVertex);
		//	}
		//}
		//else
		//{
		//	// Fallback: if the two polygons have different counts,
		//	// pick one based on t.
		//	BlendedPolygon = (t < 0.5f ? StartRing.Vertices : EndRing.Vertices);
		//}

		//// Compute the centroid of the blended polygon.
		//FVector BlendedCentroid = FVector::ZeroVector;
		//for (const FVector& V : BlendedPolygon)
		//{
		//	BlendedCentroid += V;
		//}
		//BlendedCentroid /= BlendedPolygon.Num();

		//// Use the blended polygon and its centroid for the final polygonal check.
		//return IsPointInsidePolygon(BlendedPolygon, Point, BlendedCentroid);
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

	// Winding number–based point-in-polygon test.
	bool IsPointInsidePolygonWinding(const TArray<FVector>& Polygon, const FVector& Point, const FVector& Centroid) {
		int32 NumVertices = Polygon.Num();
		if (NumVertices < 3)
		{
			return false;
		}

		// Build a local 2D coordinate system for the polygon.
		// Use the first two edges to compute a normal and then define U, V axes.
		FVector Edge1 = Polygon[1] - Polygon[0];
		FVector Edge2 = Polygon[2] - Polygon[0];
		FVector Normal = FVector::CrossProduct(Edge1, Edge2).GetSafeNormal();
		FVector U = Edge1.GetSafeNormal();
		FVector V = FVector::CrossProduct(Normal, U).GetSafeNormal();

		// Project the polygon vertices into 2D.
		TArray<FVector2D> Projected;
		for (const FVector& Vertex : Polygon)
		{
			FVector Local = Vertex - Centroid;
			Projected.Add(FVector2D(FVector::DotProduct(Local, U), FVector::DotProduct(Local, V)));
		}

		// Project the point into 2D.
		FVector LocalPoint = Point - Centroid;
		FVector2D PPoint(FVector::DotProduct(LocalPoint, U), FVector::DotProduct(LocalPoint, V));

		// Compute the winding number.
		int WindingNumber = 0;
		for (int i = 0; i < NumVertices; i++)
		{
			int j = (i + 1) % NumVertices;
			// Check if the edge crosses the horizontal line through PPoint.
			if (Projected[i].Y <= PPoint.Y)
			{
				if (Projected[j].Y > PPoint.Y)
				{
					// Compute the cross product (in 2D) to determine if PPoint is left of the edge.
					float Cross = (Projected[j].X - Projected[i].X) * (PPoint.Y - Projected[i].Y) -
						(Projected[j].Y - Projected[i].Y) * (PPoint.X - Projected[i].X);
					if (Cross > 0)
					{
						WindingNumber++;
					}
				}
			}
			else
			{
				if (Projected[j].Y <= PPoint.Y)
				{
					float Cross = (Projected[j].X - Projected[i].X) * (PPoint.Y - Projected[i].Y) -
						(Projected[j].Y - Projected[i].Y) * (PPoint.X - Projected[i].X);
					if (Cross < 0)
					{
						WindingNumber--;
					}
				}
			}
		}

		return WindingNumber != 0;
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