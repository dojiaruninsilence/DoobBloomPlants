#include "DoobGeometryUtils.h"

#include "Math/Vector.h"
#include "Math/Transform.h"
#include "Containers/Array.h"
#include "DrawDebugHelpers.h"

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

		RemoveDuplicateVertices(CombinedVertices);

		OutRingData.CombinedVertices = OrderRingVertices(CombinedVertices);

		OutRingData.Centroid = ComputeCentroid(OutRingData.CombinedVertices);

		/*FVector FirstVertex = OutRingData.CombinedVertices[0];
		OutRingData.CombinedVertices.Add(FirstVertex);*/
	}

	void GenerateHalfIntersectionRing(
		const FTubeData& TubeA,
		const FTubeData& TubeB,
		TArray<FVector>& OutRingVertices,
		float Precision
	) {
		//// Precompute Tube B bounding volume for quick checks
		//FBox TubeBBox(ForceInit);
		//for (const FRingData& Ring : TubeB.Rings) {
		//	TubeBBox += Ring.Center; 
		//}


		for (int32 LineIndexA = 0; LineIndexA < TubeA.Rings.Num() - 1; ++LineIndexA) {
			const FRingData& CurrentRingA = TubeA.Rings[LineIndexA];
			const FRingData& NextRingA = TubeA.Rings[LineIndexA + 1];

			// Check if either ring center is inside Tube B's bounding box
			//if (!TubeBBox.IsInside(CurrentRingA.Center) && !TubeBBox.IsInside(NextRingA.Center)) {
			//	continue; // Skip if both rings are outside Tube B's bounding box
			//}

			int32 NumVerticesA = CurrentRingA.Vertices.Num();

			for (int32 VertexIndexA = 0; VertexIndexA < NumVerticesA; ++VertexIndexA) {
				FVector LineStartA = CurrentRingA.Vertices[VertexIndexA];
				FVector LineEndA = NextRingA.Vertices[VertexIndexA];

				UE_LOG(LogTemp, Log, TEXT("Generate Intersection Ring Line Start: %s, Line End: %s"), *LineStartA.ToString(), *LineEndA.ToString());

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

						UE_LOG(LogTemp, Log, TEXT("V0: %s, V1: %s, V2: %s, V3: %s"), *V0.ToString(), *V1.ToString(), *V3.ToString(), *V3.ToString());

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

		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.CardinalVertices[2], TubeIntersectionData.MainTube, LowVertStartRing, LowVertEndRing);
		FindSegmentForPoint(TubeIntersectionData.IntersectionRing.CardinalVertices[0], TubeIntersectionData.MainTube, HighVertStartRing, HighVertEndRing);

		FVector LowestVertexCenter = CalculateCenterLinePoint(TubeIntersectionData.IntersectionRing.CardinalVertices[2], LowVertStartRing, LowVertEndRing);
		FVector HighestVertexCenter = CalculateCenterLinePoint(TubeIntersectionData.IntersectionRing.CardinalVertices[0], HighVertStartRing, HighVertEndRing);

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
		TArray<FVector> ReorderedIntersection = ReorderedArray(
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
			ReorderedIntersection = ReorderedArray(
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
			FVector VertexToAdd = ReorderedIntersection[i];\

			FirstLateralIntersectionRing.Vertices.Add(VertexToAdd);

			FVector SuperTempCenter = CalculateCenterLinePoint(VertexToAdd, IntersectionRingStartRing, IntersectionRingEndRing);
			FVector SuperTempDirection = (VertexToAdd - SuperTempCenter).GetSafeNormal();
		}

		int32 EndIndex = 0;

		FRingData CurrentLateralIntersectionRing;
		FRingData PreviousLateralIntersectionRing;

		for (int32 i = StartIndex; i > EndIndex; --i) {
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

				/*FVector TempCurrentVertex = TempNextVertex;*/

				FRingData CurrentMainTubeStartRing;
				FRingData CurrentMainTubeEndRing;

				FindSegmentForPoint(TempCurrentVertex, TubeIntersectionData.MainTube, CurrentMainTubeStartRing, CurrentMainTubeEndRing);

				int32 TempNextIndex = (j + 1) % PreviousLateralIntersectionRing.Vertices.Num();
				int32 TempPrevIndex = (j - 1) % PreviousLateralIntersectionRing.Vertices.Num();

				bool IsInsideFrustum = IsPointInsideFrustum(CurrentMainTubeStartRing, CurrentMainTubeEndRing, TempCurrentVertex);
				bool IsAlongLine = IsPointOnLine(PreviousLateralIntersectionRing.Vertices[j], PreviousLateralIntersectionRing.Vertices[TempNextIndex], TempCurrentVertex) ||
					IsPointOnLine(PreviousLateralIntersectionRing.Vertices[TempPrevIndex], PreviousLateralIntersectionRing.Vertices[j], TempCurrentVertex);

				if (IsInsideFrustum && !IsAlongLine) {
					TempCurrentVertices.Add(PreviousLateralIntersectionRing.Vertices[j]);
					UE_LOG(LogTemp, Log, TEXT("Vertex in main tube: %s"), *PreviousLateralIntersectionRing.Vertices[j].ToString());
				}
				else {
					TempCurrentVertices.Add(TempCurrentVertex);
					UE_LOG(LogTemp, Log, TEXT("Vertex NOT in main tube: %s"), *PreviousLateralIntersectionRing.Vertices[j].ToString());
				}

				/*TempCurrentVertices.Add(TempCurrentVertex);*/

				
				FindSegmentForPoint(PreviousLateralIntersectionRing.Vertices[TempNextIndex], TubeIntersectionData.LateralTube, TempNextStartRing, TempNextEndRing);
				FVector TempNextCenterlinePoint = CalculateCenterLinePoint(PreviousLateralIntersectionRing.Vertices[TempNextIndex], TempNextStartRing, TempNextEndRing);
				FVector TempNextDirection = (PreviousLateralIntersectionRing.Vertices[TempNextIndex] - TempNextCenterlinePoint).GetSafeNormal();

				TempNextVertex = FindIntersectionOnNewRing(CurrentLateralIntersectionRing.Vertices, TempNextDirection, CurrentLateralIntersectionRing.Center, PreviousLateralIntersectionRing.Vertices[TempNextIndex]);

				/*FVector TempLineIntersectionVertex;
				if (DoLinesIntersect(TempCurrentVertex, TempNextVertex, PreviousLateralIntersectionRing.Vertices[j], PreviousLateralIntersectionRing.Vertices[TempNextIndex], TempLineIntersectionVertex)) {
					TempCurrentVertices.Add(TempLineIntersectionVertex);
					UE_LOG(LogTemp, Log, TEXT("Lines intersect at vertex: %s"), *TempLineIntersectionVertex.ToString());
				}*/
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

		// Debug log: Output the generated intersection point
		UE_LOG(LogTemp, Log, TEXT("Generated Intersection: %s"), *OutIntersectionPoint.ToString());

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

	void ConnectIntersectionRingToSquare(FTwoTubeIntersectionData& TubeIntersectionData, int32& BaseIndex) {
		// combine arrays
		TArray<FVector> TempRing = TubeIntersectionData.IntersectionSquare.TopRightVertices;
		TempRing.Append(TubeIntersectionData.IntersectionSquare.BottomRightVertices);
		TempRing.Append(TubeIntersectionData.IntersectionSquare.BottomLeftVertices);
		TempRing.Append(TubeIntersectionData.IntersectionSquare.TopLeftVertices);
		TArray<FVector> TempSquare = TubeIntersectionData.IntersectionRing.TopRightVertices;
		TempSquare.Append(TubeIntersectionData.IntersectionRing.BottomRightVertices);
		TempSquare.Append(TubeIntersectionData.IntersectionRing.BottomLeftVertices);
		TempSquare.Append(TubeIntersectionData.IntersectionRing.TopLeftVertices);

		if (TempRing.Num() != TempSquare.Num()) {
			UE_LOG(LogTemp, Warning, TEXT("Vertex count mismatch: TempRing (%d), TempSquare (%d)"), TempRing.Num(), TempSquare.Num());
			//return;
		}

		for (int32 i = 0; i < TempRing.Num() - 1; ++i) {
			int32 CurrentA = BaseIndex + i;
			int32 NextA = BaseIndex + (i + 1);

			int32 CurrentB = BaseIndex + TempRing.Num() + i;
			int32 NextB = BaseIndex + TempRing.Num() + (i + 1);

			// Triangle 1
			TubeIntersectionData.Triangles.Add(CurrentA);
			TubeIntersectionData.Triangles.Add(CurrentB);
			TubeIntersectionData.Triangles.Add(NextA);

			// Triangle 2
			TubeIntersectionData.Triangles.Add(NextA);
			TubeIntersectionData.Triangles.Add(CurrentB);
			TubeIntersectionData.Triangles.Add(NextB);
		}

		TubeIntersectionData.AllVertices.Append(TempRing);
		TubeIntersectionData.AllVertices.Append(TempSquare);
		BaseIndex += TempRing.Num() + TempSquare.Num();
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

		//RemoveInternalVertices(TubeIntersectionData.MainTube, TubeIntersectionData.LateralTube);

		FullLateralTubeRings.Add(TubeIntersectionData.LateralTubeFirstFullRing);

		for (const FRingData Ring : TubeIntersectionData.LateralTubeRemovedVertices.Rings) {
			if (Ring.bIsComplete) {
				FullLateralTubeRings.Add(Ring);
			}
		}

		TubeIntersectionData.MainTubePartialRings = FullLateralTubeRings;

		ConnectRingArray(FullLateralTubeRings, TubeIntersectionData.AllVertices, TubeIntersectionData.Triangles, BaseIndex);
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

	bool IsPointInsideFrustum(const FRingData& StartRing, const FRingData& EndRing, const FVector& Point) {
		// compute the height vector and its squared length
		FVector LengthVector = EndRing.Center - StartRing.Center;
		float LengthSquared = LengthVector.SizeSquared();

		//UE_LOG(LogTemp, Warning, TEXT("Length Vector: %s, Length Squared: %f"), *LengthVector.ToString(), LengthSquared);

		// project the point onto the height vector
		FVector PointToStart = Point - StartRing.Center;
		float DotProduct = FVector::DotProduct(PointToStart, LengthVector);
		float t = DotProduct / LengthSquared;

		//UE_LOG(LogTemp, Warning, TEXT("Point: %s, Projection Factor (t): %f"), *Point.ToString(), t);

		// check if the projection is within the frustum's height
		if (t < 0.0f || t > 1.0f)
		{
			//UE_LOG(LogTemp, Warning, TEXT("Point is outside the frustum's height."));
			return false;
		}

		// interpolate the radius at the projected height
		float InterpolatedRadius = FMath::Lerp(StartRing.Radius, EndRing.Radius, t);
		//UE_LOG(LogTemp, Warning, TEXT("Interpolated Radius: %f"), InterpolatedRadius);

		// compute the projected point on the frustum axis
		FVector ProjectedPoint = StartRing.Center + t * LengthVector;
		//UE_LOG(LogTemp, Warning, TEXT("Projected Point: %s"), *ProjectedPoint.ToString());

		// Check the distance from the axis
		float DistanceToAxis = FVector::Dist(Point, ProjectedPoint);
		//UE_LOG(LogTemp, Warning, TEXT("Distance to Axis: %f, Interpolated Radius: %f"), DistanceToAxis, InterpolatedRadius);

		// Return whether the point is within the radius
		bool bIsInside = DistanceToAxis <= InterpolatedRadius;
		/*UE_LOG(LogTemp, Warning, TEXT("Point %s is %s the frustum."),
			*Point.ToString(),
			bIsInside ? TEXT("inside") : TEXT("outside"));*/

		return bIsInside;
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
				else {
					UE_LOG(LogTemp, Log, TEXT("-----------------------Point Removed for inside frustum------------------------"));
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

		/*for (int32 TubeRingIndexA = 0; TubeRingIndexA < TubeA.Rings.Num() - 1; ++TubeRingIndexA) {
			FRingData CurrentRingA = TubeA.Rings[TubeRingIndexA];
			FRingData NextRingA = TubeA.Rings[TubeRingIndexA + 1];

			for (int32 TubeRingIndexB = 0; TubeRingIndexB < TubeB.Rings.Num(); ++TubeRingIndexB) {
				FRingData CurrentRingB = TubeB.Rings[TubeRingIndexB];
				FRingData TempRing;

				for (int32 VertexIndex = 0; VertexIndex < CurrentRingB.Vertices.Num(); ++VertexIndex) {
					FVector CurrentVertex = CurrentRingB.Vertices[VertexIndex];

					if (!IsPointInsideFrustum(CurrentRingA, NextRingA, CurrentVertex)) {
						TempRing.Vertices.Add(CurrentVertex);
					}
					else {
						UE_LOG(LogTemp, Log, TEXT("-----------------------Point Removed for inside frustum------------------------"));
					}
				}

				if (TempRing.Vertices.Num() > 0) {
					TempRing.Center = CurrentRingB.Center;
					TempRing.Direction = CurrentRingB.Direction;
					TempRing.Normal = CurrentRingB.Normal;
					TempRing.Radius = CurrentRingB.Radius;
					TempRing.RingID = CurrentRingB.RingID;
					TempRing.UpVector = CurrentRingB.UpVector;

					if (TempRing.Vertices.Num() == CurrentRingB.Vertices.Num()) {
						TempRing.bIsClosed = CurrentRingB.bIsClosed;
						TempRing.bIsComplete = CurrentRingB.bIsComplete;
					}
					else if (TempRing.Vertices.Num() < CurrentRingB.Vertices.Num()) {
						TempRing.bIsClosed = false;
						TempRing.bIsComplete = false;
					}

					TempTube.Add(TempRing);
				}				
			}
		}*/

		TubeB.Rings = TempTube;

		// Log Updated TubeB.Rings
		UE_LOG(LogTemp, Warning, TEXT("Updated TubeB.Rings:"));
		for (const auto& Ring : TubeB.Rings) {
			for (const auto& Vertex : Ring.Vertices) {
				UE_LOG(LogTemp, Warning, TEXT("Vertex: X=%f, Y=%f, Z=%f"), Vertex.X, Vertex.Y, Vertex.Z);
			}
		}
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

		//int32 ArraySize = IntersectionRing.CombinedVertices.Num();

		//if (ArraySize == 0) {
		//	UE_LOG(LogTemp, Error, TEXT("CombinedVertices array is empty!"));
		//	return;
		//}

		//int32 LowestIndex;
		//FVector LowestVertex;
		//DoobGeometryUtils::FindClosestVertex(StartCenter, IntersectionRing.CombinedVertices, LowestVertex, LowestIndex);
		//int32 HighestIndex;
		//FVector HighestVertex;
		//DoobGeometryUtils::FindClosestVertex(EndCenter, IntersectionRing.CombinedVertices, HighestVertex, HighestIndex);

		///*int32 RightIndex = (HighestIndex + LowestIndex) / 2;
		//int32 LeftIndex = LowestIndex + RightIndex;*/

		////// Wrap around indices
		////int32 RightIndex = (HighestIndex + LowestIndex / 2) % ArraySize; // Middle index, wrapping if necessary
		////int32 LeftIndex = (LowestIndex + RightIndex) % ArraySize;

		////// Safeguard against invalid indices
		////if (RightIndex >= ArraySize || LeftIndex >= ArraySize || RightIndex < 0 || LeftIndex < 0) {
		////	UE_LOG(LogTemp, Error, TEXT("Calculated indices are out of bounds: RightIndex=%d, LeftIndex=%d"), RightIndex, LeftIndex);
		////	return;
		////}

		///*int32 LeftIndexAdjust = 0;
		//if (ArraySize % 4 == 0) {
		//	LeftIndexAdjust = 1;
		//}*/

		//TArray<FVector> ReorderedLowRing = ReorderedArray(IntersectionRing.CombinedVertices, LowestIndex);
		//TArray<FVector> ReorderedHighRing = ReorderedArray(IntersectionRing.CombinedVertices, HighestIndex);

		//int32 TempHighIndex = 0;
		//int32 TempLowIndex = FindVertexIndex(ReorderedHighRing, LowestVertex);

		//int32 TempRightIndex = (TempHighIndex + TempLowIndex) / 2 + 1;

		////TempLowIndex = FindVertexIndex(ReorderedLowRing, HighestVertex);

		//int32 TempLeftIndex = (ArraySize - TempRightIndex) /*+ LeftIndexAdjust*/;

		//FVector TempRightVertex = ReorderedHighRing[TempRightIndex];
		//FVector TempLeftVertex = ReorderedHighRing[TempLeftIndex];

		//int32 RightIndex = FindVertexIndex(IntersectionRing.CombinedVertices, TempRightVertex);
		//int32 LeftIndex = FindVertexIndex(IntersectionRing.CombinedVertices, TempLeftVertex);

		///*FVector TempRightVertex = ReorderedHighRing[ArraySize / 4];
		//FVector TempLeftVertex = ReorderedHighRing[(ArraySize / 2) + (ArraySize / 4)];

		//int32 RightIndex = FindVertexIndex(IntersectionRing.CombinedVertices, TempRightVertex);
		//int32 LeftIndex = FindVertexIndex(IntersectionRing.CombinedVertices, TempLeftVertex);*/

		//FVector RightVertex = TempRightVertex;
		//FVector LeftVertex = TempLeftVertex;

		//IntersectionRing.CardinalVertices = { HighestVertex, RightVertex, LowestVertex, LeftVertex };
		//IntersectionRing.CardinalIndices = { HighestIndex, RightIndex, LowestIndex, LeftIndex };
	}

	int32 FindVertexIndex(const TArray<FVector>& Vertices, const FVector& TargetVertex) {
		for (int32 i = 0; i < Vertices.Num(); i++) {
			if (Vertices[i] == TargetVertex) {
				return i;
			}
		}

		return -1;
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

		TArray<FVector> ReorderedIntersection = ReorderedArray(
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

		/*Algo::Reverse(TempBottomLeftSquareVertices);
		Algo::Reverse(TempBottomLeftRingVertices);*/

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

		TubeIntersectionData.IntersectionRing.TopRightVertices = TempTopRightRingVertices;
		TubeIntersectionData.IntersectionRing.BottomRightVertices = TempBottomRightRingVertices;
		TubeIntersectionData.IntersectionRing.BottomLeftVertices = TempBottomLeftRingVertices;
		TubeIntersectionData.IntersectionRing.TopLeftVertices = TempTopLeftRingVertices;

		TubeIntersectionData.IntersectionSquare.TopRightVertices = TempTopRightSquareVertices;
		TubeIntersectionData.IntersectionSquare.BottomRightVertices = TempBottomRightSquareVertices;
		TubeIntersectionData.IntersectionSquare.BottomLeftVertices = TempBottomLeftSquareVertices;
		TubeIntersectionData.IntersectionSquare.TopLeftVertices = TempTopLeftSquareVertices;
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
				indexCountFinal = indexCount;
			}
		}
		//TempTopRing.Vertices.Add(IntersectionSquare.Corners[3]);

		TArray<FVector> ReorderedTopVertices = ReorderedArray(TempTopRing.Vertices, indexCountFinal);
		TempTopRing.Vertices = ReorderedTopVertices;

		FRingData TempBottomRing = TubeIntersectionData.MTBelowIntersectionRing;
		TempBottomRing.Vertices.Empty();

		indexCount = 0;
		indexCountFinal = 0;

		// Flags to track corner insertion
		bool Corner0Added = false;
		bool Corner1Added = false;
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
				indexCountFinal = indexCount;
			}
		}

		TArray<FVector> ReorderedBottomVertices = ReorderedArray(TempBottomRing.Vertices, indexCountFinal);
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

			// skip rings outside of intersection rings
			/*if (RingProjection < BottomProjection || RingProjection > TopProjection) {
				TempTube.Add(Ring);
				continue;
			}*/

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
					indexCountFinal = indexCount;
				}
			}

			TArray<FVector> ReorderedModifiedVertices = ReorderedArray(ModifiedRing.Vertices, indexCountFinal);
			ModifiedRing.Vertices = ReorderedModifiedVertices;
			ModifiedRing.bIsClosed = false;
			ModifiedRing.bIsComplete = false;

			TempTube.Add(ModifiedRing);
		}

		TubeIntersectionData.MainTube.Rings = TempTube;
	}

	// ---------------------------------------------------------- need to add this to a container utils file ------------------------------------------------------------------//
	TArray<FVector> ReorderedArray(const TArray<FVector>& OriginalArray, int32 StartIndex) {
		// Make sure the start index is within bounds
		if (StartIndex < 0 || StartIndex >= OriginalArray.Num())
		{
			return OriginalArray;  // Return the original array if the index is invalid
		}

		// Create the reordered array
		TArray<FVector> Reordered;

		// Add the second part (from start index to the end)
		for (int32 i = StartIndex; i < OriginalArray.Num(); ++i)
		{
			Reordered.Add(OriginalArray[i]);
		}

		// Add the first part (from the beginning to the start index)
		for (int32 i = 0; i < StartIndex; ++i)
		{
			Reordered.Add(OriginalArray[i]);
		}

		return Reordered;
	}

	// -------------------------------------------------------- need to build a test and check this part -----------------------------------------------------------------------//
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
		/*TArray<FVector> UniqueVertices;

		for (const FVector& Vertex : Vertices) {
			bool bIsDuplicate = false;

			for (const FVector& UniqueVertex : UniqueVertices) {
				if (FVector::DistSquared(Vertex, UniqueVertex) < FMath::Square(Tolerance)) {
					bIsDuplicate = true;
					break;
				}
			}

			if (!bIsDuplicate) {
				UniqueVertices.Add(Vertex);
			}
		}

		Vertices = UniqueVertices;*/
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

	/*bool DoLinesIntersect(
		const FVector& Line1Start,
		const FVector& Line1End,
		const FVector& Line2Start,
		const FVector& Line2End,
		FVector& OutIntersectionPoint,
		float Tolerance
	) {
		// Calculate direction vectors for the lines
		FVector Dir1 = Line1End - Line1Start;
		FVector Dir2 = Line2End - Line2Start;

		// Check if the lines are parallel
		FVector CrossDir = FVector::CrossProduct(Dir1, Dir2);
		if (CrossDir.SizeSquared() < Tolerance) {
			return false; // Lines are parallel, no intersection
		}

		// Solve for intersection parameters (s, t) using a system of equations
		FVector LineDiff = Line2Start - Line1Start;
		float A = FVector::DotProduct(LineDiff, FVector::CrossProduct(Dir2, CrossDir));
		float B = FVector::DotProduct(Dir1, CrossDir);

		if (FMath::Abs(B) < Tolerance) {
			return false; // Lines are coplanar but do not intersect
		}

		float T = A / B;
		if (T < 0.0f || T > 1.0f) {
			return false; // Intersection is outside the bounds of Line1
		}

		FVector Intersection = Line1Start + T * Dir1;

		// Check if Intersection is within the bounds of Line2
		FVector DiffToIntersection = Intersection - Line2Start;
		float ParamOnLine2 = FVector::DotProduct(DiffToIntersection, Dir2) / Dir2.SizeSquared();

		if (ParamOnLine2 < 0.0f || ParamOnLine2 > 1.0f) {
			return false; // Intersection is outside the bounds of Line2
		}

		// Valid intersection
		OutIntersectionPoint = Intersection;
		return true;
	}*/

	// -------------------------------------------------------- need to build a test and check this part end -----------------------------------------------------------------------//

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