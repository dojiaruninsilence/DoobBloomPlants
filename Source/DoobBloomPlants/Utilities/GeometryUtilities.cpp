#include "GeometryUtilities.h"

namespace GeometryUtilities {
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
}