#include "DoobMeshUtils.h"

namespace DoobMeshUtils {
	void SmoothMeshVertices(TArray<FVector>& Vertices, const TArray<int32>& Triangles, int32 Iterations) {
		TArray<FVector> SmoothedVertices = Vertices;

		for (int32 Iteration = 0; Iteration < Iterations; ++Iteration) {
			TArray<FVector> TempVertices = SmoothedVertices;

			for (int32 VertexIndex = 0; VertexIndex < Vertices.Num(); ++VertexIndex) {
				FVector AveragePosition = FVector::ZeroVector;
				int32 ConnectedCount = 0;

				// Find neighbors for the current vertex
				for (int32 i = 0; i < Triangles.Num(); i += 3) {
					if (Triangles[i] == VertexIndex || Triangles[i + 1] == VertexIndex || Triangles[i + 2] == VertexIndex) {
						// accumulate neighboring vertices
						if (Triangles[i] != VertexIndex) { AveragePosition += TempVertices[Triangles[i]]; ConnectedCount++; }
						if (Triangles[i + 1] != VertexIndex) { AveragePosition += TempVertices[Triangles[i + 1]]; ConnectedCount++; }
						if (Triangles[i + 2] != VertexIndex) { AveragePosition += TempVertices[Triangles[i + 2]]; ConnectedCount++; }
					}
				}

				if (ConnectedCount > 0) {
					AveragePosition /= ConnectedCount;
					SmoothedVertices[VertexIndex] = FMath::Lerp(TempVertices[VertexIndex], AveragePosition, 0.5f); // Blend towards neighbors
				}
			}
		}

		Vertices = SmoothedVertices;
	}

	void RecalculateNormals(const TArray<FVector>& Vertices, const TArray<int32>& Triangles, TArray<FVector>& Normals) {
		// Ensure there are at least 3 vertices (for a triangle)
		if (Triangles.Num() < 3) return;

		// Initialize normals array
		Normals.SetNumZeroed(Vertices.Num());

		// Loop through all triangles and calculate normals
		for (int32 i = 0; i < Triangles.Num(); i += 3)
		{
			int32 Index0 = Triangles[i];
			int32 Index1 = Triangles[i + 1];
			int32 Index2 = Triangles[i + 2];

			// Get the vertices of the triangle
			FVector Vertex0 = Vertices[Index0];
			FVector Vertex1 = Vertices[Index1];
			FVector Vertex2 = Vertices[Index2];

			// Calculate two edges of the triangle
			FVector Edge1 = Vertex1 - Vertex0;
			FVector Edge2 = Vertex2 - Vertex0;

			// Calculate the normal using the cross product
			FVector TriangleNormal = FVector::CrossProduct(Edge1, Edge2).GetSafeNormal();

			// Add the normal to each of the triangle's vertices
			Normals[Index0] += TriangleNormal;
			Normals[Index1] += TriangleNormal;
			Normals[Index2] += TriangleNormal;
		}

		// Normalize the normals for each vertex
		for (int32 i = 0; i < Normals.Num(); ++i)
		{
			Normals[i] = -Normals[i].GetSafeNormal();
		}
	}

	void RemoveDegenerateTriangles(TArray<FVector>& Vertices, TArray<int32>& Triangles) {
		// tolerance for considering a triangle as degenerate
		const float DegenerateThreshold = 0.001f;

		for (int32 i = Triangles.Num() - 3; i >= 0; i -= 3) {
			int32 Index0 = Triangles[i];
			int32 Index1 = Triangles[i + 1];
			int32 Index2 = Triangles[i + 2];

			// Get the vertices of the triangle
			FVector Vertex0 = Vertices[Index0];
			FVector Vertex1 = Vertices[Index1];
			FVector Vertex2 = Vertices[Index2];

			// Check for degenerate triangle
			FVector Edge1 = Vertex1 - Vertex0;
			FVector Edge2 = Vertex2 - Vertex0;

			if (Edge1.IsNearlyZero() || Edge2.IsNearlyZero() || Edge1.Equals(Edge2, DegenerateThreshold))
			{
				// Remove the degenerate triangle
				Triangles.RemoveAt(i, 3);
			}
		}
	}

	void GenerateTubeMesh(
		const TArray<FVector>& TubeVertices,
		int32 ProfileSegments,
		TArray<int32>& OutTriangles
	) {
		int32 NumProfiles = TubeVertices.Num() / ProfileSegments;

		for (int32 i = 0; i < NumProfiles - 1; ++i) {
			for (int32 j = 0; j < ProfileSegments; ++j) {
				int32 Current = i * ProfileSegments + j;
				int32 Next = Current + ProfileSegments;

				//int32 NextInProfile = (j + 1) % ProfileSegments;

				// Vertex in the same profile but next in sequence
				int32 NextInProfile = i * ProfileSegments + (j + 1) % ProfileSegments;
				// Vertex below the next one in the next profile
				int32 NextInNextProfile = Next + (j + 1) % ProfileSegments;

				// Create two triangles forming a quad
				OutTriangles.Add(Current);
				OutTriangles.Add(Next);
				OutTriangles.Add(NextInProfile);

				OutTriangles.Add(NextInProfile);
				OutTriangles.Add(Next);
				OutTriangles.Add(NextInNextProfile);
			}
		}
	}
}