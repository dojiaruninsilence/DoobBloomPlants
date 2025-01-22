// Fill out your copyright notice in the Description page of Project Settings.


#include "DoobContainerUtils.h"

namespace DoobContainerUtils {
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
}
