#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <math.h>
#include <omp.h>
#include <random>
#include <string>
#include <vector>
#include "../characteristicMatrix.h"
#include "../miscFunctions.h"
#include "../MurmurHash2.cpp"

#define SIMILARITY_THRESHOLD 0.8
#define NUM_BINS 256

void
computeSignatureMatrix(std::vector<int> &signatureMatrix, ccsMatrix* characteristicMatrix, int numSets, int numShingles, int numBins)
{
  int binSize = numShingles/numBins;
  if (numShingles % numBins) binSize++;

  #pragma omp parallel for
  for (uint i = 0; i < signatureMatrix.size(); i++) {
    signatureMatrix[i] = INT_MAX;
  }

  #pragma omp parallel for
  for (int i = 0; i < numSets; i++) {
    int offSetCM = characteristicMatrix -> col_ptr[i];
    for (int j = offSetCM; j < characteristicMatrix -> col_ptr[i+1]; j++) {
      int shingleIdx = characteristicMatrix -> row_ind[j];
      int shingleNewIdx = MurmurHash2(&shingleIdx, 4, 0)%numShingles;
      int binIdx = shingleNewIdx/binSize;
      int offSetSM = binIdx + (i*numBins);
      signatureMatrix[offSetSM] = std::min(signatureMatrix[offSetSM], shingleNewIdx);
    }
  }
}

uint
nestedLoopJoin(std::vector<int> signatureMatrix, int rSize, int sSize, std::vector<std::string> relationRSetsID, std::vector<std::string> relationSSetsID, int numBins, std::vector<int>& result)
{
  int numSets = rSize+sSize;
  std::vector<int> separateResultSizes;
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    std::vector<int> localResult;
#pragma omp single
    {
      separateResultSizes.resize(omp_get_num_threads()+1);
    }
#pragma omp for schedule(dynamic)
    for (int i = 0; i < rSize; i++) {
      for (int j = rSize; j < numSets; j++) {
	int identicalMinhashes = 0;
	int emptyBins = 0;
	for (int k = 0; k < numBins; k++) {
	  if (signatureMatrix[k+(i*numBins)] == signatureMatrix[k+(j*numBins)]) {
	    if (signatureMatrix[k+(i*numBins)] == INT_MAX) {
	      emptyBins++;
	    } else {
	      identicalMinhashes++;
	    }
	  }
	}
	double similarity = (identicalMinhashes*1.0)/((numBins*1.0) - (emptyBins*1.0));
	if (similarity >= SIMILARITY_THRESHOLD) {
	  localResult.push_back(i);
	  localResult.push_back(j-rSize);
	}
      }
    }
    separateResultSizes[tid] = localResult.size();
#pragma omp barrier
#pragma omp single
    {
      int s = 0;
      int numThreads = omp_get_num_threads();
      for (auto i: range(numThreads + 1)) {
	int t = separateResultSizes[i];
	separateResultSizes[i] = s;
	s += t;
      }
      result.resize(separateResultSizes[numThreads]);
    }
    int offset = separateResultSizes[tid];
    std::copy(localResult.begin(), localResult.end(), result.begin() + offset);
  }
  return result.size();
}

int
main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " pathToRelation/relationRFile.data pathToRelation/relationSFile.data\n";
    return 1;
  }
 
  Timer t; t.start();

  std::multimap<int,int> shingleSetMap;
  std::map<std::string,int> shingles;
  std::vector<std::string> relationRSetsID;
  std::vector<std::string> relationSSetsID;

  Timer tReadFromDisk; tReadFromDisk.start();

  int rSize = processInputRelation(shingleSetMap, shingles, relationRSetsID, argv[1], 0);
  int sSize = processInputRelation(shingleSetMap, shingles, relationSSetsID, argv[2], rSize);

  tReadFromDisk.stop();

  int numSets = rSize + sSize;
  int numShingles = shingles.size();
  int numBins = NUM_BINS;  

  Timer tPreProc; tPreProc.start();

  ccsMatrix* characteristicMatrix = buildCharacteristicMatrix(shingleSetMap);

  tPreProc.stop();

  if (numSets != (int)num_unique_keys(shingleSetMap) || numSets+1 != (int)(characteristicMatrix -> col_ptr.size())) {
    std::cout << "Error in dataset. Possible problems: record containing only stop words or empty record." << std::endl;
  }
  
  std::vector<int> signatureMatrix (numSets*numBins);

  Timer tMinhash; tMinhash.start();

  computeSignatureMatrix(signatureMatrix, characteristicMatrix, numSets, numShingles, numBins);

  tMinhash.stop();

  Timer tJoin; tJoin.start();

  std::vector<int> result;
  uint similarPairsCount = nestedLoopJoin(signatureMatrix, rSize, sSize, relationRSetsID, relationSSetsID, numBins, result) / 2;
  // Print results
  //  std::cout << "Similar pairs: " << result.size() / 2 << std::endl;
  // for (int i = 0; i < result.size()-1 ; i+=2) {
  //   std::cout << relationRSetsID[result[i]] << " " << relationSSetsID[result[i+1]] << std::endl;
  // }  

  tJoin.stop();
  
  delete(characteristicMatrix);

  t.stop();

  std::cout << t.elapsed_time() << "," << tReadFromDisk.elapsed_time() << "," << tPreProc.elapsed_time() << "," << tMinhash.elapsed_time() << "," << tJoin.elapsed_time() << "," << 0 << "," << similarPairsCount << std::endl;

  return 0;
}
