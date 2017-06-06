#include <algorithm>
#include <chrono>
#include <climits>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include "../characteristicMatrix.h"
#include "../miscFunctions.h"
#include "../MurmurHash2.cpp"

#define SIMILARITY_THRESHOLD 0.8
#define NUM_BINS 8

void
computeSignatureMatrix(std::vector<int> &signatureMatrix, ccsMatrix* characteristicMatrix, int numSets, int numShingles, int numBins, long int currentTime)
{
  int binSize = numShingles/numBins;  
  if (numShingles % numBins) binSize++;

  for (int i = 0; i < numSets; i++) {
    int offSetCM = characteristicMatrix -> col_ptr[i];
    for (int j = offSetCM; j < characteristicMatrix -> col_ptr[i+1]; j++) {
      int shingleIdx = characteristicMatrix -> row_ind[j];
      int shingleNewIdx = MurmurHash2(&shingleIdx, 4, currentTime)%numShingles;      
      int binIdx = shingleNewIdx/binSize;
      int offSetSM = binIdx + (i*numBins);
      signatureMatrix[offSetSM] = std::min(signatureMatrix[offSetSM], shingleNewIdx);
    }
  }
}

int
nestedLoopJoin(std::vector<int> signatureMatrix, int rSize, int sSize, std::vector<std::string> relationRSetsID, std::vector<std::string> relationSSetsID, int numBins, std::vector<int>& result)
{
  int similarPairs = 0, numSets = rSize+sSize;

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
//	  std::cout << relationRSetsID[i] << " " << relationSSetsID[j-rSize] << " " << similarity << "\n";
	  //Results
	  result.push_back(i);
	  result.push_back(j-rSize);
	  similarPairs++;
      }
    }
  }
  return similarPairs;
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
  long int currentTime = std::chrono::system_clock::now().time_since_epoch().count();

  int numSets = rSize + sSize;
  int numShingles = shingles.size();
  int numBins = NUM_BINS;  

  Timer tPreProc; tPreProc.start();
  //  std::cout << "Size of full matrix: " << numShingles*numSets << std::endl;
  ccsMatrix* characteristicMatrix = buildCharacteristicMatrix(shingleSetMap);
  //  std::cout << "Size of charMatrix: " << ccsMatrixSize(characteristicMatrix) << std::endl;
  //  std::cout << "Compact/full: " << 100-((ccsMatrixSize(characteristicMatrix)*1.0 / (numShingles*numSets)*1.0)*100) << "%" << std::endl;
  tPreProc.stop();
  
  if (numSets != (int)num_unique_keys(shingleSetMap) || numSets+1 != (int)(characteristicMatrix -> col_ptr.size())) {
    std::cout << "Error in dataset. Possible problems: record containing only stop words or empty record." << std::endl;
  }
  
  std::vector<int> signatureMatrix (numSets*numBins, INT_MAX);
  Timer tMin; tMin.start();
  computeSignatureMatrix(signatureMatrix, characteristicMatrix, numSets, numShingles, numBins, currentTime);
  tMin.stop();

  // Print signature matrix
  // for (int i = 0; i < numSets; i++) {
  //     for (int j = 0; j < numBins; j++) {
  // 	  std::cout << signatureMatrix[i*numBins + j] << " ";
  //     }
  //     std::cout << std::endl;
  // }

  Timer tJoin; tJoin.start();
  int similarPairsCount = 0;
  std::vector<int> result;
  similarPairsCount = nestedLoopJoin(signatureMatrix, rSize, sSize, relationRSetsID, relationSSetsID, numBins, result);

  // Print results
  // std::cout << "Similar pairs: " << result.size() / 2 << std::endl;
  // for (uint i = 0; i < result.size()-1 ; i+=2) {
  //   std::cout << relationRSetsID[result[i]] << " " << relationSSetsID[result[i+1]] << std::endl;
  // }
  tJoin.stop();
  
  delete(characteristicMatrix);

  t.stop();

  std::cout << t.elapsed_time() << "," << tReadFromDisk.elapsed_time() << "," << tPreProc.elapsed_time() << "," << tMin.elapsed_time() << "," << tJoin.elapsed_time() << "," << 0 << "," << similarPairsCount << std::endl;

  return 0;
}
