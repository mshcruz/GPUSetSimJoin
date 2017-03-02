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

// Receives signature matrix as input and group the relation based on the sum of the least significant bits of each part of the signature
std::unordered_map<int, std::vector<int> >
lengthFiltering(std::vector<int> &signatureMatrix, int numSets, int numBins)
{
    std::unordered_map<int, std::vector<int> > candidatePairs;
    // Group sets of R
    for (int i = 0; i < numSets; i++) {
	int sumBits = 0;
	// Sum the LSB of each part of the signature
	for (int j = 0; j < numBins; j++) {
	    sumBits += (signatureMatrix[i*numBins + j] % 2);
	}
	// Check if the sum value is was found before
	if (candidatePairs.find(sumBits) != candidatePairs.end()) {
	    // Add the set ID to the the of sets with that sum
	    candidatePairs[sumBits].push_back(i);
	} else {
	    // Add the first set ID with that sum
	    std::vector<int> sameGroupIDs;	    
	    sameGroupIDs.push_back(i);
	    candidatePairs.insert({sumBits, sameGroupIDs});
	}
    }

    // for (auto bucket: candidatePairs) {
    // 	std::vector<int> sameGroupIDs = bucket.second;
    // 	std::cout << bucket.first << " contains " << sameGroupIDs.size() << std::endl;
    // 	for (uint i = 0; i < sameGroupIDs.size(); i++) {
    // 	    std::cout << bucket.first << " : " <<  sameGroupIDs[i] << std::endl;	    
    // 	}
    // }
    
    
    return candidatePairs;
}
int areFromDifferentRelations(int id1, int id2, int rSize, int sSize) {
    int id1Flag = 0;
    int id2Flag = 0;
    
    if (id1 < rSize) {
	id1Flag = 0; // id1 belongs to R
    } else {
	id1Flag = 1; // id1 belongs to S
    }

    if (id2 < rSize) {
	id2Flag = 0; // id2 belongs to R
    } else {
	id2Flag = 1; // id2 belongs to S
    }
    
    if (id1Flag != id2Flag) {
	return 1;
    } else {
	return 0;
    }
}

int
lengthFilteringJoin(std::vector<int> signatureMatrix, std::unordered_map<int, std::vector<int> > candidatePairs, int rSize, int sSize, std::vector<std::string> relationRSetsID, std::vector<std::string> relationSSetsID, int numBins, std::vector<int>& result)
{
    int similarPairs = 0;
// go through each bucket and calculate the nlj of the pairs inside the bucket
    // uint count = 0;
    for (auto bucket : candidatePairs) {
	// std::cout << "bucket " <<  count << std::endl;
	// count++;
	std::vector<int> sameGroupIDs = bucket.second;
	for (uint i = 0; i < sameGroupIDs.size() && sameGroupIDs[i] < rSize; i++) {
	    for (uint j = 0; j < sameGroupIDs.size(); j++) {
		if (areFromDifferentRelations(sameGroupIDs[i], sameGroupIDs[j], rSize, sSize)) {
		    // std::cout << "Comparing: " << sameGroupIDs[i] << " and " << sameGroupIDs[j] << std::endl;
		    int identicalMinhashes = 0;
		    int emptyBins = 0;
		    for (int k = 0; k < numBins; k++) {
		    	if (signatureMatrix[k+(sameGroupIDs[i]*numBins)] == signatureMatrix[k+(sameGroupIDs[j]*numBins)]) {
		    	    if (signatureMatrix[k+(sameGroupIDs[i]*numBins)] == INT_MAX) {
		    		emptyBins++;
		    	    } else {
		    		identicalMinhashes++;
		    	    }
		    	}
		    }
		    double similarity = (identicalMinhashes*1.0)/((numBins*1.0) - (emptyBins*1.0));
		    if (similarity >= SIMILARITY_THRESHOLD) {
//		    	std::cout << relationRSetsID[sameGroupIDs[i]] << " " << relationSSetsID[sameGroupIDs[j]-rSize] << " " << similarity << "\n";
		    	//Results
		    	result.push_back(sameGroupIDs[i]);
		    	result.push_back(sameGroupIDs[j]-rSize);
		    	similarPairs++;
		    }
		}
	    }
	}
    }
    return similarPairs;
}

void
computeSignatureMatrix(std::vector<int> &signatureMatrix, ccsMatrix* characteristicMatrix, int numSets, int numShingles, int numBins, long int currentTime)
{
    // int binSize = 65535/numBins; // Chosen for faster mod operation
  int binSize = numShingles/numBins;  
  if (numShingles % numBins) binSize++;

  // for (uint i = 0; i < signatureMatrix.size(); i++) {
  //   signatureMatrix[i] = INT_MAX;
  // }
  for (int i = 0; i < numSets; i++) {
    int offSetCM = characteristicMatrix -> col_ptr[i];
    for (int j = offSetCM; j < characteristicMatrix -> col_ptr[i+1]; j++) {
      int shingleIdx = characteristicMatrix -> row_ind[j];
      // int shingleNewIdx = fastMod(MurmurHash2(&shingleIdx, 4, currentTime));
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
  std::cout << "Not using length filtering: " << std::endl;
  similarPairsCount = nestedLoopJoin(signatureMatrix, rSize, sSize, relationRSetsID, relationSSetsID, numBins, result);


  // std::cout << "Using length filtering:" << std::endl;
  // std::unordered_map<int, std::vector<int> > candidatePairs = lengthFiltering(signatureMatrix, rSize+sSize, numBins);
  // similarPairsCount = lengthFilteringJoin(signatureMatrix, candidatePairs, rSize, sSize, relationRSetsID, relationSSetsID, numBins, result);

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
