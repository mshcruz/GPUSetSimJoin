#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <vector>
#include "../characteristicMatrix.h"
#include "kernel.h"
#include "../miscFunctions.h"

int
main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << "pathToRelation/relationRFile.data pathToRelation/relationSFile.data\n";
    return 1;
  }

  Timer t; t.start();

  std::multimap<int,int> h_shingleSetMap;
  std::map<std::string,int> h_shingles;
  std::vector<std::string> relationRSetsID;
  std::vector<std::string> relationSSetsID;  

  Timer tReadFromDisk; tReadFromDisk.start();
  int rSize = processInputRelation(h_shingleSetMap, h_shingles, relationRSetsID, argv[1], 0);
  int sSize = processInputRelation(h_shingleSetMap, h_shingles, relationSSetsID, argv[2], rSize);
  tReadFromDisk.stop();

  Timer tPreProc; tPreProc.start();
  int numSets = rSize + sSize;
  int numShingles = h_shingles.size();
  int numBins = NUM_BINS;  
  ccsMatrix* h_characteristicMatrix = buildCharacteristicMatrix(h_shingleSetMap);
  tPreProc.stop();

  if (numSets != (int)num_unique_keys(h_shingleSetMap) || numSets+1 != (int)(h_characteristicMatrix -> col_ptr.size())) {
    std::cout << "Error in dataset. Possible problems: record containing only stop words or empty record." << std::endl;
    exit(-1);
  }
  
  std::vector<int> h_signatureMatrix (numSets*numBins);
  double tMinhash, tJoin, tMemoryTransfer;
  std::vector<int> h_resultPairs = kernelManager(h_signatureMatrix, h_characteristicMatrix, numShingles, rSize, sSize, numBins, tMinhash, tJoin, tMemoryTransfer);

  // Printing results
  // for (int i = 0; i < (int)h_resultPairs.size()-1; i+=2) {
  //   std::cout << relationRSetsID[h_resultPairs[i]] << " " << relationSSetsID[h_resultPairs[i+1]] << std::endl;
  // }
  // std::cout << h_resultPairs.size() / 2 << std::endl;

  delete(h_characteristicMatrix);

  t.stop();

  std::cout << t.elapsed_time() << "," << tReadFromDisk.elapsed_time() << "," << tPreProc.elapsed_time() << "," << tMinhash << "," << tJoin << "," << tMemoryTransfer << "," << h_resultPairs.size()/2 << std::endl;

  return 0;
}

