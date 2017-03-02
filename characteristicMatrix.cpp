#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include "characteristicMatrix.h"

ccsMatrix* 
buildCharacteristicMatrix(std::multimap<int,int> shingleSetMap)
{
  ccsMatrix *characteristicMatrix = new ccsMatrix();
  int newColumnFlag;
  int addedElements = 0;
  for (auto it = shingleSetMap.begin(); it != shingleSetMap.end();){
    //    std::cout << "*first: " << it->first << "\n";
    newColumnFlag = 1;
    std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> ii = shingleSetMap.equal_range(it->first);
    std::multimap<int,int>::iterator i;
    for (i = ii.first; i != ii.second; ++i) { //Checks what words are contained in the current set
      //      std::cout << "second: " << i -> second << "\n";
      characteristicMatrix -> row_ind.push_back(i -> second);
      if (newColumnFlag) {
	characteristicMatrix -> col_ptr.push_back(addedElements);
	newColumnFlag = 0;
      }
      addedElements++;
    }
    it = i; //Goes to the next set
  }
  characteristicMatrix -> col_ptr.push_back(addedElements);
  return characteristicMatrix;
}

void 
printCharacteristicMatrixVectors (ccsMatrix *characteristicMatrix)
{
  std::cout << "row_ind vector - Size: " << characteristicMatrix -> row_ind.size() << "\n";
  for (auto value : characteristicMatrix -> row_ind) {
    std::cout << value << " ";
  }
  std::cout << "\n\n";

  std::cout << "col_ptr vector - Size: " << characteristicMatrix -> col_ptr.size() << "\n";
  for (auto value : characteristicMatrix -> col_ptr) {
    std::cout << value << " ";
  }
  std::cout << "\n\n";
}

int
ccsMatrixSize (ccsMatrix *characteristicMatrix)
{
  return (characteristicMatrix -> row_ind.size()) + characteristicMatrix -> col_ptr.size(); 
}
