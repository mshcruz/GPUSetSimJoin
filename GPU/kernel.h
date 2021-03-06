#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include "../characteristicMatrix.h"

#define SIMILARITY_THRESHOLD 0.8
#define NUM_BINS 256

std::vector<int> kernelManager(std::vector<int> &h_signatureMatrix, ccsMatrix* h_characteristicMatrix, int numShingles, int rSize, int sSize, int numBins, double &tMinhash, double &tJoin, double &tMemoryTransfer);
