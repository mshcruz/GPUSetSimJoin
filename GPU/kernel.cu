#include "stdio.h"
#include "kernel.h"
#include "curand_kernel.h"
#include <thrust/scan.h>
#include <thrust/execution_policy.h>

#define THREADS_PER_BLOCK 256

__device__ int d_similarPairsCount(0);

__device__ unsigned int
MurmurHash2 ( const void * key, int len, unsigned int seed )
{
  // 'm' and 'r' are mixing constants generated offline.
  // They're not really 'magic', they just happen to work well.

  const unsigned int m = 0x5bd1e995;
  const int r = 24;

  // Initialize the hash to a 'random' value

  unsigned int h = seed ^ len;

  // Mix 4 bytes at a time into the hash

  const unsigned char * data = (const unsigned char *)key;

  while(len >= 4)
    {
      unsigned int k = *(unsigned int *)data;

      k *= m;
      k ^= k >> r;
      k *= m;

      h *= m;
      h ^= k;

      data += 4;
      len -= 4;
    }

  // Handle the last few bytes of the input array

  switch(len)
    {
    case 3: h ^= data[2] << 16;
    case 2: h ^= data[1] << 8;
    case 1: h ^= data[0];
      h *= m;
    };

  // Do a few final mixes of the hash to ensure the last few
  // bytes are well-incorporated.

  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;

  return h;
}

// __global__ void
// computeSignatureMatrix_kernel(int *d_signatureMatrix, int *d_cmRowIdx, int* d_cmColPtr, int numShingles, int numSets, int numBins)
// {
//   const int tid = threadIdx.x + blockDim.x * blockIdx.x;
//   int binSize = numShingles/numBins;
//   if (numShingles % numBins) binSize++;

//   if (tid < numSets) {
//     for (int i = 0; i < numBins; i++) {
//       d_signatureMatrix[i + (tid*numBins)] = INT_MAX; //NOT Coalesced
//     }
//     for (int i = d_cmColPtr[tid]; i < d_cmColPtr[tid+1]; i++) { //Coalesced
//       int shingleIdx = d_cmRowIdx[i]; //Not coalesced
//       int shingleNewIdx = MurmurHash2(&shingleIdx, 4, 0)%numShingles; //To do: remove mod
//       int binIdx = shingleNewIdx/binSize;
//       int offSetSM = binIdx + (tid*numBins);
//        d_signatureMatrix[offSetSM] = min(d_signatureMatrix[offSetSM], shingleNewIdx); //Not coalesced
//     }
//   }
// }

__global__ void
computeSignatureMatrix_kernel(int *d_signatureMatrix, int *d_cmRowIdx, int* d_cmColPtr, int numShingles, int rSize, int sSize)
{
  __shared__ int s_signatures[NUM_BINS];
  // __shared__ int s_initPosition;
  // __shared__ int s_endPosition;
  int binSize = numShingles/NUM_BINS;
  if (numShingles % NUM_BINS) binSize++;

  for (int t = blockIdx.x; t < (rSize+sSize); t += gridDim.x) { //Each block computes the signature for one set
    __syncthreads(); // Added by Kozawa
                     // Without this synchronization, s_signatures may be updated
                     // before the copy to device memory completes.
    //1 - Initialize the signatures in the shared memory
    for (int i = threadIdx.x; i < NUM_BINS; i += blockDim.x) {
      s_signatures[i] = INT_MAX;
    }
    __syncthreads();

    //2 - Get the positions of the elements for this set in the rowInd array
    int begin = d_cmColPtr[t];
    int end = d_cmColPtr[t + 1];

    /*** __syncthreads is necessary if you use shared memory. (by Kozawa) ***/
    // if (threadIdx.x == 0) {
    //   s_initPosition = d_cmColPtr[t];
    // }
    // if (threadIdx.x == 1) {
    //   s_endPosition = d_cmColPtr[t+1];
    // }
    // int numShinglesInRecord = s_endPosition - s_initPosition;

    int numShinglesInRecord = end - begin;

    //3 - Apply the hash function on the element and get its new position
    for (int i = threadIdx.x; i < numShinglesInRecord; i += blockDim.x) {
      int shingleIdx = d_cmRowIdx[begin + i]; //Coalesced
      int hashedShingleIdx = MurmurHash2(&shingleIdx, 4, 0)%numShingles;
      int binIdx = hashedShingleIdx / binSize;
      //Keep the minimum value for that position
      atomicMin(&s_signatures[binIdx], hashedShingleIdx); //Atomic
    }
    __syncthreads();

    //4 - Copy signature to signature matrix in the device memory
    for (int i = threadIdx.x; i < NUM_BINS; i += blockDim.x) {
      d_signatureMatrix[t * NUM_BINS + i] = s_signatures[i]; //Coalesced
    }
  }
}

__global__ void transpose  // M-by-N to N-by-M
(const int *src, int *dst, int M, int N)
{
  __shared__ int s_buf[32][32 + 1];
    const int nwarps = blockDim.x >> 5;
    const int warp_id = threadIdx.x >> 5;
    const int lane = threadIdx.x & 31;

    // Transpose
    int xmax = (N + 32 - 1) / 32;
    int ymax = (M + 32 - 1) / 32;
    for (int ti = blockIdx.x; ti < xmax * ymax; ti += gridDim.x) {
        int tx = ti % xmax;
        int ty = ti / xmax;
        __syncthreads();

        for (int y = warp_id; y < 32; y += nwarps) {
            int my_x = tx * 32 + lane;
            int my_y = ty * 32 + y;
            if (my_x < N) {
                s_buf[y][lane] = src[N * my_y + my_x];
            }
        }
        __syncthreads();

        for (int x = warp_id; x < 32; x += nwarps) {
            int my_x = ty * 32 + lane;
            int my_y = tx * 32 + x;
            if (my_y < N) {
                dst[M * my_y + my_x] = s_buf[lane][x];
            }
        }
    }
}


#define BATCH 4

__global__ void
nestedLoopJoinCounter_kernel(int* d_signatureMatrix, int rSize, int sSize, int numBins, int *d_resultPositions)
{
    __shared__ int s_signatures[BATCH][NUM_BINS];
    __shared__ int s_counts[BATCH];
    // rSsize needs to be divisible by BATCH
    for (int r = BATCH * blockIdx.x; r < rSize; r += BATCH * gridDim.x) {
        __syncthreads();
        if (threadIdx.x < BATCH) {
            s_counts[threadIdx.x] = 0;
        }
        for (int i = threadIdx.x; i < NUM_BINS; i += blockDim.x) {
#pragma unroll
            for (int j = 0; j < BATCH; ++j) {
	      s_signatures[j][i] = d_signatureMatrix[(r + j) * NUM_BINS + i];
            }
        }
        __syncthreads();

        int my_counts[BATCH];
#pragma unroll
        for (int i = 0; i < BATCH; ++i) {
            my_counts[i] = 0;
        }
        for (int s = threadIdx.x; s < sSize; s += blockDim.x) {
            int identicalMinhashes[BATCH];
            int emptyBins[BATCH];
            for (int i = 0; i < BATCH; ++i) {
                identicalMinhashes[i] = 0;
                emptyBins[i] = 0;
            }
            for (int i = 0; i < NUM_BINS; ++i) {
	      // int ss = d_signatureMatrix[i * (rSize + sSize) + rSize + s];
	      int ss = d_signatureMatrix[rSize * NUM_BINS + i * sSize + s];
#pragma unroll
                for (int j = 0; j < BATCH; ++j) {
                    if (s_signatures[j][i] == ss) {
                        if (s_signatures[j][i] == INT_MAX) {
                            ++emptyBins[j];
                        } else {
                            ++identicalMinhashes[j];
                        }
                    }
                }
            }
#pragma unroll
            for (int i = 0; i < BATCH; ++i) {
                float similarity = (identicalMinhashes[i] * 1.0f) / ((numBins * 1.0f) - (emptyBins[i] * 1.0f));
                if (similarity >= SIMILARITY_THRESHOLD) {
                    ++my_counts[i];
                }
            }
        }
#pragma unroll
        for (int i = 0; i < BATCH; ++i) {
            atomicAdd(&s_counts[i], my_counts[i]);
        }
        __syncthreads();

        if (threadIdx.x < BATCH) {
            atomicAdd(&d_similarPairsCount, s_counts[threadIdx.x]);
            d_resultPositions[r + threadIdx.x] = s_counts[threadIdx.x] * 2;
        }
    }
}

__global__ void
nestedLoopJoinOutputter_kernel(int* d_signatureMatrix, int rSize, int sSize, int numBins, int *d_resultPairs, int *d_resultPositions)
{
    __shared__ int s_signatures[BATCH][NUM_BINS];
    __shared__ int s_offsets[BATCH];
    for (int r = BATCH * blockIdx.x; r < rSize; r += BATCH * gridDim.x) {
        __syncthreads();
        if (threadIdx.x < BATCH) {
            s_offsets[threadIdx.x] = d_resultPositions[r + threadIdx.x];
        }
        for (int i = threadIdx.x; i < NUM_BINS; i += blockDim.x) {
#pragma unroll
            for (int j = 0; j < BATCH; ++j) {
	      // s_signatures[j][i] = d_signatureMatrix[i * (rSize + sSize) + r + j];
	      s_signatures[j][i] = d_signatureMatrix[(r + j) * NUM_BINS + i];
            }
        }
        __syncthreads();

        for (int s = threadIdx.x; s < sSize; s += blockDim.x) {
            int identicalMinhashes[BATCH];
            int emptyBins[BATCH];
#pragma unroll
            for (int i = 0; i < BATCH; ++i) {
                identicalMinhashes[i] = 0;
                emptyBins[i] = 0;
            }
            for (int i = 0; i < NUM_BINS; ++i) {
              //  int ss = d_signatureMatrix[i * (rSize + sSize) + rSize + s]
	      int ss = d_signatureMatrix[rSize * NUM_BINS + i * sSize + s];
#pragma unroll
                for (int j = 0; j < BATCH; ++j) {
                    if (s_signatures[j][i] == ss) {
                        if (s_signatures[j][i] == INT_MAX) {
                            ++emptyBins[j];
                        } else {
                            ++identicalMinhashes[j];
                        }
                    }
                }
            }
#pragma unroll
            for (int i = 0; i < BATCH; ++i) {
                float similarity = (identicalMinhashes[i] * 1.0f) / ((numBins * 1.0f) - (emptyBins[i] * 1.0f));
                if (similarity >= SIMILARITY_THRESHOLD) {
                    int offset = atomicAdd(&s_offsets[i], 2);
                    d_resultPairs[offset] = r + i;
                    d_resultPairs[offset + 1] = s;
                }
            }
        }
    }
}

std::vector<int>
kernelManager(std::vector<int> &h_signatureMatrix, ccsMatrix* h_characteristicMatrix, int numShingles, int sSize, int rSize, int numBins, double &tMinhash, double &tJoin, double &tMemoryTransfer)
{
  int numberOfBlocks, numSets = rSize + sSize, h_similarPairsCount;

  //Device variables
  int *d_signatureMatrix, *d_cmRowIdx, *d_cmColPtr, *d_resultPairs, *d_resultPositions;

  //Size of data structures
  int cmRowIdxSize = h_characteristicMatrix -> row_ind.size();
  int cmColPtrSize = h_characteristicMatrix -> col_ptr.size();
  int smSize = h_signatureMatrix.size();

  //Characteristic matrix
  std::vector<int> h_cmRowIdx = h_characteristicMatrix -> row_ind;
  std::vector<int> h_cmColPtr = h_characteristicMatrix -> col_ptr;

//Timer initialization
  cudaEvent_t startMinhash, stopMinhash, startJoin, stopJoin, startInitialTransfer, stopInitialTransfer, startFinalTransfer, stopFinalTransfer;
  float timeMinhash, timeJoin, timeInitialTransfer, timeFinalTransfer;

  //Memory Allocation and Memory transfer CPU -> GPU
  cudaEventCreate(&startInitialTransfer);
  cudaEventCreate(&stopInitialTransfer);
  cudaEventRecord(startInitialTransfer, 0);

  cudaMalloc(&d_cmRowIdx, sizeof(int) * cmRowIdxSize);
  cudaMalloc(&d_cmColPtr, sizeof(int) * cmColPtrSize);
  cudaMalloc(&d_signatureMatrix, sizeof(int) * smSize);
  cudaMalloc(&d_resultPositions, sizeof(int) * rSize);

  cudaMemcpy(d_cmRowIdx, &h_cmRowIdx[0], sizeof(int) * cmRowIdxSize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_cmColPtr, &h_cmColPtr[0], sizeof(int) * cmColPtrSize, cudaMemcpyHostToDevice);

  cudaDeviceSynchronize();
  cudaEventRecord(stopInitialTransfer, 0);
  cudaEventSynchronize(stopInitialTransfer);
  cudaEventElapsedTime(&timeInitialTransfer, startInitialTransfer, stopInitialTransfer);

  //Minhash building signature matrix
  numberOfBlocks = numSets/THREADS_PER_BLOCK;
  if (numSets % THREADS_PER_BLOCK) numberOfBlocks++;
  cudaEventCreate(&startMinhash);
  cudaEventCreate(&stopMinhash);
  cudaEventRecord(startMinhash, 0);

  //computeSignatureMatrix_kernel<<<numberOfBlocks, THREADS_PER_BLOCK>>>(d_signatureMatrix, d_cmRowIdx, d_cmColPtr, numShingles, numSets, numBins);
  computeSignatureMatrix_kernel<<<14 * 7 * 3 * 5, THREADS_PER_BLOCK>>>(d_signatureMatrix, d_cmRowIdx, d_cmColPtr, numShingles, rSize, sSize);

  // cudaMemcpy(&h_signatureMatrix[0], d_signatureMatrix, sizeof(int)*smSize, cudaMemcpyDeviceToHost);

  // for (int i = 0; i < rSize+sSize; i++) {
  //   for (int j = 0; j < NUM_BINS; j++ ) {
  //     std::cout << h_signatureMatrix[i*NUM_BINS+j] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  cudaDeviceSynchronize();
  cudaEventRecord(stopMinhash, 0);
  cudaEventSynchronize(stopMinhash);
  cudaEventElapsedTime(&timeMinhash, startMinhash, stopMinhash);
  tMinhash = timeMinhash*1e-3;

  //Delete the characteristic after the signature matrix constructed
  cudaFree(d_cmRowIdx);
  cudaFree(d_cmColPtr);

  //Nested Loop Join
  cudaEventCreate(&startJoin);
  cudaEventCreate(&stopJoin);
  cudaEventRecord(startJoin, 0);

  int *d_signatureMatrixT;
  cudaMalloc(&d_signatureMatrixT, sizeof(int) * smSize);
  cudaMemcpy(d_signatureMatrixT, d_signatureMatrix, sizeof(int) * rSize * NUM_BINS, cudaMemcpyDeviceToDevice);
  transpose<<<256, 256>>>(d_signatureMatrix + rSize * NUM_BINS, d_signatureMatrixT + rSize * NUM_BINS, sSize, NUM_BINS);
  nestedLoopJoinCounter_kernel<<<14 * 7 * 3 * 5, THREADS_PER_BLOCK>>>(d_signatureMatrixT, rSize, sSize, numBins, d_resultPositions);
  // nestedLoopJoinCounter_kernel<<<numberOfBlocks, THREADS_PER_BLOCK>>>(d_signatureMatrix, rSize, sSize, numBins, d_resultPositions);
  // std::vector<int> h_resultPositions(rSize);
  // cudaMemcpy(&h_resultPositions[0], d_resultPositions, sizeof(int)*rSize, cudaMemcpyDeviceToHost);
  // for (int i = 0; i < h_resultPositions.size(); i++) {
  //   std::cout << h_resultPositions[i] << " ";
  // }
  // std::cout << std::endl;
  // std::cout << std::endl;
  // std::cout << std::endl;

  thrust::exclusive_scan(thrust::device, d_resultPositions, d_resultPositions + rSize, d_resultPositions); // in-place scan

  // cudaMemcpy(&h_resultPositions[0], d_resultPositions, sizeof(int)*rSize, cudaMemcpyDeviceToHost);
  // for (int i = 0; i < h_resultPositions.size(); i++) {
  //   std::cout << h_resultPositions[i] << " ";
  // }
  // std::cout << std::endl;

  cudaMemcpyFromSymbol(&h_similarPairsCount, d_similarPairsCount, sizeof(int), 0, cudaMemcpyDeviceToHost);
  cudaMalloc(&d_resultPairs, sizeof(int) * h_similarPairsCount * 2);
  nestedLoopJoinOutputter_kernel<<<14 * 7 * 3 * 5, THREADS_PER_BLOCK>>>(d_signatureMatrixT, rSize, sSize, numBins, d_resultPairs, d_resultPositions);

  cudaDeviceSynchronize();
  cudaEventRecord(stopJoin, 0);
  cudaEventSynchronize(stopJoin);
  cudaEventElapsedTime(&timeJoin, startJoin, stopJoin);
  tJoin = timeJoin*1e-3;

  //Transfer memory back to the CPU and free GPU memory
  cudaEventCreate(&startFinalTransfer);
  cudaEventCreate(&stopFinalTransfer);
  cudaEventRecord(startFinalTransfer, 0);

  std::vector<int> h_resultPairs (h_similarPairsCount*2);
  cudaMemcpy(&h_resultPairs[0], d_resultPairs, sizeof(int) * h_similarPairsCount * 2, cudaMemcpyDeviceToHost);
  //  cudaMemcpy(&h_signatureMatrix[0], d_signatureMatrix, sizeof(int)*smSize, cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();
  cudaEventRecord(stopFinalTransfer, 0);
  cudaEventSynchronize(stopFinalTransfer);
  cudaEventElapsedTime(&timeFinalTransfer, startFinalTransfer, stopFinalTransfer);
  tMemoryTransfer = timeInitialTransfer*1e-3 + timeFinalTransfer*1e-3;

  cudaFree(d_signatureMatrix);
  cudaFree(d_resultPairs);
  cudaFree(d_resultPositions);

  return h_resultPairs;
}
