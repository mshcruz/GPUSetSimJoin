COMPILER = $(CXX)
CPPFLAGS = -Wall -O3 -std=c++11
CUFLAGS = -O3 -use_fast_math -c #-Xptxas -arch=sm_20 #-v 
TARGET = simJoinGPU
LIBS = -L/usr/local/cuda/lib64 -lcudart
OBJS = ../characteristicMatrix.o kernel.o
ARGS = ../Datasets/trecAbstractsR10Random.txt ../Datasets/trecAbstractsS10Random.txt

.PHONY: clean all run

all: $(TARGET)

../characteristicMatrix.o: ../characteristicMatrix.cpp ../characteristicMatrix.h
	$(COMPILER) $(CPPFLAGS) -c $< -o $@

kernel.o: kernel.cu kernel.h
	nvcc $(CUFLAGS) -o $@ $<

$(TARGET): main.cpp ../miscFunctions.h $(OBJS)
	$(COMPILER) $(CPPFLAGS) $(OBJS) $< -o $@  $(LIBS)

run: $(TARGET)
	./$(TARGET) $(ARGS)

clean:
	rm -rf $(OBJS) $(TARGET)
