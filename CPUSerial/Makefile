COMPILER = $(CXX)
CPPFLAGS = -Wall -O3 -std=c++11
TARGET = simJoinCPUSerial
OBJS = ../characteristicMatrix.o
ARGS = ../Datasets/trecAbstractsR10Random.txt ../Datasets/trecAbstractsS10Random.txt

.PHONY: all clean run

all: $(TARGET)

../characteristicMatrix.o: ../characteristicMatrix.cpp ../characteristicMatrix.h 
	$(COMPILER) $(CPPFLAGS) -c $< -o $@

$(TARGET): main.cpp ../miscFunctions.h $(OBJS)
	$(COMPILER) $(CPPFLAGS) $(OBJS) $< -o $@

run: $(TARGET)
	./$(TARGET) $(ARGS)

clean:
	rm -rf $(TARGET) $(OBJS)
