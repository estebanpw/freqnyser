CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: freqgen massive-gen compare-gen consecutive-finder one-motif-finder read-massgen-file



freqgen: freqgen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm freqgen.c -lm -o $(BIN)/freqgen

massive-gen: massive-gen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm massive-gen.c -lm -o $(BIN)/massive-gen

compare-gen: compare-gen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm compare-gen.c -lm -o $(BIN)/compare-gen

consecutive-finder: consecutive-finder.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm consecutive-finder.c -lm -o $(BIN)/consecutive-finder

one-motif-finder: one-motif-finder.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm one-motif-finder.c -lm -o $(BIN)/one-motif-finder

read-massgen-file: read-massgen-file.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm read-massgen-file.c -lm -o $(BIN)/read-massgen-file

clean:
	rm -rf $(BIN)/freqgen $(BIN)/massive-gen $(BIN)/compare-gen $(BIN)/consecutive-finder $(BIN)/one-motif-finder $(BIN)/read-massgen-file

