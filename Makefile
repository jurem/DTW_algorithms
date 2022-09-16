BIN = bin
SRC = src

CC = gcc
CFLAGS = -std=c99 -O3
# -fno-vectorize -fno-slp-vectorize -fno-unroll-loops
# -arch=native


HDRS = common.h algs.h algs_rect.h

all: $(BIN) $(BIN)/dtw $(BIN)/genseq

$(BIN):
	@mkdir -p $(BIN)


$(BIN)/dtw.o: $(SRC)/dtw.c $(SRC)/common.h $(SRC)/barrier.h $(SRC)/algs.h $(SRC)/algs_rect.h $(SRC)/algs_diag.h $(SRC)/algs_skew.h
	$(CC) $(CFLAGS) -c $< -o $@


$(BIN)/dtw: $(BIN)/dtw.o
	$(CC) $(CFLAGS) $< -o $@


$(BIN)/genseq: $(SRC)/genseq.c
	$(CC) $(CFLAGS) $< -o $@


# ***** other

clean:
	@rm -rf $(BIN)

@PHONY: all clean


# benchmarks

fullbench:
	bench/bench.sh

minibench:
	rm -f results-Minis/rect_{fw,bw,fwbw,fwbw_par}.txt
	bench/bench.sh
	python3 bench/plot.py results-Minis

@PHONY: bench minibench
