BIN = bin
SRC = src

CC = gcc
CFLAGS = -std=c99 -lpthread -O3
#-fno-vectorize
#-fno-slp-vectorize
#-fno-unroll-loops
#-arch=native


ALG_H := algs.h algs_blocks.h algs_rect.h algs_diag.h algs_skew.h algs_rect_mem.h algs_skew_mem.h
ALG_H := $(addprefix $(SRC)/,$(ALG_H))
HDRS = common.h algs.h algs_rect.h
HDRS := $(addprefix $(SRC)/,$(HDRS))

all: $(BIN) $(BIN)/dtw $(BIN)/genseq

$(BIN):
	@mkdir -p $(BIN)


$(BIN)/dtw.o: $(SRC)/dtw.c $(HDRS) $(ALG_H)
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
