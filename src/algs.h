#ifndef ALGS_H
#define ALGS_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include "common.h"

#include "algs_rect.h"
#include "algs_diag.h"
#include "algs_skew.h"

typedef val_t (*DTWFun)(seq_t a, size_t n, seq_t b, size_t m);

typedef struct DTWInfo {
	const char* name;
	DTWFun fun;
} DTWInfo;


const DTWInfo algs[] = {
    { "rect_fw", dtw_rect_fw },
    { "rect_bw", dtw_rect_bw },
    { "rect_fwbw", dtw_rect_fwbw },
    { "rect_fwbw_par", dtw_rect_fwbw_par },
    { "rect_fwrev", dtw_rect_fwrev },
    { "rect_fwfwrev", dtw_rect_fwfwrew },
    { "rect_fwfwrev_par", dtw_rect_fwfwrew_par },
    { "rect_fw_strides", dtw_rect_fw_strides },
    { "diag_fw", dtw_diag_fw },
    { "diag_bw", dtw_diag_bw1 },
    { "diag_fwbw", dtw_diag_fwbw },
    { "diag_fwbw_par", dtw_diag_fwbw_par },
    { "skew_fw", dtw_skew_fw },
    { "skew_bw", dtw_skew_bw },
    { "skew_fwbw", dtw_skew_fwbw },
    { "skew_fwbw_par", dtw_skew_fwbw_par },
    { "skew_fw_strides", dtw_skew_fw_strides },
    { 0, 0 }
};

const char* algName(int i) {
	return algs[i].name;
}


const DTWInfo* findAlg(const char* name) {
	for (int i = 0; algs[i].name != 0; i++)
		if (strcmp(name, algs[i].name) == 0) return &algs[i];
	return 0;
}


void printAlgs(const char *prefix) {
    printf("%s", prefix);
    for (int i = 0; algName(i) != NULL; i++) printf("%s ", algName(i));
    printf("\n");
}

#endif
