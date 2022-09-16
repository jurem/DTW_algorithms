#include <stdlib.h>

int thread_count = 2;

#include "algs.h"
#include "barrier.h"

int main(int argc, char* argv[]) {
    const char* algname = argv[1];
    const char* filenamea = argv[2];
    const char* filenameb = argv[3];
    thread_count = argc >= 5 ? atoi(argv[4]) : 2;

    const DTWInfo* info = findAlg(algname);
    if (!info) error(1, "Invalid algorithm: '%s'", algname);

    size_t n;
    seq_t a = load_seq(filenamea, &n);
    size_t m;
    seq_t b = load_seq(filenameb, &m);
    // printf("n=%ld, m=%ld\n", n, m);
    // print_seq(a, n);
    // print_seq(b, m);

    timing_t timer;
    timer_start(&timer);
    val_t val = info->fun(a, n, b, m);
    timer_stop(&timer);
    printf("%ld %ld %"VALFMT" %ld %ld\n", n, m, val, timer.realtime, timer.cputime);
    // printf("%ld\t%ld\t%"VALFMT"\t%ld\t%ld\n", n, m, val, timer.realtime, timer.cputime);
    return 0;
}