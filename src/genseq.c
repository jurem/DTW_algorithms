#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
    // args: seed start delta len
    int seed = atoi(argv[1]);
    int start = atoi(argv[2]);
    int delta = atoi(argv[3]);
    int len = atoi(argv[4]);

    srandom(seed);
    while (len-- > 0) {
        start += random() % (2 * delta + 1) - delta;
        printf("%d ", start);
    }
    printf("\n");

    return 0;
}