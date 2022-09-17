// ********** value macros **********

//#define ABS(a)                  (((a) >= 0) ? (a) : (-(a)))
#define ABS(a)                  abs(a)
#define DIST(a, b)              ABS((a) - (b))
#define MIN(a, b)               ((a) <= (b) ? (a) : (b))
#define MIN3(a, b, c)           MIN((a), MIN((b), (c)))
#define MAX(a, b)               ((a) >= (b) ? (a) : (b))


// ********** table creation **********

#define TNEW(n, m)              (val_t*) calloc((n) * (m), sizeof(val_t))
#define TFREE(t)                free(t)


