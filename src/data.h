#ifndef DATA_H
#define DATA_H

// ********** value representation **********

#ifdef USE_FLOAT
    typedef float val_t;
    #define VALINF      INFINITY
    #define VALFMT      "f"
#else
    typedef int val_t;
    #define VALINF      INT_MAX
    #define VALFMT      "d"
#endif

// ********** value calculations **********

//#define ABS(a)                  (((a) >= 0) ? (a) : (-(a)))
#define ABS(a)                  abs(a)
#define DIST(a, b)              ABS((a) - (b))
#define MIN(a, b)               ((a) <= (b) ? (a) : (b))
#define MIN3(a, b, c)           MIN((a), MIN((b), (c)))
#define MAX(a, b)               ((a) >= (b) ? (a) : (b))

// ********** sequences and tables **********

typedef val_t* seq_t;
typedef val_t* tab_t;

#define SNEW(n)              	(val_t*) calloc((n), sizeof(val_t))
#define SFREE(t)                free(t)

#define TNEW(n, m)              (val_t*) calloc((n) * (m), sizeof(val_t))
#define TFREE(t)                free(t)

// ********** table access **********

// element of table t with dimension n x m, row i, col j
#define TPTR(t, m, i, j)        (&(t)[(i) * (m) + (j)])
// #define T(t, m, i, j)           (t)[(i) * (m) + (j)]
#define T(t, m, i, j)           (*(TPTR(t, m, i, j)))
#define TGET(t, m, i, j)        T((t), (m), (i), (j))
#define TPUT(t, m, i, j, v)     T((t), (m), (i), (j)) = (v)

// shortcuts
#define T_(i, j)                T(t, m, i, j)

// ********** printing and loading **********

void print_seq(seq_t a, size_t n) {
    for (int i = 0; i < n; i++)
        printf("%"VALFMT" ", a[i]);
    printf("\n");
}

void print_tab(tab_t t, size_t n, size_t m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf("%"VALFMT" ", T_(i, j));
        printf("\n");
    }
    printf("\n");
}

seq_t read_seq(FILE* file, size_t *length) {
    if (*length <= 0) *length = 1000;
    seq_t a = malloc(*length * sizeof(val_t));
    size_t count = 0;
    val_t val;
    while (fscanf(file, "%"VALFMT, &val) > 0) {
        if (count >= *length) {
            *length *= 2;
            a = realloc(a, *length * sizeof(val_t));
        }
        a[count++] = val;
    }
    char buf[16];
    fscanf(file, "%s", buf); // read one more char (stdin delimiter of sequences)
    *length = count;
    return realloc(a, count * sizeof(val_t));
}

seq_t load_seq(const char* filename, size_t *length) {
    FILE *file = strcmp("-", filename) == 0 ? stdin : fopen(filename, "r");
    if (!file) error(1, "Cannot open file: '%s'", filename);
    seq_t seq = read_seq(file, length);
    if (strcmp("-", filename) != 0) fclose(file);
    return seq;
}

// ********** threads **********

#define MAX_THREADS 16
pthread_t threads[MAX_THREADS];

typedef struct {
    seq_t a;
    size_t n;
    seq_t b;
    size_t m;
    tab_t t;
    size_t h;
} dtw_thread_data_t;

dtw_thread_data_t dtw_thread_data;

// shortcuts for setting and getting
#define DTW_DATA_SET(data) \
    (data).a = (a); \
    (data).n = (n); \
    (data).b = (b); \
    (data).m = (m); \
    (data).t = (t); \
    (data).h = (h);

#define DTW_DATA_GET(data) \
    tab_t t = data.t; \
    seq_t a = data.a; \
    seq_t b = data.b; \
    size_t n = data.n; \
    size_t m = data.m; \
    size_t h = data.h;

// ********** stride data **********

typedef struct {
    int from_j;
    int to_j;
    int stride;  // stride len
    int start_i;
    int start_j;
    volatile int line;      // current line
    volatile int waitline;  // which line to wait for
    pthread_cond_t cond;
} stride_data_t;

stride_data_t stride_data[MAX_THREADS];

pthread_mutex_t stride_mutex = PTHREAD_MUTEX_INITIALIZER;

//inline
void wait_for_stride_line(int id, int which, int line) {
    if (which < 0 || stride_data[which].line >= line) return;
    pthread_mutex_lock(&stride_mutex);
    while (stride_data[which].line < line) {
        // printf("wait for %d line %d (%d) \n", which, line, stride_data[which].line); fflush(stdout);
        stride_data[which].waitline = line;
        pthread_cond_wait(&stride_data[which].cond, &stride_mutex);
    }
    stride_data[which].waitline = INT_MAX;
    pthread_mutex_unlock(&stride_mutex);
}


#endif