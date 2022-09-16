#include "common.h"
#include <unistd.h>

#define RELAX_FW(t, m, i, j)    MIN3(T(t, m, i - 1, j), T(t, m, i, j - 1), T(t, m, i - 1, j - 1))
#define RELAX_BW(t, m, i, j)    MIN3(T(t, m, i + 1, j), T(t, m, i, j + 1), T(t, m, i + 1, j + 1))
#define RELAX_FW_(i, j)         RELAX_FW(t, m, i, j)
#define RELAX_BW_(i, j)         RELAX_BW(t, m, i, j)


// ********** forward **********

val_t dtw_rect_fw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    // init first row
    T_(0, 0) = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    // calculate
    for (int i = 1; i < n; i++) {
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0); // first col
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_FW_(i, j);
    }
    val_t r = T_(n - 1, m - 1);
    TFREE(t);
    return r;
}


// ********** backward **********

val_t dtw_rect_bw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    // init last row
    T_(n - 1, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int j = m - 2; j >= 0; j--)
        T_(n - 1, j) = DIST(a[n - 1], b[j]) + T_(n - 1, j + 1);
    // calculate
    for (int i = n - 2; i >= 0; i--) {
        T_(i, m - 1) = DIST(a[i], b[m - 1]) + T_(i + 1, m - 1); // first col
        for (int j = m - 2; j >= 0; j--)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_BW_(i, j);
    }
    val_t r = T_(0, 0);
    TFREE(t);
    return r;
}


// ********** forward & backward **********

val_t dtw_rect_fwbw(seq_t a, size_t n, seq_t b, size_t m) {
// assumes: n > 1
    tab_t t = TNEW(n, m);
    size_t h = (n + 1) / 2;  // + 1 for rounding up (first half may be one line longer)
    // top half
    T_(0, 0) = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    for (int i = 1; i < h; i++) {
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
        for (int j = 1; j < m; j++)
           T_(i, j) = DIST(a[i], b[j]) + RELAX_FW_(i, j);
   }
     // bottom half
    T_(n - 1, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int j = m - 2; j >= 0; j--)
        T_(n - 1, j) = DIST(a[n - 1], b[j]) + T_(n - 1, j + 1);
    for (int i = n - 2; i >= h; i--) {
        T_(i, m - 1) = DIST(a[i], b[m - 1]) + T_(i + 1, m - 1);
        for (int j = m - 2; j >= 0; j--)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_BW_(i, j);
    }
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - 1, j) + MIN(T_(h, j), T_(h, j + 1)));
    r = MIN(r, T_(h - 1, m - 1) + T_(h, m - 1));
    TFREE(t);
    return r;
}


// ********** forward & backward in parallel **********

void* dtw_rect_fwbw_tophalf(void *args) {
    DTW_DATA_GET(args);
    T_(0, 0) = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    for (int i = 1; i < h; i++) {
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
        for (int j = 1; j < m; j++)
           T_(i, j) = DIST(a[i], b[j]) + RELAX_FW_(i, j);
    }
    pthread_exit(NULL);
}

void* dtw_rect_fwbw_bottomhalf(void *args) {
    DTW_DATA_GET(args);
    T_(n - 1, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int j = m - 2; j >= 0; j--)
        T_(n - 1, j) = DIST(a[n - 1], b[j]) + T_(n - 1, j + 1);
    for (int i = n - 2; i >= h; i--) {
        T_(i, m - 1) = DIST(a[i], b[m - 1]) + T_(i + 1, m - 1);
        for (int j = m - 2; j >= 0; j--)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_BW_(i, j);
    }
    pthread_exit(NULL);
}

val_t dtw_rect_fwbw_par(seq_t a, size_t n, seq_t b, size_t m) {
// assumes: n > 1
    tab_t t = TNEW(n, m);
    size_t h = (n + 1) / 2;
    // thread data
    DTW_DATA_SET(dtw_data);
    // create and run the top and bottom thread
    pthread_create(&threads[0], NULL, dtw_rect_fwbw_tophalf, &dtw_data);
    pthread_create(&threads[1], NULL, dtw_rect_fwbw_bottomhalf, &dtw_data);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - 1, j) + MIN(T_(h, j), T_(h, j + 1)));
    r = MIN(r, T_(h - 1, m - 1) + T_(h, m - 1));
    TFREE(t);
    return r;
}

// ********** forward with strings reversed **********

val_t dtw_rect_fwrev(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    // init first row
    T_(0, 0) = DIST(a[n - 1], b[m - 1]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[n - 1], b[m - 1 - j]) + T_(0, j - 1);
    // calculate
    for (int i = 1; i < n; i++) {
        T_(i, 0) = DIST(a[n - 1 - i], b[m - 1]) + T_(i - 1, 0); // first col
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[n - 1 - i], b[m - 1 - j]) + RELAX_FW_(i, j);

    }
    val_t r = T_(n - 1, m - 1);
    TFREE(t);
    return r;
}

// ********** forward & forward with strings reversed **********

val_t dtw_rect_fwfwrew(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    size_t h = (n + 1) / 2;  // + 1 for rounding up (first half may be one line longer)
    // top half
    T_(0, 0) = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    for (int i = 1; i < h; i++) {
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_FW_(i, j);
    }
    // bottom half
    T_(h, 0) = DIST(a[n - 1], b[m - 1]);
    for (int j = 1; j < m; j++)
        T_(h, j) = DIST(a[n - 1], b[m - 1 - j]) + T_(h, j - 1);
    for (int i = h + 1; i < n; i++) {
        T_(i, 0) = DIST(a[n - 1 + h - i], b[m - 1]) + T_(i - 1, 0);
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[n - 1 + h - i], b[m - 1 - j]) + RELAX_FW_(i, j);
    }
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - 1, j) + MIN(T_(n - 1, m - 1 - j), T_(n - 1, m - 1 - j - 1)));
    r = MIN(r, T_(h - 1, m - 1) + T_(n - 1, m - 1));
    //
    TFREE(t);
    return r;
}

// ********** forward & forward with strings reversed in parallel **********

void* dtw_rect_fwfwrev_tophalf(void *args) {
    DTW_DATA_GET(args);
    T_(0, 0) = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    for (int i = 1; i < h; i++) {
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_FW_(i, j);
    }
    pthread_exit(NULL);
}

void* dtw_rect_fwfwrev_bottomhalf(void *args) {
    DTW_DATA_GET(args);
    T_(h, 0) = DIST(a[n - 1], b[m - 1]);
    for (int j = 1; j < m; j++)
        T_(h, j) = DIST(a[n - 1], b[m - 1 - j]) + T_(h, j - 1);
    for (int i = h + 1; i < n; i++) {
        T_(i, 0) = DIST(a[n - 1 + h - i], b[m - 1]) + T_(i - 1, 0);
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[n - 1 + h - i], b[m - 1 - j]) + RELAX_FW_(i, j);
    }
    pthread_exit(NULL);
}

val_t dtw_rect_fwfwrew_par(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    size_t h = (n + 1) / 2;
    // thread data
    DTW_DATA_SET(dtw_data);
    // create and run the top and bottom thread
    pthread_create(&threads[0], NULL, dtw_rect_fwfwrev_tophalf, &dtw_data);
    pthread_create(&threads[1], NULL, dtw_rect_fwfwrev_bottomhalf, &dtw_data);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - 1, j) + MIN(T_(n - 1, m - 1 - j), T_(n - 1, m - 1 - j - 1)));
    r = MIN(r, T_(h - 1, m - 1) + T_(n - 1, m - 1));
    TFREE(t);
    return r;
}

// ********** forward parallel strides **********

void* dtw_rect_fw_strides_stride(void* args) {
    int id = (long) args;
    int from_j = stride_data[id].from_j;
    int to_j = stride_data[id].to_j;
    DTW_DATA_GET(&dtw_data);
    // calculate    
    for (int i = 1; i < n; i++) {
        // wait_for_stride_line(id, id - 1, i);
        if (id > 0) {
            pthread_mutex_lock(&stride_mutex);
            while (stride_data[id - 1].line < i) {
                stride_data[id - 1].waitline = i;
                // fprintf(stderr, "%d wait for %d\n", id, stride_data[id - 1].waitline);
                pthread_cond_wait(&stride_data[id - 1].cond, &stride_mutex);
            }
            pthread_mutex_unlock(&stride_mutex);
        }
        for (int j = from_j; j < to_j; j++)
            T_(i, j) = DIST(a[i], b[j]) + RELAX_FW_(i, j);
        if (i >= stride_data[id].waitline) {
            stride_data[id].line = i;
            // fprintf(stderr, "%d finished %d\n", id, stride_data[id].line);
            stride_data[id].waitline = INT_MAX;
            pthread_cond_signal(&stride_data[id].cond);
        }
    }
    stride_data[id].line = n - 1;
    pthread_cond_signal(&stride_data[id].cond);
    pthread_exit(NULL);
}

val_t dtw_rect_fw_strides(seq_t a, size_t n, seq_t b, size_t m) {
    int stride = m / thread_count;
    tab_t t = TNEW(n, m);
    int h = 0; // dummy
    // init thread data
    DTW_DATA_SET(dtw_data);
    for (int i = 0; i < thread_count; i++) {
        stride_data[i].from_j = 1 + i * stride;
        stride_data[i].to_j = MIN(m, 1 + (i + 1) * stride);
        stride_data[i].line = 0;
        stride_data[i].waitline = INT_MAX;
        pthread_cond_init(&stride_data[i].cond, 0);
    }
    // init first row and col
    T_(0, 0) = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    for (int i = 1; i < n; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
    // create and run threads
    pthread_attr_t attr;
    pthread_attr_init(&attr);   // privzeti atributi
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for (int i = 0; i < thread_count; i++)
        pthread_create(&threads[i], &attr, dtw_rect_fw_strides_stride, (void*) (long) i);
    for (int i = 0; i < thread_count; i++) {
        int rc = pthread_join(threads[i], NULL);
        if (rc > 0) fprintf(stderr, "error: pthread_join: %d\n", rc);
    }
    val_t r = T_(n - 1, m - 1);
    TFREE(t);
    pthread_attr_destroy(&attr);
    return r;
}

#undef RELAX_FW
#undef RELAX_BW
