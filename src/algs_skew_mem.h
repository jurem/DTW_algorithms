static inline seq_t skew_mem_fw_old(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i, int lower) {
    seq_t x = &t[0];
    seq_t y = &t[m];
    seq_t z = &t[2*m];
    // the first & the second row
    x[0] = DIST(a[0], b[0]);
    y[0] = x[0] + DIST(a[1], b[0]);
    y[1] = x[0] + DIST(a[0], b[1]);
    // upper triangle
    for (int i = 2; i < m; i++) {
        z[0] = DIST(a[i], b[0]) + y[0];
        for (int j = 1; j < i; j++)
            z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        z[i] = DIST(a[0], b[i]) + y[i - 1];  // the last element
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    // middle part
    for (int i = m; i < end_i; i++) {
        z[0] = DIST(a[i], b[0]) + y[0];
        for (int j = 1; j < m; j++)
            z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    // lower triangle
    if (lower) {
        for (int i = n; i < n + m - 1; i++) {
            for (int j = i - n + 1; j < m; j++)
                z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
            // rotate rows
        seq_t t = x; x = y; y = z; z = t;
        }
    }
    return y;
}

static inline seq_t skew_mem_fw(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i, int lower) {
    seq_t x = &t[0];
    seq_t y = &t[m];
    seq_t z = &t[2*m];
    // the first & the second row
    x[0] = DIST(a[0], b[0]);
    y[0] = x[0] + DIST(a[1], b[0]);
    y[1] = x[0] + DIST(a[0], b[1]);
    // the rest
    for (int i = 2; i < end_i; i++) {
        if (i < n) z[0] = DIST(a[i], b[0]) + y[0];
        int left = MAX(1, i - (int)n + 1);
        int right = MIN(i, m);
        for (int j = left; j < right; j++)
            z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        if (i < m) z[i] = DIST(a[0], b[i]) + y[i - 1];  // the diagonal element
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    return y;
}

static inline seq_t skew_mem_fr(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i) {
    seq_t x = &t[0];
    seq_t y = &t[m];
    seq_t z = &t[2*m];
    // the first & the second row
    x[0] = DIST(a[n - 1], b[m - 1]);
    y[0] = x[0] + DIST(a[n - 2], b[m - 1]);
    y[1] = x[0] + DIST(a[n - 1], b[m - 2]);
    // upper triangle
    for (int i = 2; i < m; i++) {
        z[0] = DIST(a[n - 1 - i], b[m - 1]) + y[0];
        for (int j = 1; j < i; j++)
            z[j] = DIST(a[n - 1 - i + j], b[m - 1 - j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        z[i] = DIST(a[n - 1], b[m - 1 - i]) + y[i - 1];  // the last element
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    // middle part
    for (int i = m; i < end_i; i++) {
        z[0] = DIST(a[n - 1 - i], b[m - 1]) + y[0];
        for (int j = 1; j < m; j++)
            z[j] = DIST(a[n - 1 - i + j], b[m - 1 - j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        // print_seq(d, m);
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    // lower triangle
    if (n <= end_i) {
        for (int i = n; i < n + m - 1; i++) {
            for (int j = i - n + 1; j < m; j++)
                z[j] = DIST(a[n - 1 - i + j], b[m - 1 - j]) + MIN3(x[j - 1], y[j - 1], y[j]);
            // print_seq(d, m);
            // rotate rows
            seq_t t = x; x = y; y = z; z = t;
        }
    }
    return y;
}

static inline val_t skew_rev_merge(seq_t x, seq_t y, seq_t w, seq_t z, size_t m) {
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, x[j] + w[m - 1 - j - 1]);
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, y[j] + MIN3(w[m - 1 - j], w[m - 1 - j - 1], z[m - 1 - j - 1]));
    r = MIN(r, y[m - 1] + z[0]);
    return r;
}

// ********** basic algorithms **********

val_t dtw_skew_mem_fw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(3, m);
    seq_t z = skew_mem_fw(a, n, b, m, t, n + m - 1, 1);
    val_t r = z[m - 1];
    TFREE(t);
    return r;
}

val_t dtw_skew_mem_fr(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(3, m);
    seq_t z = skew_mem_fr(a, n, b, m, t, n);
    val_t r = z[m - 1];
    TFREE(t);
    return r;
}

// ********** combined algorithms **********

val_t dtw_skew_mem_fwfr(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t1 = TNEW(3, m);
    tab_t t2 = TNEW(3, m);
    size_t half = (n + m) / 2;
    seq_t y = skew_mem_fw(a, n, b, m, t1, half, 0);
    seq_t w = skew_mem_fr(a, n, b, m, t2, n + m - 1 - half);
    seq_t x = &t1[((half - 2) % 3) * m];
    seq_t z = &t2[((n + m - 1 - half - 2) % 3) * m];
    val_t r = skew_rev_merge(x, y, w, z, m);
    TFREE(t1);
    TFREE(t2);
    return r;
}

void* skew_mem_fw_tophalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t d = skew_mem_fw(a, n, b, m, &t[0], half, 0);
    pthread_exit(d);
}

void* skew_mem_fr_bottomhalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t d = skew_mem_fr(a, n, b, m, &t[3 * m], n + m - 1 - half);
    pthread_exit(d);
}

val_t dtw_skew_mem_fwfr_par(seq_t a, size_t n, seq_t b, size_t m) {
    seq_t t = TNEW(6, m);
    size_t half = (n + m) / 2;
    DTW_DATA_SET(dtw_thread_data, half, 0);
    // create and run the top and the bottom thread
    pthread_attr_t attr;
    pthread_attr_init(&attr);   // privzeti atributi
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_create(&threads[0], &attr, skew_mem_fw_tophalf, 0);
    pthread_create(&threads[1], &attr, skew_mem_fr_bottomhalf, 0);
    pthread_attr_destroy(&attr);
    // join threads and merge the results
    seq_t y;
    seq_t w;
    pthread_join(threads[0], (void*)&y);
    pthread_join(threads[1], (void*)&w);
    seq_t x = &t[((half - 2) % 3) * m];
    seq_t z = &t[(3 + (n + m - 1 - half - 2) % 3) * m];
    val_t r = skew_rev_merge(x, y, w, z, m);
    TFREE(t);
    return r;
}

// ********** strides **********

void* stride0(void* args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t x = &t[0];
    seq_t y = &t[m];
    seq_t z = &t[2*m];
    // the first & the second row
    x[0] = DIST(a[0], b[0]);
    y[0] = x[0] + DIST(a[1], b[0]);
    y[1] = x[0] + DIST(a[0], b[1]);
    // the rest
    for (int i = 2; i < half; i++) {
        if (i < n) z[0] = DIST(a[i], b[0]) + y[0];
        int left = MAX(1, i - (int)n + 1);
        int right = MIN(i, stride);
        for (int j = left; j < right; j++)
            z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        if (i < m) z[i] = DIST(a[0], b[i]) + y[i - 1];  // the diagonal element
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    return NULL;
}

void* stride1(void* args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t x = &t[0];
    seq_t y = &t[m];
    seq_t z = &t[2*m];
    // the first & the second row
    x[0] = DIST(a[0], b[0]);
    y[0] = x[0] + DIST(a[1], b[0]);
    y[1] = x[0] + DIST(a[0], b[1]);
    // the rest
    for (int i = 2; i < half; i++) {
        if (i < n) z[0] = DIST(a[i], b[0]) + y[0];
        int left = MAX(1, i - (int)n + 1);
        int right = MIN(i, m);
        for (int j = left; j < right; j++)
            z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
        if (i < m) z[i] = DIST(a[0], b[i]) + y[i - 1];  // the diagonal element
        // rotate rows
        seq_t t = x; x = y; y = z; z = t;
    }
    return NULL;
}


#include <unistd.h>

volatile int finished_count = 0;
volatile int finished = 0;

// void* t2(void* args) {
//     int id = (int)args;
//     DTW_DATA_GET(dtw_thread_data);
//     while (!finished) {
//         pthread_mutex_lock(&mem_stride_mutex);
//         while (finished_count & (1 << (id - 2)))
//             pthread_cond_wait(&mem_stride_cond, &mem_stride_mutex);
//         pthread_mutex_unlock(&mem_stride_mutex);
//         seq_t x = mem_strides[id].x;
//         seq_t y = mem_strides[id].y;
//         seq_t z = mem_strides[id].z;
//         int i = mem_strides[id].i;
//         int left = mem_strides[id].left;
//         int right = mem_strides[id].right;
//         // printf("T%d: %d %d\n", id, left, right);
//         for (int j = left; j < right; j++)
//             z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
//         pthread_mutex_lock(&mem_stride_mutex);
//             finished_count = finished_count | (0x01 << (id-2));
//             if (finished_count == 0xFF)
//                 pthread_cond_signal(&mem_stride_cond_main);            
//         pthread_mutex_unlock(&mem_stride_mutex);
//     }
//     pthread_cond_broadcast(&mem_stride_cond);
//     pthread_cond_signal(&mem_stride_cond_main);
//     pthread_exit(NULL);
// }

// static inline seq_t skew_mem_fw_s(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i, int lower) {
//     seq_t x = &t[0];
//     seq_t y = &t[m];
//     seq_t z = &t[2*m];
//     // the first & the second row
//     x[0] = DIST(a[0], b[0]);
//     y[0] = x[0] + DIST(a[1], b[0]);
//     y[1] = x[0] + DIST(a[0], b[1]);
//     pthread_attr_t attr;
//     pthread_attr_init(&attr);
//     pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
//     finished_count = 2;
//     pthread_create(&threads[2], NULL, t2, (void*)2);
//     pthread_create(&threads[3], NULL, t2, (void*)3);
//     pthread_create(&threads[4], NULL, t2, (void*)4);
//     pthread_create(&threads[5], NULL, t2, (void*)5);
//     pthread_create(&threads[6], NULL, t2, (void*)6);
//     pthread_create(&threads[7], NULL, t2, (void*)7);
//     pthread_create(&threads[8], NULL, t2, (void*)8);
//     pthread_create(&threads[9], NULL, t2, (void*)9);
//     // the rest
//     for (int i = 2; i < end_i; i++) {
//         if (i < n) z[0] = DIST(a[i], b[0]) + y[0];
//         int left = MAX(1, i - (int)n + 1);
//         int right = MIN(i, m);
//         if (right - left < 100) {
//             for (int j = left; j < right; j++)
//                 z[j] = DIST(a[i - j], b[j]) + MIN3(x[j - 1], y[j - 1], y[j]);
//         } else {
//             int h = (right - left) / 8;
//             mem_strides[2].x = x;
//             mem_strides[2].y = y;
//             mem_strides[2].z = z;
//             mem_strides[2].i = i;
//             mem_strides[2].left = left;
//             mem_strides[2].right = left + h;
//             mem_strides[3].x = x;
//             mem_strides[3].y = y;
//             mem_strides[3].z = z;
//             mem_strides[3].i = i;
//             mem_strides[3].left = left + h;
//             mem_strides[3].right = left + 2*h;
//             mem_strides[4].x = x;
//             mem_strides[4].y = y;
//             mem_strides[4].z = z;
//             mem_strides[4].i = i;
//             mem_strides[4].left = left + 2*h;
//             mem_strides[4].right = left + 3*h;
//             mem_strides[5].x = x;
//             mem_strides[5].y = y;
//             mem_strides[5].z = z;
//             mem_strides[5].i = i;
//             mem_strides[5].left = left + 3*h;
//             mem_strides[5].right = left + 4*h;
//             mem_strides[6].x = x;
//             mem_strides[6].y = y;
//             mem_strides[6].z = z;
//             mem_strides[6].i = i;
//             mem_strides[6].left = left + 4*h;
//             mem_strides[6].right = left + 5*h;
//             mem_strides[7].x = x;
//             mem_strides[7].y = y;
//             mem_strides[7].z = z;
//             mem_strides[7].i = i;
//             mem_strides[7].left = left + 5*h;
//             mem_strides[7].right = left + 6*h;
//             mem_strides[8].x = x;
//             mem_strides[8].y = y;
//             mem_strides[8].z = z;
//             mem_strides[8].i = i;
//             mem_strides[8].left = left + 6*h;
//             mem_strides[8].right = left + 7*h;
//             mem_strides[9].x = x;
//             mem_strides[9].y = y;
//             mem_strides[9].z = z;
//             mem_strides[9].i = i;
//             mem_strides[9].left = left + 7*h;
//             mem_strides[9].right = right;
//             pthread_mutex_lock(&mem_stride_mutex);
//                 finished_count = 0;
//                 pthread_cond_broadcast(&mem_stride_cond);
//                 while (finished_count != 0xFF)
//                     pthread_cond_wait(&mem_stride_cond_main, &mem_stride_mutex);
//             pthread_mutex_unlock(&mem_stride_mutex);
//         }
//         if (i < m) z[i] = DIST(a[0], b[i]) + y[i - 1];  // the diagonal element
//         // rotate rows
//         seq_t t = x; x = y; y = z; z = t;

//     }
//     finished = 1;
//     finished_count = 0;
//     pthread_cond_broadcast(&mem_stride_cond);
//     pthread_join(threads[2], NULL);
//     pthread_join(threads[3], NULL);
//     pthread_join(threads[4], NULL);
//     pthread_join(threads[5], NULL);
//     pthread_join(threads[6], NULL);
//     pthread_join(threads[7], NULL);
//     pthread_join(threads[8], NULL);
//     pthread_join(threads[9], NULL);
//     pthread_attr_destroy(&attr);
//     return y;
// }

val_t dtw_skew_mem_fw_s(seq_t a, size_t n, seq_t b, size_t m) {
    int thread_count = 2;
    int stride = m / thread_count;
    //
    tab_t t = TNEW(3, m);
    int h = n + m - 1;
    // init thread data
    DTW_DATA_SET(dtw_thread_data, h, m);
    stride0(NULL);
    //seq_t s = skew_mem_fw_s(a, n, b, m, t, n + m - 1, 1);
    seq_t s = &t[((h-1) % 3) * m];
    val_t r = s[m - 1];
    TFREE(t);
    return r;
}


