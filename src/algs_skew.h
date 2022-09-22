// N.B. Currently a non-working mess of several versions.

#define RELAX_FW(t, m, i, j)  MIN3(T(t, m, i - 1, j - 1), T(t, m, i - 1, j), T(t, m, i - 2, j - 1))
#define RELAX_BW(t, m, i, j)  MIN3(T(t, m, i + 1, j + 1), T(t, m, i + 1, j), T(t, m, i + 2, j + 1))
#define RELAX_FW_(i, j)       RELAX_FW(t, m, i, j)
#define RELAX_BW_(i, j)       RELAX_BW(t, m, i, j)


// ********** forward **********

val_t dtw_skew_fw(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n >= m
    tab_t t = TNEW(n + m - 1, m);
    // init first col and skewed row
    T_(0, 0) = DIST(a[0], b[0]);
    for (int i = 1; i < n; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
    for (int j = 1; j < m; j++)
        T_(j, j) = DIST(a[0], b[j]) + T_(j - 1, j - 1);
    // upper triangle
    for (int i = 2; i < m; i++)
        for (int j = 1; j < i; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    // middle part
    for (int i = m; i < n; i++)
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    // lower triangle
    for (int i = n; i < n + m - 1; i++)
        for (int j = i - n + 1; j < m; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    val_t r = T_(n + m - 2, m - 1);
    TFREE(t);
    return r;
}


// ********** backward **********

val_t dtw_skew_bw(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n >= m
    tab_t t = TNEW(n + m - 1, m);
    // init last col and skewed row
    T_(n + m - 2, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int i = n + m - 3; i >= m - 1; i--)
        T_(i, m - 1) = DIST(a[i - m + 1], b[m - 1]) + T_(i + 1, m - 1);
    for (int j = m - 2; j >= 0; j--)
        T_(n - 1 + j, j) = DIST(a[n - 1], b[j]) + T_(n + j, j + 1);
    // lower triangle
    for (int i = n + m - 4; i >= n - 1; i--)
        for (int j = i - n + 2; j <= m - 2; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    // middle part
    for (int i = n - 2; i >= m - 1; i--)
        for (int j = 0; j <= m - 2; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    // uppper triangle
    for (int i = m - 2; i >= 0; i--)
        for (int j = 0; j <= i; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    val_t r = T_(0, 0);
    TFREE(t);
    return r;
}

// ********** forward reversed **********

val_t dtw_skew_fr(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n >= m
    tab_t t = TNEW(n + m - 1, m);
    // init first col and skewed row
    T_(0, 0) = DIST(a[n - 1], b[m - 1]);
    for (int i = 1; i < n; i++)
        T_(i, 0) = DIST(a[n - 1 - i], b[m - 1]) + T_(i - 1, 0);
    for (int j = 1; j < m; j++)
        T_(j, j) = DIST(a[n - 1], b[m - 1 - j]) + T_(j - 1, j - 1);
    // upper triangle
    for (int i = 2; i < m; i++)
        for (int j = 1; j < i; j++)
            T_(i, j) = DIST(a[n - 1 - i + j], b[m - 1 - j]) + RELAX_FW_(i, j);
    // middle part
    for (int i = m; i < n; i++)
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[n - 1 - i + j], b[m - 1 - j]) + RELAX_FW_(i, j);
    // lower triangle
    for (int i = n; i < n + m - 1; i++)
        for (int j = i - n + 1; j < m; j++)
            T_(i, j) = DIST(a[n - 1 - i + j], b[m - 1 - j]) + RELAX_FW_(i, j);
    val_t r = T_(n + m - 2, m - 1);
    TFREE(t);
    return r;
}

// ********** forward & backward **********

val_t dtw_skew_fwbw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n + m - 1, m);
    size_t h = (n + m) / 2;
    // top half
    // init first col and skewed row
    T_(0, 0) = DIST(a[0], b[0]);
    for (int i = 1; i < h; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
    for (int j = 1; j < m; j++)
        T_(j, j) = DIST(a[0], b[j]) + T_(j - 1, j - 1);
    // upper triangle
    for (int i = 2; i < m; i++)
        for (int j = 1; j < i; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    // middle part
    for (int i = m; i < h; i++)
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    // bottom half
    int ofs = n == m ? 1 : 0;
    // init last col and skewed row
    T_(n + m - 2, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int i = n + m - 3; i >= h; i--)
        T_(i, m - 1) = DIST(a[i - m + 1], b[m - 1]) + T_(i + 1, m - 1);
    for (int j = m - 2; j >= ofs; j--)
        T_(n - 1 + j, j) = DIST(a[n - 1], b[j]) + T_(n + j, j + 1);
    // lower triangle
    for (int i = n + m - 4; i >= n - 1 + ofs; i--)
        for (int j = i - n + 2; j <= m - 2; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    // // middle part
    for (int i = n - 2; i >= h; i--)
        for (int j = 0; j <= m - 2; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - 2, j) + T_(h, j + 1));
    for (int j = ofs; j < m - 1; j++)
        r = MIN(r, T_(h - 1, j) + MIN3(T_(h, j), T_(h, j + 1), T_(h + 1, j + 1)));
    if (ofs)
        r = MIN(r, T_(h - 1, 0) + T_(h, 1));
    r = MIN(r, T_(h - 1, m - 1) + T_(h, m - 1));
    TFREE(t);
    return r;
}


// ********** forward & backward in parallel **********

void* dtw_skew_fwbw_par_tophalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    // init first col and skewed row
    T_(0, 0) = DIST(a[0], b[0]);
    for (int i = 1; i < h; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
    for (int j = 1; j < m; j++)
        T_(j, j) = DIST(a[0], b[j]) + T_(j - 1, j - 1);
    // upper triangle
    for (int i = 2; i < m; i++)
        for (int j = 1; j < i; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    // middle part
    for (int i = m; i < h; i++)
        for (int j = 1; j < m; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
    pthread_exit(NULL);
}

void* dtw_skew_fwbw_par_bottomhalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    int ofs = n == m ? 1 : 0;
    // init last col and skewed row
    T_(n + m - 2, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int i = n - 2; i > n - h; i--)
        T_(i + m - 1, m - 1) = DIST(a[i], b[m - 1]) + T_(i + m, m - 1);
    for (int j = m - 2; j >= ofs; j--)
        T_(n - 1 + j, j) = DIST(a[n - 1], b[j]) + T_(n + j, j + 1);
    // lower triangle
    for (int i = n + m - 4; i >= n - 1 + ofs; i--)
        for (int j = i - n + 2; j <= m - 2; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    // middle part
    for (int i = n - 2; i >= h; i--)
        for (int j = 0; j <= m - 2; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_BW_(i, j);
    pthread_exit(NULL);
}


val_t dtw_skew_fwbw_par(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n + m - 1, m);
    size_t h = (n + m) / 2;
    // thread data
    DTW_DATA_SET(dtw_thread_data);
    // create and run the top and bottom thread
    pthread_create(&threads[0], NULL, dtw_skew_fwbw_par_tophalf, 0);
    pthread_create(&threads[1], NULL, dtw_skew_fwbw_par_bottomhalf, 0);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    // merge results
    int ofs = n == m ? 1 : 0;
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - 2, j) + T_(h, j + 1));
    for (int j = ofs; j < m - 1; j++)
        r = MIN(r, T_(h - 1, j) + MIN3(T_(h, j), T_(h, j + 1), T_(h + 1, j + 1)));
    if (ofs)
        r = MIN(r, T_(h - 1, 0) + T_(h, 1));
    r = MIN(r, T_(h - 1, m - 1) + T_(h, m - 1));
    TFREE(t);
    return r;
}


// ********** forward parallel strides **********

void* dtw_skew_fw_strides_stride(void* args) {
    int id = (long) args;
    int stride = stride_data[id].stride;
    int stride_buf_size = MAX(20, stride / 5);
    int start_i = stride_data[id].start_i;
    int start_j = stride_data[id].start_j;
    DTW_DATA_GET(dtw_thread_data);
    // printf("thread %d: %d %d %d\n", id, start_i, start_j, stride);
    // upper triangle
    int i = stride_data[id].line = start_i;
    while (i < start_i + stride - 1) {
        // printf("id %d at %d\n", id, i);
        wait_for_stride_line(id, id - 1, i - 1);
        for (int j = start_j; j <= start_j + i - start_i; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
        i++;
    }
    // middle part
    while (i < start_i + n - 2) {
        wait_for_stride_line(id, id - 1, i - 1);
        for (int j = start_j; j < start_j + stride; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
        if (i >= stride_data[id].waitline + stride_buf_size) {
            stride_data[id].line = i;
            pthread_cond_signal(&stride_data[id].cond);
        }
        i++;
    }
    // lower triangle
    wait_for_stride_line(id, id - 1, i - 1);
    while (i < start_i + n - 2 + stride) {
        for (int j = start_j + i - start_i - n + 2; j < start_j + stride; j++)
            T_(i, j) = DIST(a[i - j], b[j]) + RELAX_FW_(i, j);
        stride_data[id].line = i;
        if (i >= stride_data[id].waitline + stride_buf_size) {
            stride_data[id].line = i;
            pthread_cond_signal(&stride_data[id].cond);
        }
        i++;
    }
    pthread_cond_signal(&stride_data[id].cond);
    pthread_exit(NULL);
}

val_t dtw_skew_fw_strides(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n >= m
    int stride = m / thread_count;
    tab_t t = TNEW(n + m - 1, m);
    int h = 0; // dummy
    // init thread data
    DTW_DATA_SET(dtw_thread_data);
    for (int i = 0; i < thread_count; i++) {
        stride_data[i].start_i = 2 + i * stride;
        stride_data[i].start_j = MIN(m - 1, 1 + i * stride);
        stride_data[i].stride = stride_data[i].start_j + stride < m ? stride : m - stride_data[i].start_j;
        stride_data[i].waitline = INT_MAX;
        pthread_cond_init(&stride_data[i].cond, 0);
    }
    // init first col and skewed row
    T_(0, 0) = DIST(a[0], b[0]);
    for (int i = 1; i < n; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
    for (int j = 1; j < m; j++)
        T_(j, j) = DIST(a[0], b[j]) + T_(j - 1, j - 1);
    // create and run threads
    pthread_attr_t attr;
    pthread_attr_init(&attr);   // privzeti atributi
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for (int i = 0; i < thread_count; i++)
        pthread_create(&threads[i], &attr, dtw_skew_fw_strides_stride, (void*) (long) i);
    for (int i = 0; i < thread_count; i++) {
        int rc = pthread_join(threads[i], NULL);
        if (rc > 0) fprintf(stderr, "error: pthread_join: %d\n", rc);
    }
    // print_tab(t, n + m - 1, m);
    val_t r = T_(n + m - 2, m - 1);
    TFREE(t);
    pthread_attr_destroy(&attr);
    return r;
}


#undef RELAX_FW
#undef RELAX_BW
