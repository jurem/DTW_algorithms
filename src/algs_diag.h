// N.B. Currently a non-working mess of several versions.

#define RELAX_FW(t, m, i, j)    MIN3(T(t, m, i - 1, j), T(t, m, i, j - 1), T(t, m, i - 1, j - 1))
#define RELAX_BW(t, m, i, j)    MIN3(T(t, m, i + 1, j), T(t, m, i, j + 1), T(t, m, i + 1, j + 1))
#define RELAX_FW_(i, j)         RELAX_FW(t, m, i, j)
#define RELAX_BW_(i, j)         RELAX_BW(t, m, i, j)


// TODO: check merging for diag_fwbw/par when n=m

// ********** forward **********

val_t dtw_diag_fw(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n <= m
    tab_t t = TNEW(n, m);
    fw_init(a, n, b, m, t, n);
    // left triangle
    for (int j = 2; j < n; j++)
        for (int i = 1; i < j; i++)
            T_(i, j - i) = DIST(a[i], b[j - i]) + RELAX_FW_(i, j - i);
    // middle part
    for (int j = n; j <= m; j++)
        for (int i = 1; i < n; i++)
            T_(i, j - i) = DIST(a[i], b[j - i]) + RELAX_FW_(i, j - i);
    // right triangle
    for (int j = n - 1; j > 1; j--)
        for (int i = 1; i < j; i++)
            T_(n - j + i, m - i) = DIST(a[n - j + i], b[m - i]) + RELAX_FW_(n - j + i, m - i);
    val_t r = T_(n - 1, m - 1);
    TFREE(t);
    return r;
}


// ********** backward **********

val_t dtw_diag_bw1(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n <= m
    tab_t t = TNEW(n, m);
    bw_init(t, a, b, n, m, 0);
    // right triangle
    for (int j = 2; j <= n - 1; j++)
        for (int i = 2; i <= j; i++)
            T_(n - (j - i + 2), m - i) = DIST(a[n - j + i - 2], b[m - i]) + RELAX_BW_(n - j + i - 2, m - i);
    // middle part
    for (int j = m - 2; j >= n - 2; j--)
        for (int i = 0; i <= n - 2; i++)
            T_(i, j - i) = DIST(a[i], b[j - i]) + RELAX_BW_(i, j - i);
    // left triangle
    for (int j = n - 3; j >= 0; j--)
        for (int i = 0; i <= j; i++)
            T_(i, j - i) = DIST(a[i], b[j - i]) + RELAX_BW_(i, j - i);
    val_t r = T_(0, 0);
    TFREE(t);
    return r;
}

val_t dtw_diag_bw(seq_t a, size_t n, seq_t b, size_t m) {
// assume: n >= m
    tab_t t = TNEW(n, m);
    bw_init(t, a, b, n, m, 0);
    // right triangle
    for (int i = 2; i <= m - 1; i++)
        for (int j = 1; j < i; j++)
            T_(n - i + j, m - j) = DIST(a[n - i + j], b[m - j]) + RELAX_BW_(n - i + j, m - j +1);
    // middle part
    for (int i = m - 1; i < n; i++)
        for (int j = 0; j < m - 1; j++)
            T_(n - i + j - 1, m - j - 1) = DIST(a[n - i + j - 1], b[m - j - 1]) + RELAX_BW_(n - i + j - 1, m - j - 1);
    // left triangle
    for (int i = m - 2; i >= 1; i--)
        for (int j = 0; j < i; j++)
            T_(i - j - 1, j) = DIST(a[i - j - 1], b[j]) + RELAX_BW_(i - j - 1, j);
    val_t r = T_(0, 0);
    // print_tab(t, n, m);
    TFREE(t);
    return r;
}


// ********** forward & backward **********

val_t dtw_diag_fwbw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    size_t h = n / 2 + (m - 1) / 2;
    // top half
    fw_init(a, n, b, m, t, h);
    // upper triangle
    for (int i = 1; i <= m - 2; i++)
        for (int j = 0; j < i; j++)
            T_(i - j, j + 1) = DIST(a[i], b[j]) + RELAX_FW_(i - j, j + 1);
    // middle part
    for (int i = m - 1; i < h; i++)
        for (int j = 0; j < m - 1; j++)
            T_(i - j, j + 1) = DIST(a[i], b[j]) + RELAX_FW_(i - j, j + 1);
    // bottom half
    bw_init(t, a, b, n, m, h);
    // lower triangle
    for (int i = 0; i <= m - 3; i++)
        for (int j = 0; j <= i; j++)
            // T_(n - i - 2 - j, m - j - 2) = RELAX_BW_(n - i - 2 - j, m - j - 2);
            T_(n - 2 - j, m - 2 - i + j) = DIST(a[i], b[j]) + RELAX_BW_(n - 2 - j, m - 2 - i + j);
    // middle part
    for (int i = m - 1; i < h; i++)
        for (int j = 0; j < m - 1; j++)
            T_(n - i + j - 1, m - j - 2) = DIST(a[i], b[j]) + RELAX_BW_(n - i + j - 1, m - j - 2);
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - j, j) + MIN(T_(h - j + 1, j), T_(h - j + 1, j + 1)));
    r = MIN(r, T_(h - m + 1, m - 1) + T_(h - m + 2, m - 1));
    // print_tab(t, n, m);
    TFREE(t);
    return r;
}


// ********** forward & backward in parallel **********

void* dtw_diag_fwbw_par_tophalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    // init first row and col
    T_(0, 0) = DIST(a[0], b[0]);
    for (int i = 1; i <= h; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    // upper triangle
    for (int i = 1; i <= m - 2; i++)
        for (int j = 0; j < i; j++)
            T_(i - j, j + 1) = DIST(a[i], b[j]) + RELAX_FW_(i - j, j + 1);
    // middle part
    for (int i = m - 1; i < h; i++)
        for (int j = 0; j < m - 1; j++)
            T_(i - j, j + 1) = DIST(a[i], b[j]) + RELAX_FW_(i - j, j + 1);
    pthread_exit(NULL);
}

void* dtw_diag_fwbw_par_bottomhalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    // bottom half
    // init last row and col
    T_(n - 1, m - 1) = DIST(a[n - 1], b[m - 1]);
    for (int i = n - 2; i > h - m + 1; i--)
        T_(i, m - 1) = DIST(a[i], b[m - 1]) + T_(i + 1, m - 1);
    for (int j = m - 2; j >= 0; j--)
        T_(n - 1, j) = DIST(a[n - 1], b[j]) + T_(n - 1, j + 1);
    // lower triangle
    for (int i = 0; i <= m - 3; i++)
        for (int j = 0; j <= i; j++)
            // T_(n - i - 2 - j, m - j - 2) = RELAX_BW_(n - i - 2 - j, m - j - 2);
            T_(n - 2 - j, m - 2 - i + j) = DIST(a[i], b[j]) + RELAX_BW_(n - 2 - j, m - 2 - i + j);
    // middle part
    for (int i = m - 1; i < h; i++)
        for (int j = 0; j < m - 1; j++)
            T_(n - i + j - 1, m - j - 2) = DIST(a[i], b[j]) + RELAX_BW_(n - i + j - 1, m - j - 2);
    pthread_exit(NULL);
}

val_t dtw_diag_fwbw_par(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(n, m);
    size_t h = n / 2 + (m - 1) / 2;
    // thread data
    DTW_DATA_SET(dtw_thread_data);
    // create and run the top and bottom thread
    pthread_create(&threads[0], NULL, dtw_diag_fwbw_par_tophalf, 0);
    pthread_create(&threads[1], NULL, dtw_diag_fwbw_par_bottomhalf, 0);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    // merge results
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, T_(h - j, j) + MIN(T_(h - j + 1, j), T_(h - j + 1, j + 1)));
    r = MIN(r, T_(h - m + 1, m - 1) + T_(h - m + 2, m - 1));
    TFREE(t);
    return r;
}

#undef RELAX_FW
#undef RELAX_BW
