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
