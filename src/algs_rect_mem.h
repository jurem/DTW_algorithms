static inline seq_t rect_mem_fw(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i) {
    seq_t s = &t[0];
    seq_t d = &t[m];
    // first row
    s[0] = DIST(a[0], b[0]);
    for (int j = 1; j < m; j++)
        s[j] = DIST(a[0], b[j]) + s[j - 1];
    // compute    
    for (int i = 1; i < end_i; i++) {
        // next row
        d[0] = DIST(a[i], b[0]) + s[0];
        for (int j = 1; j < m; j++)
            d[j] = DIST(a[i], b[j]) + MIN3(d[j - 1], s[j - 1], s[j]);
        // swap rows
        seq_t t = s; s = d; d = t;
    }
    return s;
}

static inline seq_t rect_mem_fr(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i) {
    seq_t s = &t[0];
    seq_t d = &t[m];
    // first row
    s[0] = DIST(a[n - 1], b[m - 1]);
    for (int j = 1; j < m; j++)
        s[j] = DIST(a[n - 1], b[m - 1 - j]) + s[j - 1];
    // compute    
    for (int i = 1; i < end_i; i++) {
        // next row
        d[0] = DIST(a[n - 1 - i], b[m - 1]) + s[0];
        for (int j = 1; j < m; j++)
            d[j] = DIST(a[n - 1 - i], b[m - 1 - j]) + MIN3(d[j - 1], s[j - 1], s[j]);
        // swap rows
        seq_t t = s; s = d; d = t;
    }
    return s;
}

// ********** basic algorithms **********

val_t dtw_rect_mem_fw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(2, m);
    seq_t s = rect_mem_fw(a, n, b, m, t, n);
    val_t r = s[m - 1];
    TFREE(t);
    return r;
}

val_t dtw_rect_mem_fr(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(2, m);
    seq_t s = rect_mem_fr(a, n, b, m, t, n);
    val_t r = s[m - 1];
    TFREE(t);
    return r;
}

// ********** combined algorithms **********

val_t dtw_rect_mem_fwfr(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t1 = TNEW(2, m);
    tab_t t2 = TNEW(2, m);
    size_t h = (n + 1) / 2;
    seq_t s = rect_mem_fw(a, n, b, m, t1, h);
    seq_t d = rect_mem_fr(a, n, b, m, t2, n - h);
    // merge results
    val_t r = rect_rev_merge(s, d, m);
    //
    TFREE(t1);
    TFREE(t2);
    return r;
}

void* rect_mem_fw_tophalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t d = rect_mem_fw(a, n, b, m, &t[0], half);
    pthread_exit(d);
}

void* rect_mem_fr_bottomhalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t d = rect_mem_fr(a, n, b, m, &t[2*m], n - half);
    pthread_exit(d);
}

val_t dtw_rect_mem_fwfr_par(seq_t a, size_t n, seq_t b, size_t m) {
    size_t h = (n + 1) / 2;
    // thread data
    seq_t t = TNEW(4, m);
    DTW_DATA_SET(dtw_thread_data, h, 0);
    // create and run the top and bottom thread
    pthread_create(&threads[0], NULL, rect_mem_fw_tophalf, 0);
    pthread_create(&threads[1], NULL, rect_mem_fr_bottomhalf, 0);
    seq_t s;
    seq_t d;
    pthread_join(threads[0], (void*)&s);
    pthread_join(threads[1], (void*)&d);
    // merge results
    val_t r = rect_rev_merge(s, d, m);
    //
    TFREE(t);
    return r;
}
