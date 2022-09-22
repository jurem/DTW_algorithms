static inline seq_t skew_mem_fw(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i, int lower) {
    seq_t s1 = &t[0];
    seq_t s2 = &t[m];
    seq_t d = &t[2*m];
    // the first & the second row
    s1[0] = DIST(a[0], b[0]);
    s2[0] = s1[0] + DIST(a[1], b[0]);
    s2[1] = s1[0] + DIST(a[0], b[1]);
    // upper triangle
    for (int i = 2; i < m; i++) {
        d[0] = DIST(a[i], b[0]) + s2[0];
        for (int j = 1; j < i; j++)
            d[j] = DIST(a[i - j], b[j]) + MIN3(s1[j - 1], s2[j - 1], s2[j]);
        d[i] = DIST(a[0], b[i]) + s2[i - 1];  // the last element
        // rotate rows
        seq_t t = s1; s1 = s2; s2 = d; d = t;
    }
    // middle part
    for (int i = m; i < end_i; i++) {
        d[0] = DIST(a[i], b[0]) + s2[0];
        for (int j = 1; j < m; j++)
            d[j] = DIST(a[i - j], b[j]) + MIN3(s1[j - 1], s2[j - 1], s2[j]);
        // rotate rows
        seq_t t = s1; s1 = s2; s2 = d; d = t;
    }
    // lower triangle
    if (lower) {
        for (int i = n; i < n + m - 1; i++) {
            for (int j = i - n + 1; j < m; j++)
                d[j] = DIST(a[i - j], b[j]) + MIN3(s1[j - 1], s2[j - 1], s2[j]);
            // rotate rows
            seq_t t = s1; s1 = s2; s2 = d; d = t;
        }
    }
    return s2;
}

static inline seq_t skew_mem_fr(seq_t a, size_t n, seq_t b, size_t m, tab_t t, size_t end_i) {
    seq_t s1 = &t[0];
    seq_t s2 = &t[m];
    seq_t d = &t[2*m];
    // the first & the second row
    s1[0] = DIST(a[n - 1], b[m - 1]);
    s2[0] = s1[0] + DIST(a[n - 2], b[m - 1]);
    s2[1] = s1[0] + DIST(a[n - 1], b[m - 2]);
    // upper triangle
    for (int i = 2; i < m; i++) {
        d[0] = DIST(a[n - 1 - i], b[m - 1]) + s2[0];
        for (int j = 1; j < i; j++)
            d[j] = DIST(a[n - 1 - i + j], b[m - 1 - j]) + MIN3(s1[j - 1], s2[j - 1], s2[j]);
        d[i] = DIST(a[n - 1], b[m - 1 - i]) + s2[i - 1];  // the last element
        // print_seq(d, m);
        // rotate rows
        seq_t t = s1; s1 = s2; s2 = d; d = t;
    }
    // middle part
    for (int i = m; i < end_i; i++) {
        d[0] = DIST(a[n - 1 - i], b[m - 1]) + s2[0];
        for (int j = 1; j < m; j++)
            d[j] = DIST(a[n - 1 - i + j], b[m - 1 - j]) + MIN3(s1[j - 1], s2[j - 1], s2[j]);
        // print_seq(d, m);
        // rotate rows
        seq_t t = s1; s1 = s2; s2 = d; d = t;
    }
    // lower triangle
    if (n <= end_i) {
        for (int i = n; i < n + m - 1; i++) {
            for (int j = i - n + 1; j < m; j++)
                d[j] = DIST(a[n - 1 - i + j], b[m - 1 - j]) + MIN3(s1[j - 1], s2[j - 1], s2[j]);
            // print_seq(d, m);
            // rotate rows
            seq_t t = s1; s1 = s2; s2 = d; d = t;
        }
    }
    return s2;
}

// ********** basic algorithms **********

val_t dtw_skew_mem_fw(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(3, m);
    seq_t s = skew_mem_fw(a, n, b, m, t, n, 1);
    val_t r = s[m - 1];
    TFREE(t);
    return r;
}

val_t dtw_skew_mem_fr(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t = TNEW(3, m);
    seq_t s = skew_mem_fr(a, n, b, m, t, n);
    val_t r = s[m - 1];
    TFREE(t);
    return r;
}

val_t dtw_skew_mem_fwfr(seq_t a, size_t n, seq_t b, size_t m) {
    tab_t t1 = TNEW(3, m);
    tab_t t2 = TNEW(3, m);
    size_t h = (n + m) / 2;
    seq_t s = skew_mem_fw(a, n, b, m, t1, h, 0);
    seq_t d = skew_mem_fr(a, n, b, m, t2, n + m - 1 - h);
    // merge
    seq_t s2;
    if (s == &t1[0]) s2 = &t1[2*m];
    else if (s == &t1[m]) s2 = &t1[0];
    else s2 = &t1[m];
    seq_t d2;
    if (d == &t2[0]) d2 = &t2[2*m];
    else if (d == &t2[m]) d2 = &t2[0];
    else d2 = &t2[m];
    // printf("-----------\n");
    // print_seq(s2, m);
    // print_seq(s, m);
    // print_seq(d, m);
    // print_seq(d2, m);
    //
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, s2[j] + d[m - 1 - j - 1]);
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, s[j] + MIN3(d[m - 1 - j], d[m - 1 - j - 1], d2[m - 1 - j - 1]));
    r = MIN(r, s[m - 1] + d[m - 1]);
    // free
    TFREE(t1);
    TFREE(t2);
    return r;
}

// ********** combined algorithms **********

void* skew_mem_fw_tophalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t s = skew_mem_fw(a, n, b, m, &t[0], h, 0);
    pthread_exit(s);
}

void* skew_mem_fr_bottomhalf(void *args) {
    DTW_DATA_GET(dtw_thread_data);
    seq_t d = skew_mem_fr(a, n, b, m, &t[3*m], n + m - 1 - h);
    pthread_exit(d);
}

val_t dtw_skew_mem_fwfr_par(seq_t a, size_t n, seq_t b, size_t m) {
    size_t h = (n + m) / 2;
    // thread data
    seq_t t = TNEW(6, m);
    DTW_DATA_SET(dtw_thread_data);
    // create and run the top and bottom thread
    pthread_create(&threads[0], NULL, skew_mem_fw_tophalf, 0);
    pthread_create(&threads[1], NULL, skew_mem_fr_bottomhalf, 0);
    seq_t s;
    seq_t d;
    pthread_join(threads[0], (void*)&s);
    pthread_join(threads[1], (void*)&d);
    // merge
    seq_t t1 = &t[0];
    seq_t t2 = &t[3*m];
    seq_t s2;
    if (s == &t1[0]) s2 = &t1[2*m];
    else if (s == &t1[m]) s2 = &t1[0];
    else s2 = &t1[m];
    seq_t d2;
    if (d == &t2[0]) d2 = &t2[2*m];
    else if (d == &t2[m]) d2 = &t2[0];
    else d2 = &t2[m];
    // printf("-----------\n");
    // print_seq(s2, m);
    // print_seq(s, m);
    // print_seq(d, m);
    // print_seq(d2, m);
    //
    val_t r = VALINF;
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, s2[j] + d[m - 1 - j - 1]);
    for (int j = 0; j < m - 1; j++)
        r = MIN(r, s[j] + MIN3(d[m - 1 - j], d[m - 1 - j - 1], d2[m - 1 - j - 1]));
    r = MIN(r, s[m - 1] + d[m - 1]);
    // free
    TFREE(t);
    return r;
}
