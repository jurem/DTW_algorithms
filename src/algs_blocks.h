// ********** forward **********

static inline void fw_init(tab_t t, seq_t a, seq_t b, size_t n, size_t m, size_t end_i) {
    T_(0, 0) = DIST(a[0], b[0]);
    // the first row
    for (int j = 1; j < m; j++)
        T_(0, j) = DIST(a[0], b[j]) + T_(0, j - 1);
    // the first col
    for (int i = 1; i < end_i; i++)
        T_(i, 0) = DIST(a[i], b[0]) + T_(i - 1, 0);
}

static inline void fw_block(tab_t t, seq_t a, seq_t b, size_t n, size_t m, size_t start_i, size_t end_i, size_t start_j, size_t end_j) {
    for (int i = start_i; i < end_i; i++) {
        val_t * p = TPTR(t, m, i, start_j);
        for (int j = start_j; j < end_j; j++) {
            *p = DIST(a[i], b[j]) + MIN3(*(p-m-1), *(p-m), *(p-1));
            p++;
        }
    }
}

// ********** backward **********

static inline void bw_init(tab_t t, seq_t a, seq_t b, size_t n, size_t m, size_t end_i) {
    T_(n - 1, m - 1) = DIST(a[n - 1], b[m - 1]);
    // the last row
    for (int j = m - 2; j >= 0; j--)
        T_(n - 1, j) = DIST(a[n - 1], b[j]) + T_(n - 1, j + 1);
    // the last col
    for (int i = n - 2; i >= (int)end_i; i--)
        T_(i, m - 1) = DIST(a[i], b[m - 1]) + T_(i + 1, m - 1);
}

static inline void bw_block(tab_t t, seq_t a, seq_t b, size_t n, size_t m, size_t start_i, size_t end_i, size_t start_j, size_t end_j) {
    for (int i = start_i; i >= (int)end_i; i--) {
        val_t * p = TPTR(t, m, i, start_j);
        for (int j = start_j; j >= (int)end_j; j--) {
            T_(i, j) = DIST(a[i], b[j]) + MIN3(*(p+m+1), *(p+m), *(p+1));
            p--;
        }
    }
}


// ********** forward with reversed sequences **********

static inline void fr_init(tab_t t, seq_t a, seq_t b, size_t n, size_t m, size_t end_i, size_t offset_i) {
    T_(offset_i, 0) = DIST(a[n - 1], b[m - 1]);
    // the first row
    for (int j = 1; j < m; j++)
        T_(offset_i, j) = DIST(a[n - 1], b[m - 1 - j]) + T_(offset_i, j - 1);
    // the first col
    for (int i = 1; i < end_i; i++)
        T_(offset_i + i, 0) = DIST(a[n - 1 - i], b[m - 1]) + T_(offset_i + i - 1, 0);
}

static inline void fr_block(tab_t t, seq_t a, seq_t b, size_t n, size_t m, size_t start_i, size_t end_i, size_t start_j, size_t end_j, size_t offset_i) {
    for (int i = start_i; i < end_i; i++) {
        val_t * p = TPTR(t, m, offset_i + i, start_j);
        for (int j = start_j; j < end_j; j++) {
            T_(offset_i + i, j) = DIST(a[n - 1 - i], b[m - 1 - j]) + MIN3(*(p-m-1), *(p-m), *(p-1));
            p++;
        }
    }
}
