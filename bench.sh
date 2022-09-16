#!/bin/bash

function err {
    echo $1 >&1
    exit 42
}

function small {
    seq 1000 1000 10000
}

function halfmedium {
    seq 10000 10000 50000
}

function medium {
    seq 10000 10000 100000
}

function big {
    seq 100000 100000 20000
}

function bench_algs {
    local algs="$1"
    local lens="$2"
    #
    echo 'Generating inputs ...'
    for len in $lens; do
        bin/genseq $SEED $START $DELTA $((len-0)) >"$DST/a-$len.in"
        bin/genseq $((SEED+1)) $START $DELTA $len >"$DST/b-$len.in"
    done
    #
    for alg in $algs; do
        local results="$DST/$alg.csv"
        $OVERWRITE || { test -f "$results" && continue; }
        # <full_name>=<algorithm>-<optional_argument>
        arg=${alg##*-}
        alg=${alg%%-*}
        test -n "$arg" && local prefix=" with argument $arg"
        echo "Benchmarking $alg$prefix"
        {
            printf "%s %s %s %s %s\n" n m value realtime cputime
            for len in $lens; do
                bin/dtw $alg "$DST/a-$len.in" "$DST/b-$len.in" $thread_count
            done
        } | tee "$results"
    done
    #
    rm -f "$DST/"{a,b}*.in
}

function show_help {
    echo "Run as: $(basename $0) <lengths> <algorithms> <options>"
    echo

    echo Algorithms
    echo "  $RECT"
}

# Settings for sequence generator
SEED=42
START=100
DELTA=2

# Destination folder
DST=results-$(hostname)
mkdir -p $DST


function prefix {
    local p=$1
    shift
    eval "echo $p{$*}"
}

function postfix {
    local p=$1
    shift
    eval "echo {$*}$p"
}


# Algorithms
FAMILIES="rect,diag,skew"
STRIDES="2,4,8"
#
RECT_BASE=$(prefix rect_ fw,bw,fwrev)
RECT_COMB=$(prefix rect_ fwbw,fwbw_par,fwfwrev,fwfwrev_par)
RECT_STRIDES=$(prefix rect_fw_strides- $STRIDES)
RECT="$RECT_BASE $RECT_COMB"
#
DIAG=$(prefix diag_ fw,bw,fwbw,fwbw_par)
SKEW=$(prefix skew_ fw,bw,fwbw,fwbw_par)
#
FW=$(postfix _fw rect,diag,skew)
BW=$(postfix _bw rect,diag,skew)
FWBW=$(postfix _fwbw rect,diag,skew)
PAR=$(postfix _fwbw_par rect,diag,skew)

# Options
COMPILE=false
OVERWRITE=false
PLOT=false
while test -n "$1"; do
    case $1 in
        # determine sizes of inputs
        small|halfmedium|medium|big)
            lens="$lens $($1)"
        ;;
        # determine algorithms
        rect_base|rect_comb|rect_strides|rect|\
        diag|\
        skew|skew_strides|\
        fw|bw|fwbw|par)
            var=$(echo $1 | tr "a-z" "A-Z")
            algs="$algs ${!var}"
        ;;
        # options
        compile) COMPILE=true ;;
        overwrite) OVERWRITE=true ;;
        plot) PLOT=true ;;
        #
        help) show_help ;;
        *) err "Invalid argument '$1'" ;;
    esac
    shift
done
# sort lens
lens=$(printf "%s\n" $lens | sort -g | uniq)


# compile
test -d bin || COMPILE=true
if $COMPILE; then
    echo 'Compiling ...'
    make clean && make || err "Compile failed"
fi

# bench
bench_algs "$algs" "$lens"

# plot
if $PLOT; then
    echo 'Plotting ...' &
    python3 plot.py $DST $algs
fi
