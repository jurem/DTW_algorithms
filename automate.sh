#!/bin/bash
#SBATCH --job-name=dtw
#SBATCH --output=dtw.out
#SBATCH --ntasks=5
#SBATCH --nodes=1
# --- SaBATCH --time:00:01:00

function is_sling {
    which sinfo >/dev/null
}

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
    echo 'Benchmarking ...'
    for alg in $algs; do
        local results="$DST/$alg.csv"
        echo -n "Executing $alg ..."
        $OVERWRITE || { test -f "$results" && echo " skipped" && continue; }
        echo
        # <full_name>=<algorithm>-<optional_argument>
        arg=${alg##*-}
        alg=${alg%%-*}
        {
            printf "%s %s %s %s %s\n" n m value realtime cputime
            for len in $lens; do
                if is_sling; then
                    srun --ntasks=1 bin/dtw $alg "$DST/a-$len.in" "$DST/b-$len.in" $arg
                else
                    bin/dtw $alg "$DST/a-$len.in" "$DST/b-$len.in" $arg
                fi
            done
        } | tee "$results"
    done
    #
    rm -f "$DST/"{a,b}*.in
}

function show_help {
    echo "Run as: $(basename $0) <lengths> <algorithms> <options>"
    echo
    echo '  Lengths:'
    echo '    small halfmedium medium big'
    echo '  Algorithms:'
    echo "    $RECT"
    echo "    $DIAG"
    echo "    $SKEW"
    exit 0
}

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

# Settings for sequence generator
SEED=42
START=100
DELTA=2

# Algorithms
FAMILIES="rect,diag,skew"
STRIDES="2,4,8"
#
RECT_BASE=$(prefix rect_ fw,bw,fr)
RECT_COMB=$(prefix rect_ fwbw,fwfr,fwbw_par,fwfr_par)
RECT_STRIDES=$(prefix rect_fw_strides- $STRIDES)
RECT="$RECT_BASE $RECT_COMB $RECT_STRIDES"
#
DIAG_BASE=$(prefix diag_ fw,bw)
DIAG_COMB=$(prefix diag_ fwbw,fwbw_par)
DIAG="$DIAG_BASE $DIAG_COMB"
#
SKEW_BASE=$(prefix skew_ fw,bw)
SKEW_COMB=$(prefix skew_ fwbw,fwbw_par)
SKEW_STRIDES=$(prefix skew_fw_strides- $STRIDES)
SKEW="$SKEW_BASE $SKEW_COMB $SKEW_STRIDES"
#
FW=$(postfix _fw rect,diag,skew)
BW=$(postfix _bw rect,diag,skew)
FWBW=$(postfix _fwbw rect,diag,skew)
PAR=$(postfix _fwbw_par rect,diag,skew)
#
ALLALGS="$RECT $DIAG $SKEW"

# Options
COMPILE=false
BENCH=false
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
        fig1)
            algs="$algs rect_fw rect_bw rect_fwbw rect_fwbw rect_fwbw_par" ;;
        optimize)
            algs="$algs rect_fw rect_fw1" ;;            
        # destination folder
        results-*)
            DST="$1" ;;
        # options
        compile) COMPILE=true ;;
        bench) BENCH=true ;;
        overwrite) OVERWRITE=true ;;
        plot) PLOT=true ;;
        #
        help) show_help ;;
        *) err "Invalid argument '$1'" ;;
    esac
    shift
done

# Sort lens and algs
lens=$(printf "%s\n" $lens | sort -g | uniq)
algs=$(printf "%s\n" $algs | sort | uniq | tr "\n" ' ')

# Set destination folder
test -z "$DST" && DST=results-$(hostname)
mkdir -p $DST

# Print settings
echo Selected options:
echo "  Destination folder: $DST"
echo "  Lengths: $(echo $lens | tr "\n" ' ')"
echo "  Algorithms: $algs"

# compile
test -d bin || COMPILE=true
if $COMPILE; then
    echo 'Compiling ...'
    make clean && make || err "Compile failed"
fi

# bench
if $BENCH; then
    bench_algs "$algs" "$lens"
fi

# plot
if $PLOT; then
    echo 'Plotting ...' &
    python3 plot.py $DST $algs
fi
