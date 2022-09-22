#!/bin/bash
#SBATCH --job-name=dtw
#SBATCH --output=dtw.out
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --exclusive

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

function medium {
    seq 10000 5000 50000
}

function big {
    seq 50000 10000 100000
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

function check_algs {
    local algs="$1"
    local len="$2"
    local len2=$((len/2))
    bin/genseq $SEED $START $DELTA $len >"$DST/a-$len.in"
    bin/genseq $((SEED+1)) $START $DELTA $len2 >"$DST/b-$len2.in"
    echo 'Checking ...'
    #
    for alg in $algs; do
        arg=${alg##*-}
        alg=${alg%%-*}
        [[ "$arg" == "$alg" ]] && unset arg
        local res=$(bin/dtw $alg "$DST/a-$len.in" "$DST/b-$len2.in" $arg)
        echo "$res $alg ($arg)"
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
RECT="$RECT_BASE $RECT_COMB"
RECT_MEM_BASE=$(prefix rect_mem_ fw,fr)
RECT_MEM_COMB=$(prefix rect_mem_ fwfr,fwfr_par)
RECT_MEM="$RECT_MEM_BASE $RECT_MEM_COMB"
RECT_PAR=$(prefix rect_ fwbw_par,fwfr_par,mem_fwfr_par)
#
RECT_PAR=$(prefix rect_ fwbw_par,fwfr_par,mem_fwfr_par)
RECT_STRIDES=$(prefix rect_fw_strides- $STRIDES)
RECT_ALL="$RECT $RECT_MEM"
#
DIAG_BASE=$(prefix diag_ fw,bw)
DIAG_COMB=$(prefix diag_ fwbw,fwbw_par)
DIAG_ALL="$DIAG_BASE $DIAG_COMB"
#
SKEW_BASE=$(prefix skew_ fw,bw,fr)
SKEW_COMB=$(prefix skew_ fwbw,fwbw_par)
SKEW="$SKEW_BASE $SKEW_COMB"
SKEW_MEM_BASE=$(prefix skew_mem_ fw,fr)
SKEW_MEM_COMB=$(prefix skew_mem_ fwfr)
SKEW_MEM="$SKEW_MEM_BASE $SKEW_MEM_COMB"
#
SKEW_PAR=$(prefix skew_ fwbw_par,mem_fwfr_par)
SKEW_STRIDES=$(prefix skew_fw_strides- $STRIDES)
SKEW_ALL="$SKEW $SKEW_MEM $SKEW_STRIDES"
#
FW=$(postfix _fw rect,skew)
BW=$(postfix _bw rect,skew)
FWBW=$(postfix _fwbw rect,skew)
PAR=$(postfix _fwbw_par rect,skew)
STRIDES="$RECT_STRIDES $SKEW_STRIDES"
MEM="$RECT_MEM $SKEW_MEM"
#
ALL="$RECT_ALL $DIAG_ALL $SKEW_ALL"

# Options
COMPILE=false
CHECK=false
BENCH=false
OVERWRITE=false
PLOT=false
PLOTACC=false
while test -n "$1"; do
    case $1 in
        # determine sizes of inputs
        small|halfmedium|medium|big)
            lens="$lens $($1)"
        ;;
        # determine algorithms
        rect_base|rect_comb|rect|\
        rect_mem_base|rect_mem_comb|rect_mem|\
        rect_strides|rect_par|rect_all|\
        diag|\
        skew_base|skew_comb|skew|\
        skew_mem_base|skew_mem_comb|skew_mem|\
        skew_strides|skew_par|skew_all|\
        fw|bw|fwbw|par|strides|mem|\
        all)
            var=$(echo $1 | tr "a-z" "A-Z")
            algs="$algs ${!var}"
        ;;

        alg-*)
            algs="$algs ${1#alg-}"
        ;;
        best)
            algs="$algs skew_fwbw_par"
        ;;    
        fig1)
            algs="$algs rect_fw rect_bw rect_fwbw rect_fwbw rect_fwbw_par" ;;
        optimize)
            algs="$algs rect_fw rect_fw1"
        ;;
        # destination folder
        results-*)
            DST="$1" ;;
        # options
        compile) COMPILE=true ;;
        check) CHECK=true ;;
        bench) BENCH=true ;;
        overwrite) OVERWRITE=true ;;
        plot) PLOT=true ;;
        plotacc) PLOTACC=true ;;
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

# check
if $CHECK; then
    check_algs "$algs" 1000
fi

# bench
if $BENCH; then
    bench_algs "$algs" "$lens"
fi

# plot
if $PLOT; then
    echo 'Plotting ...' &
    python3 plot.py $DST time $algs
elif $PLOTACC; then
    echo 'Plotting ...' &
    python3 plot.py $DST acc $algs
fi
