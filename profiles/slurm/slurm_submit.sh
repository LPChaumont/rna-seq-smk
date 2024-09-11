#!/usr/bin/env bash

# https://github.com/snakemake/snakemake/issues/2494#issuecomment-1783711397

sbatch_params=$1
rule=$2
wildcards=$3
dependencies=$4
script=$5

#wildcards=${wildcards//\//--}
#logs="--output=logs/${rule}/${rule}-${wildcards}-%j.out --error=logs/${rule}/${rule}-${wildcards}-%j.err"

if [[ ! -z $dependencies ]]
then
    depend="--depend="
    for dependency in $dependencies
    do
        depend="$depend"afterok:"$dependency"","
    done
    depend=${depend%,}
fi

#echo "$dependencies" >&2

sbr="$(sbatch $sbatch_params $depend $script )"
if [[ "$sbr" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]}"
    exit 0
else
    echo "sbatch failed"
    exit 1
fi
