#!/usr/bin/env bash

cd "$(dirname "$0")/../"

workdir="$(grep ^workdir config.yaml | cut -d ' ' -f 2 | tr -d "'")" # oh no, very bad
mkdir -p "./$workdir/cluster_logs"

snakemake \
    -pr \
    --use-conda \
    --restart-times 3 \
    --cores 100 \
    --local-cores 1 \
    --latency-wait 30 \
    --cluster-config "./cluster/cluster.json" \
    --cluster "$(realpath ./cluster/custom_bsub.sh) \
        {cluster.extra_args} \
        -J {cluster.name} \
        -R 'rusage[mem={cluster.resources}]' \
        -n {cluster.nCPUs} \
        -W {cluster.time} \
        -oo {cluster.output} \
        -eo {cluster.error}" \
    --cluster-status "$(realpath ./cluster/cluster_status.py)" \
    "$@"
