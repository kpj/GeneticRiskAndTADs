#!/usr/bin/env bash

cd "$(dirname "$0")/../"
mkdir -p ./pipeline_run/cluster_logs

snakemake \
    -pr \
    -j 999 \
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
