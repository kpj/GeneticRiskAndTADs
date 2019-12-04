#!/usr/bin/env bash

mkdir -p ./pipeline_run/cluster_logs

snakemake \
    -pr \
    -j 999 \
    --latency-wait 30 \
    --cluster-config cluster.json \
    --cluster "$(realpath ./custom_bsub.sh) \
        {cluster.extra_args} \
        -J {cluster.name} \
        -R {cluster.resources} \
        -n {cluster.nCPUs} \
        -W {cluster.time} \
        -oo {cluster.output} \
        -eo {cluster.error}" \
    --cluster-status "$(realpath ./cluster_status.py)"
