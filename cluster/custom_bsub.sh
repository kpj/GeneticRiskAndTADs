#!/usr/bin/env bash
# define submission function which only returns job id

if [ $# -eq 0 ]; then
    exit 1
fi

bsub_output=$(bsub "$@")
regex="Job <([0-9]+)> is submitted .*"

if [[ $bsub_output =~ $regex ]]; then
    echo "${BASH_REMATCH[1]}"
    exit $?
else
    exit 1
fi
