#!/usr/bin/env python3

"""
Adapted from https://github.com/Snakemake-Profiles/snakemake-lsf/
"""

import re
import sys
import time
import subprocess

from typing import Tuple


WAIT_BETWEEN_TRIES = 5
TRY_TIMES = 3
SUCCESS = 'success'
RUNNING = 'running'
FAILED = 'failed'
STATUS_TABLE = {
    'PEND': RUNNING,
    'RUN': RUNNING,
    'DONE': SUCCESS,
    'PSUSP': RUNNING,
    'USUSP': RUNNING,
    'SSUSP': RUNNING,
    'WAIT': RUNNING,
    'EXIT': FAILED,
    'POST_DONE': SUCCESS,
    'POST_ERR': FAILED,
    'UNKWN': RUNNING,
}


def get_status_for_snakemake(job_status: str) -> str:
    status = STATUS_TABLE.get(job_status, 'unknown')

    if status == 'unknown':
        print(
            f'Got an unknown job status {job_status}. Defaulting to \'{FAILED}\'...',
            file=sys.stderr,
        )
        status = FAILED

    return status


def query_status(jobid: int) -> Tuple[str, str]:
    cmd = f'bjobs -o \'stat\' -noheader {jobid}'
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    out, err = proc.communicate()

    return out.decode().strip(), err.decode().strip()


def query_status_failed(stderr: str, jobid: int) -> bool:
    return stderr.startswith(f'Job <{jobid}> is not found')


def main(job_str):
    # get job id
    match = re.search(r'Job <(\d+)> is submitted', job_str)
    job_id = match.group(1)

    # query status
    stdout, stderr = query_status(job_id)

    tries = 0
    while query_status_failed(stderr, job_id) and tries < TRY_TIMES:
        time.sleep(WAIT_BETWEEN_TRIES)
        stdout, stderr = query_status(job_id)
        tries += 1

    # parse status
    job_status = stdout
    status_for_snakemake = get_status_for_snakemake(job_status)

    # tell snakemake
    print(status_for_snakemake)


if __name__ == '__main__':
    main(sys.argv[1])
