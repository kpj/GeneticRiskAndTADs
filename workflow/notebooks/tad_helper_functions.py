import sys
import collections

import numpy as np
import pandas as pd

from tqdm.auto import tqdm


class EmptyTAD(Exception):
    pass


class TADTooSmall(Exception):
    pass


class OverlappingTADs(Exception):
    pass


class RangeDict(dict):
    """Optimized for ranges with step==1."""
    def __getitem__(self, item):
        if type(item) != range:
            for key in self:
                if key.step == 1:
                    if key.start <= item < key.stop:
                        return self[key]
                else:
                    if item in key:
                        return self[key]
            raise KeyError(item)
        else:
            return super().__getitem__(item)


def get_tad_lengths(row, border_range):
    """Get TAD and border lengths depending on type."""
    tad_start = row.tad_start
    tad_stop = row.tad_stop

    bp_in, bp_out = border_range

    # rescale border length if TAD is too small
    tad_len = tad_stop - tad_start
    if tad_len <= 0:
        raise EmptyTAD(row)

    if tad_len < 2*bp_in:
        bp_in = tad_len // 4
        # raise TADTooSmall(row)

    # assert that TADs are not overlapping
    # rescale borders if they would overlap
    dist_prev = (tad_start - row.prev_tad_stop) \
        if (
            not np.isnan(row.prev_tad_stop)
            and row.chrname == row.prev_tad_chr
        ) else float('inf')
    dist_next = (row.next_tad_start - tad_stop) \
        if (
            not np.isnan(row.next_tad_start)
            and row.chrname == row.next_tad_chr
        ) else float('inf')
    dist_min = min(dist_prev, dist_next)
    if dist_min < 0:
        raise OverlappingTADs(row)

    if dist_min < 2*bp_out:
        bp_out = dist_min // 4

    # final sanity checks
    bp_in = int(bp_in)
    bp_out = int(bp_out)
    assert bp_in >= 0 and bp_out >= 0

    # return appropriate ranges
    return (
        range(tad_start - bp_out, tad_start + bp_in),
        range(tad_start + bp_in, tad_stop - bp_in),
        range(tad_stop - bp_in, tad_stop + bp_out)
    )


def parse_tad_annotations(border_range, fname):
    print(f' > Parsing TADs ({border_range})', file=sys.stderr)
    df_tad = pd.read_csv(fname)
    df_tad['prev_tad_stop'] = df_tad.tad_stop.shift(1)
    df_tad['next_tad_start'] = df_tad.tad_start.shift(-1)
    df_tad['prev_tad_chr'] = df_tad.chrname.shift(1)
    df_tad['next_tad_chr'] = df_tad.chrname.shift(-1)

    error_counter = collections.defaultdict(int)
    res = collections.defaultdict(RangeDict)

    tqdm.pandas()
    for row in tqdm(df_tad.itertuples(), total=df_tad.shape[0]):
        try:
            rb1, rt, rb2 = get_tad_lengths(row, border_range)
        except EmptyTAD as ex:
            error_counter['empty_tad'] += 1
            continue
        except TADTooSmall as ex:
            error_counter['small_tad'] += 1
            continue
        except OverlappingTADs as ex:
            error_counter['overlapping_tads'] += 1
            continue

        # normalize chromosome name
        chrom = row.chrname
        if chrom.startswith('chr'):
            chrom = chrom[3:]

        # store range-associations
        res[chrom][rb1] = 'border'
        res[chrom][rt] = 'tad'
        res[chrom][rb2] = 'border'

    if error_counter:
        print('TAD errors:')
        for k, v in sorted(error_counter.items()):
            print(f' > {k}: {v}')

    return dict(res)
