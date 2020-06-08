import argparse
import logging
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyranges

BIN_SIZE = 100
OUTPUT_FILE = f'scores_by_virus_anova_pvals_{BIN_SIZE}.csv'
MAX_LENGTH = 31686
PALS_DIR = './pal'
COLUMNS = [
    'Start',
    'End',
    'Stem_len',
    'Loop_len',
    'Stem1',
    'Stem2',
    'Loop',
    'Representation',
    'Full_sequence',
]

logger = logging.getLogger('overlaps_by_virus_anova')
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


def get_bins_borders(bin_size, chr_size):
    return np.arange(0, chr_size + bin_size - 1, bin_size / 2, dtype=int,)  # Include right border


def get_ranges_scores_by_virus(inputs):
    '''
    Split each virus into bins of fixed length and count intersections between stem-loop ranges
    and each bin range.
    '''

    virus, pal_df = inputs
    logger.info(f'PID: {os.getpid()} - Processing {virus}.')
    pal_df['Chromosome'] = np.repeat(0, pal_df.shape[0])  # Workaround for non-chromosome data
    ranges = pyranges.PyRanges(pal_df)
    bin_borders = get_bins_borders(BIN_SIZE, 31686)
    bin_ranges = list(zip(bin_borders[:-2], bin_borders[2:]))
    bins_stats = []
    for range_ in bin_ranges:
        intersections_df = ranges.intersect(
            pyranges.PyRanges(
                starts=[range_[0]], ends=[range_[1]], chromosomes=[0]
            )  # Chromosomes param added as a workaround for non-chromosome data
        ).df
        if intersections_df.empty:
            bins_stats.append(dict(virus=virus, Start=range_[0], End=range_[1], Score=0,))
            continue
        score = intersections_df.shape[0]
        bins_stats.append(dict(Start=range_[0], End=range_[1], virus=virus, Score=score,))
    result_df = pd.DataFrame(bins_stats)
    return result_df


def main(processes):
    '''
    '''
    logger.info(f'Running script with up to {processes} processes.')
    logger.info('Reading input data.')

    filenames = os.listdir(PALS_DIR)
    input_data = {
        filename: pd.read_csv(os.path.join(PALS_DIR, filename), sep='\t', header=None).rename(
            columns=dict(enumerate(COLUMNS))
        )
        for filename in filenames
    }

    logger.info('Computing bins pvals.')
    with Pool(processes=min(processes, len(filenames))) as pool:
        results = pool.map(
            get_ranges_scores_by_virus, ((virus, pal_df) for virus, pal_df in input_data.items()),
        )

    logger.info('Concatenating result.')

    result_df = pd.concat(results, axis=0)

    logger.info('Writing result to file.')
    result_df.to_csv(OUTPUT_FILE)

    logger.info(f'Finished computing. Outputfile is {OUTPUT_FILE}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--processes", type=int, help="Count of processes to run in parallel.", default=os.cpu_count()
    )
    args = parser.parse_args()
    main(args.processes)
