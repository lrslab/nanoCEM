import numpy as np
OUTLIER_THRESH = 5.0
SEQ_MIN = np.array(["A"], dtype="S1").view(np.uint8)[0]
SEQ_TO_INT_ARR = np.full(26, -1, dtype=int)
SEQ_TO_INT_ARR[0] = 0
SEQ_TO_INT_ARR[2] = 1
SEQ_TO_INT_ARR[6] = 2
SEQ_TO_INT_ARR[19] = 3
MAX_POINTS_FOR_THEIL_SEN = 1000
def rescale_lstsq(dacs, levels, shift, scale):
    norm_sig = (dacs - shift) / scale
    shift_est, scale_est = np.linalg.lstsq(
        np.column_stack([np.ones_like(norm_sig), norm_sig]),
        levels,
        rcond=None,
    )[0]
    if scale_est == 0:
        return shift, scale
    new_shift = shift - (scale * shift_est / scale_est)
    new_scale = scale / scale_est
    return new_shift, new_scale

def index_from_int_kmer(int_kmer, kmer_len):
    idx = 0
    for kmer_pos in range(kmer_len):
        idx += int_kmer[kmer_len - kmer_pos - 1] * (4 ** kmer_pos)
    return idx

def seq_to_int(seq):
    """Convert string sequence to integer encoded array

    Args:
        seq (str): Nucleotide sequence

    Returns:
        np.array containing integer encoded sequence
    """
    return SEQ_TO_INT_ARR[
        np.array(list(seq), dtype="c").view(np.uint8) - SEQ_MIN
    ]

def extract_levels(seq, int_kmer_levels, kmer_len, center_idx):
    int_seq=seq_to_int(seq).astype(np.int32)
    levels = np.zeros(int_seq.shape[0], dtype=np.float32)
    levels_mv = levels
    for pos in range(int_seq.shape[0] - kmer_len + 1):
        center_pos = pos + center_idx
        kmer_idx = index_from_int_kmer(int_seq[pos: pos + kmer_len], kmer_len)
        levels_mv[center_pos] = int_kmer_levels[kmer_idx]
    return levels

def rescale(
    seq,
    dacs,
    shift,
    scale,
    seq_to_sig_map,
    dwell_filter_pctls=(10, 90),
    min_abs_level=0.2,
    edge_filter_bases=10,
    min_levels=10,
):
    """Estimate new scaling parameters base on current signal mapping to
    estimated levels.

    Args:
        levels (np.array): Estimated reference levels for each base
        dacs (np.ndarray): Unnormalized DAC signal
        shift (float): Shift from dac to normalized signal. via formula:
            norm = (dac - shift) / scale
        scale (float): Scale from dac to normalized signal
        seq_to_sig_map (np.ndarray): Position within signal array assigned
            to each base in seq
        dwell_filter_pctls (tuple): Lower and upper percentile values to
            filter short and stalled bases from estimation
        min_abs_level (float): Minimum (absolute values) level to include
            in computaiton
        edge_filter_bases (int): Number of bases at the edge of a read to
            remove from optimization computation
        min_levels (int): Minimum number of levels to perform re-scaling
    """
    levels = extract_levels(seq)
    with np.errstate(invalid="ignore"):
        dacs_cumsum = np.empty(dacs.size + 1)
        dacs_cumsum[0] = 0
        dacs_cumsum[1:] = np.cumsum(dacs)
        dwells = np.diff(seq_to_sig_map)
        dac_means = np.diff(dacs_cumsum[seq_to_sig_map]) / dwells

    # filter bases base on lower and upper percentile of dwell to remove
    # regions of poor signal assignment
    # estimate dwell limits (these limits are treated as exclusive
    # boundaries for filtering as the lower boundary may be 1)
    dwells = np.diff(seq_to_sig_map)
    dwell_min, dwell_max = np.percentile(dwells, dwell_filter_pctls)
    edge_filter = np.full(dwells.size, True, dtype=np.bool)
    if edge_filter_bases > 0:
        edge_filter[:edge_filter_bases] = False
        edge_filter[-edge_filter_bases:] = False
    valid_bases = np.logical_and.reduce(
        (
            dwells > dwell_min,
            dwells < dwell_max,
            np.abs(levels - np.mean(levels)) > min_abs_level,
            np.logical_not(np.isnan(dac_means)),
            edge_filter,
        )
    )
    filt_levels = levels[valid_bases]
    filt_dacs = dac_means[valid_bases]
    if filt_levels.size < min_levels:
        raise Exception("Too few positions")

    return rescale_lstsq(filt_dacs, filt_levels, shift, scale)

def normalize_signal_with_lim(signal, lower_lim=-OUTLIER_THRESH, upper_lim=OUTLIER_THRESH):
    shift = np.median(signal)
    # note factor to scale MAD to SD
    scale = np.median(np.absolute(signal - shift)) * 1.4826
    norm_signal = (signal - shift) / scale

    # norm_signal = c_apply_outlier_thresh(
    #         norm_signal, lower_lim, upper_lim)

    norm_signal = np.where(norm_signal > upper_lim, upper_lim, norm_signal)
    norm_signal = np.where(norm_signal < lower_lim, lower_lim, norm_signal)
    return norm_signal

def normalize_signal(signal):
    shift = np.median(signal)
    # note factor to scale MAD to SD
    scale = np.median(np.absolute(signal - shift)) * 1.4826
    norm_signal = (signal - shift) / scale

    # norm_signal = c_apply_outlier_thresh(
    #         norm_signal, lower_lim, upper_lim)
    return norm_signal


