import numpy as np
OUTLIER_THRESH = 5.0

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


