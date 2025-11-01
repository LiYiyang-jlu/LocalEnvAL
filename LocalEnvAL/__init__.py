try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False