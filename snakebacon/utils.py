import numpy as np
import scipy.stats as stats


def suggest_accumulation_rate(chron):
    """From core age-depth data, suggest mean accumulation rate (cm/y)
    """
    # Follow's Bacon's method @ Bacon.R ln 30 - 44
    # Suggested round vals.
    sugg = np.tile([1, 2, 5], (4, 1)) * np.reshape(np.repeat([0.1, 1.0, 10, 100], 3), (4, 3))
    # Get ballpark accumulation rates, uncalibrated dates.
    ballpacc = stats.linregress(x=chron.depth, y=chron.age * 1.1).slope
    ballpacc = np.abs(sugg - ballpacc)
    sugg = sugg.flat[ballpacc.argmin()]  # Suggest rounded acc.rate with lowest abs diff.
    return sugg

