import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

def clean_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    keep = np.isfinite(x) & np.isfinite(y)

    return x[keep], y[keep]

def first_value(result):
    """
    Handles both older and newer SciPy return formats.
    Older SciPy returns a tuple: (statistic, pvalue)
    Newer SciPy returns an object with .statistic
    """
    if hasattr(result, "statistic"):
        return result.statistic

    return result[0]

def correlation_metrics(proxy, observed):
    proxy, observed = clean_pair(proxy, observed)

    if len(proxy) < 3:
        return {
            "n": len(proxy),
            "pearson": np.nan,
            "spearman": np.nan,
        }

    pearson_value = first_value(pearsonr(proxy, observed))
    spearman_value = first_value(spearmanr(proxy, observed))

    return {
        "n": len(proxy),
        "pearson": float(pearson_value),
        "spearman": float(spearman_value),
    }

def polarity_agreement(predicted, observed):
    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)

    keep = np.isfinite(predicted) & np.isfinite(observed)
    predicted = predicted[keep]
    observed = observed[keep]

    if len(predicted) == 0:
        return {
            "n": 0,
            "agreement_fraction": np.nan,
        }

    predicted_sign = np.sign(predicted)
    observed_sign = np.sign(observed)

    nonzero = (predicted_sign != 0) & (observed_sign != 0)

    if nonzero.sum() == 0:
        return {
            "n": 0,
            "agreement_fraction": np.nan,
        }

    agreement = predicted_sign[nonzero] == observed_sign[nonzero]

    return {
        "n": int(nonzero.sum()),
        "agreement_fraction": float(agreement.mean()),
    }

def high_speed_arrival_error(predicted_times, observed_times):
    predicted_times = pd.to_datetime(predicted_times)
    observed_times = pd.to_datetime(observed_times)

    if len(predicted_times) == 0 or len(observed_times) == 0:
        return np.nan

    errors_hours = []

    for predicted in predicted_times:
        deltas = np.abs((observed_times - predicted).total_seconds()) / 3600
        errors_hours.append(deltas.min())

    return float(np.mean(errors_hours))

# Synthetic test data.
# This is not science output.
proxy_speed = np.array([350, 380, 420, 500, 620, 590, 450, np.nan, 390])
omni_speed = np.array([340, 400, 430, 520, 610, 570, 460, 410, 395])

predicted_polarity = np.array([1, 1, -1, -1, 1, 1, np.nan])
observed_polarity = np.array([1, -1, -1, -1, 1, 1, 1])

predicted_streams = ["2024-08-03 12:00", "2024-08-11 06:00"]
observed_streams = ["2024-08-03 18:00", "2024-08-12 00:00"]

corr = correlation_metrics(proxy_speed, omni_speed)
pol = polarity_agreement(predicted_polarity, observed_polarity)
timing = high_speed_arrival_error(predicted_streams, observed_streams)

print("Metric sanity test")
print()
print("Speed proxy correlation:")
print(corr)
print()
print("Polarity agreement:")
print(pol)
print()
print("Mean high-speed-stream timing error in hours:")
print(timing)

assert corr["n"] == 8
assert np.isfinite(corr["pearson"])
assert pol["n"] == 6
assert np.isfinite(pol["agreement_fraction"])
assert np.isfinite(timing)

print()
print("Status: pass")
