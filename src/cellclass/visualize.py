from matplotlib import pyplot as plt
import numpy as np

def plot_acg_bar_raw_line_smooth(
    lags_ms: np.ndarray,
    acg_counts: np.ndarray,
    acg_smooth: np.ndarray,
    cell_ids: np.ndarray | None = None,
    cell_idx: int = 0,
    *,
    xlim: tuple[float, float] = (0, 50),
    show_negative: bool = False,
):
    """
    Raw ACG as bar plot, smoothed as line.
    """
    lags_ms = np.asarray(lags_ms)
    raw = np.asarray(acg_counts)[:, cell_idx]
    sm = np.asarray(acg_smooth)[:, cell_idx]

    if show_negative:
        m = (lags_ms >= -xlim[1]) & (lags_ms <= xlim[1])
    else:
        m = (lags_ms >= xlim[0]) & (lags_ms <= xlim[1])

    x = lags_ms[m]
    y = raw[m]
    ys = sm[m]

    # infer bin width from centers
    if x.size >= 2:
        bin_ms = float(np.median(np.diff(x)))
    else:
        bin_ms = 1.0

    title = f"ACG cell_idx={cell_idx}"
    if cell_ids is not None:
        title += f", cell_id={cell_ids[cell_idx]}"

    plt.figure()
    # Bars centered on bin centers, width ~ bin size
    plt.bar(x, y, width=0.9 * bin_ms, align="center", alpha=0.6, label="raw counts")
    # Smoothed line
    plt.plot(x, ys, label="smoothed", linewidth=2)

    if not show_negative:
        plt.axvline(0, linestyle="--", linewidth=1)

    plt.xlabel("Lag (ms)")
    plt.ylabel("ACG counts")
    plt.title(title)
    plt.legend()
    plt.show()