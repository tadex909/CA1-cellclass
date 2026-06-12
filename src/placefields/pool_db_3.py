import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, kruskal, mannwhitneyu

# -----------------------------
# Fonctions utilitaires
# -----------------------------
class TSD:
    """Classe pour stocker des séries temporelles avec métadonnées."""
    def __init__(self, tvec: np.ndarray, data: np.ndarray):
        self.tvec = tvec
        self.data = data
        self.usr: Dict = {}
        self.cfg: Dict = {'history': {'func_name': [], 'cfg': []}}
def group_age(age):
    if 15 <= age <= 17:
        return "15-17"
    elif 18 <= age <= 20:
        return "18-20"
    elif 21 <= age <= 22:
        return "21-22"
    elif 23 <= age <= 25:
        return "23-25"
    elif 29 <= age <= 30:
        return "29-30"
    else:
        return None

def cliffs_delta(x, y):
    x, y = np.array(x), np.array(y)
    return (np.sum(x[:, None] > y) - np.sum(x[:, None] < y)) / (len(x) * len(y))

def normalize_posterior(p):
    if np.nansum(p) == 0:
        return np.full_like(p, np.nan)
    return p / np.nansum(p)

def in_object_zone(pos, object_zones):
    mask = np.zeros_like(pos, dtype=bool)
    for zmin, zmax in object_zones:
        mask |= (pos >= zmin) & (pos <= zmax)
    return mask

# -----------------------------
# Metrics
# -----------------------------
def posterior_metrics(posterior, real_pos, bin_centers):
    """Calcule MASS, OFFDIAG, WIDTH, ERROR, CORR."""
    decoded_pos = bin_centers[np.nanargmax(posterior, axis=1)]
    error = np.abs(decoded_pos - real_pos)

    diag_idx = np.digitize(real_pos, bin_centers) - 1
    diag_idx = np.clip(diag_idx, 0, posterior.shape[1]-1)
    diag_mass = posterior[np.arange(len(posterior)), diag_idx]
    total_mass = np.nansum(posterior, axis=1)
    mass = diag_mass / total_mass
    offdiag = 1 - mass

    mean_pos = np.nansum(posterior * bin_centers, axis=1)
    width = np.sqrt(np.nansum(posterior * (bin_centers - mean_pos[:, None])**2, axis=1))

    corr = np.corrcoef(real_pos, decoded_pos)[0,1] if len(real_pos) > 1 else np.nan

    return {"MASS": mass, "OFFDIAG": offdiag, "WIDTH": width, "ERROR": error, "CORR": corr}

# -----------------------------
# Boxplots et stats
# -----------------------------
def boxplot_metric_by_age(results, metric, direction, age_groups_order):
    data = [results[metric][direction][g] for g in age_groups_order]

    # ---------- PRINT MOYENNES ----------
    print(f"\n{metric.upper()} | {direction} — Moyennes par groupe d'âge")
    for g, d in zip(age_groups_order, data):
        if len(d) > 0:
            print(f"  {g} : mean = {np.nanmean(d):.4f} ± {np.nanstd(d):.4f} (n={len(d)})")
        else:
            print(f"  {g} : aucune donnée")

    data_nonempty = [d for d in data if len(d) > 0]
    if len(data_nonempty) == 0:
        print(f"Aucune donnée pour {metric} | {direction}")
        return

    plt.figure(figsize=(8,5))
    plt.boxplot(data, labels=age_groups_order, showfliers=False)
    plt.ylabel(metric.upper())
    plt.xlabel("Groupe d'âge")
    plt.title(f"{metric.upper()} – {direction}")
    plt.grid(True, axis='y')

    # ---------- KRUSKAL-WALLIS ----------
    stat, p_global = kruskal(*data_nonempty)
    print(f"{metric.upper()} | {direction} - Kruskal-Wallis p = {p_global:.4e}")

    # ---------- COMPARAISONS SUCCESSIVES ----------
    y_max = np.nanmax([np.nanmax(d) for d in data_nonempty])
    step = 0.07 * y_max

    for i in range(len(data)-1):
        g1, g2 = data[i], data[i+1]
        if len(g1) > 0 and len(g2) > 0:
            _, p = mannwhitneyu(g1, g2)
            if p < 0.001: label = "***"
            elif p < 0.01: label = "**"
            elif p < 0.05: label = "*"
            else: label = "n.s."
            plt.text(i+1.5, y_max + step*i, label,
                     ha='center', color='red', fontsize=12)

    plt.tight_layout()
    plt.show()

def plot_box_zone_comparison(object_data, non_object_data, ylabel, title, paired=False):
    object_data = [v for v in object_data if not np.isnan(v)]
    non_object_data = [v for v in non_object_data if not np.isnan(v)]

    if len(object_data) == 0 or len(non_object_data) == 0:
        print(f"Aucune donnée pour {title}")
        return

    # ---------- PRINT MOYENNES ----------
    print(f"{title}")
    print(f"  Zone objet : mean = {np.mean(object_data):.4f} ± {np.std(object_data):.4f} (n={len(object_data)})")
    print(f"  Hors zone  : mean = {np.mean(non_object_data):.4f} ± {np.std(non_object_data):.4f} (n={len(non_object_data)})")

    plt.figure(figsize=(6,5))
    plt.boxplot([object_data, non_object_data],
                labels=["Zone objet", "Hors zone"])
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, axis='y')

    y_max = max(max(object_data), max(non_object_data)) * 1.1

    if paired and len(object_data) == len(non_object_data):
        stat, pval = wilcoxon(object_data, non_object_data)
        test_name = "Wilcoxon"
    else:
        stat, pval = mannwhitneyu(object_data, non_object_data)
        test_name = "Mann-Whitney"

    delta = cliffs_delta(object_data, non_object_data)

    plt.text(1.5, y_max, f"{test_name} p = {pval:.3e}",
             ha='center', fontsize=12, color='red')
    plt.text(1.5, y_max*1.05, f"Cliff's δ = {delta:.2f}",
             ha='center', fontsize=12, color='blue')

    plt.show()

    print(f"  → p = {pval:.4e}, Cliff's δ = {delta:.2f}")

# -----------------------------
# Heatmap metrics
# -----------------------------
def diagonal_width(heatmap):
    widths = []
    for i in range(heatmap.shape[0]):
        p = heatmap[i,:]
        if np.all(np.isnan(p)): continue
        p = normalize_posterior(p)
        pos = np.arange(len(p))
        mu = np.nansum(pos * p)
        sigma = np.sqrt(np.nansum((pos - mu)**2 * p))
        widths.append(sigma)
    return np.nanmean(widths)

def diagonal_mass(heatmap, band=1):
    vals = [heatmap[i,j] for i in range(heatmap.shape[0]) for j in range(heatmap.shape[1]) if abs(i-j)<=band]
    return np.nanmean(vals)

def off_diagonal_energy(heatmap, band=2):
    vals = [heatmap[i,j] for i in range(heatmap.shape[0]) for j in range(heatmap.shape[1]) if abs(i-j)>band]
    return np.nanmean(vals)

def diagonal_error(heatmap):
    errors = []
    for i in range(heatmap.shape[0]):
        p = heatmap[i,:]
        if np.all(np.isnan(p)): continue
        p = normalize_posterior(p)
        decoded = np.nansum(np.arange(len(p))*p)
        errors.append(abs(decoded - i))
    return np.nanmean(errors)

def diagonal_correlation(heatmap):
    real, decoded = [], []
    for i in range(heatmap.shape[0]):
        p = heatmap[i,:]
        if np.all(np.isnan(p)): continue
        p = normalize_posterior(p)
        real.append(i)
        decoded.append(np.nansum(np.arange(len(p))*p))
    if len(real) < 2: return np.nan
    return np.corrcoef(real, decoded)[0,1]

# -----------------------------
# Chargement des données
# -----------------------------
data = np.load(r'C:\Users\suire\Documents\Python Scripts\saved_sessions\pooled_data.npz', allow_pickle=True)
session_names = data["session_names"]
decoded_positions = data["decoded_positions"]
posteriors = data["posteriors"]
posterior_tvecs = data["posterior_tvecs"]
trajs_all_dirs = data["trajs_all_dirs"]
ages = data["age_list"]
object_zones = [(20,55),(120,135)]
age_groups_order = ["15-17","18-20","21-22","23-25","29-30"]

# -----------------------------
# Stockage metrics par session, par direction, par zone
# -----------------------------
results = {}
xmin, xmax = 0, 145

for i in range(len(session_names)):
    age_group = group_age(ages[i])
    if age_group is None:
        continue
    posterior = posteriors[i]
    tvec = posterior_tvecs[i]
    trajs = trajs_all_dirs[i]
    n_bins = posterior.shape[1]
    bin_centers = np.linspace(xmin, xmax, n_bins)
    bin_edges = np.linspace(xmin, xmax, n_bins+1)
    for direction in trajs.keys():
        for idx_traj, traj in enumerate(trajs[direction]):

            direction_name = "Aller" if idx_traj % 2 == 0 else "Retour"


            real_pos = np.interp(tvec, traj.tvec, traj.data, left=np.nan, right=np.nan)
            # après interpolation
            mask = ~np.isnan(real_pos)
            real_pos = real_pos[mask]
            post = posterior[mask, :]

            # 🔒 nouveau masque
            valid_post = ~np.all(np.isnan(post), axis=1)
            real_pos = real_pos[valid_post]
            post = post[valid_post, :]

            if post.shape[0] == 0:
                continue

            # Heatmap
            heatmap = np.full((n_bins, n_bins), np.nan)
            for j in range(len(bin_edges)-1):
                idx_bin = (real_pos >= bin_edges[j]) & (real_pos < bin_edges[j+1])
                if np.any(idx_bin):
                    heatmap[j,:] = np.nanmean(post[idx_bin,:], axis=0)

            # Metrics
            metrics = posterior_metrics(post, real_pos, bin_centers)
            obj_mask = in_object_zone(real_pos, object_zones)
            nonobj_mask = ~obj_mask

            for zone, mask_zone in zip(["OBJ","NONOBJ"], [obj_mask, nonobj_mask]):
                key = (age_group, direction_name, zone)
                if key not in results:
                    results[key] = {k: [] for k in metrics.keys()}
                for k in ["MASS","OFFDIAG","WIDTH","ERROR"]:
                    if np.sum(mask_zone) > 0:
                        results[key][k].extend(metrics[k][mask_zone])
                results[key]["CORR"].append(metrics["CORR"])

# -----------------------------
# Boxplots zones objets
# -----------------------------
for metric in ["MASS","OFFDIAG","WIDTH","ERROR"]:
    for direction_name in ["Aller","Retour"]:
        print(f"\n=== {metric} | {direction_name} ===")
        for age_group in age_groups_order:
            k_obj = (age_group,direction_name,"OBJ")
            k_non = (age_group,direction_name,"NONOBJ")
            obj_data = results.get(k_obj,{}).get(metric,[])
            non_data = results.get(k_non,{}).get(metric,[])
            plot_box_zone_comparison(obj_data, non_data, ylabel=metric, title=f"{metric} | {direction_name} | {age_group}")

# -----------------------------
# Boxplots par âge
# -----------------------------
metrics_for_age = ["MASS","OFFDIAG","WIDTH","ERROR","CORR"]
results_by_age = {m: {"Aller": {g: [] for g in age_groups_order},
                      "Retour": {g: [] for g in age_groups_order}} for m in metrics_for_age}

for key, metric_dict in results.items():
    age_group, direction, zone = key
    for m in metrics_for_age:
        results_by_age[m][direction][age_group].extend(metric_dict[m])

for m in metrics_for_age:
    for direction_name in ["Aller","Retour"]:
        boxplot_metric_by_age(results_by_age, metric=m, direction=direction_name, age_groups_order=age_groups_order)
plt.show()