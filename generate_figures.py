"""
Generate 3 publication-quality figures (PDF + 300-DPI PNG) for the Giuga's Conjecture amsart paper.
"""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\Admin\.gemini\antigravity\scratch\giuga_paper\figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
})

# =====================================================================
# Figure 1: Dual-LP Certified Lower Bounds vs Cutoff X (k and log10(n))
# =====================================================================
X_vals = np.array([50000, 100000, 200000, 300000, 321000]) / 1000.0  # in thousands
k_vals = np.array([3914, 5714, 7538, 8122, 8138])
dig_X  = np.array([100000, 200000, 300000, 321000]) / 1000.0
dig_v  = np.array([27045.46, 37237.38, 40639.82, 40730.99])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6))

ax1.plot(X_vals, k_vals, marker="o", lw=2.0, color="#1f4e79", label="Dual-LP Certified $k$")
ax1.axhline(4771, color="#b22222", ls="--", lw=1.3, label="Borwein–Maitland–Skerritt (2013): 4,771")
ax1.axhline(3459, color="#7f7f7f", ls=":", lw=1.2, label="Borwein et al. (1996): 3,459")
ax1.annotate("8,138", xy=(321, 8138), xytext=(255, 8450),
             arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0),
             fontweight="bold", color="#1f4e79")
ax1.set_title("(a) Certified Prime Factors $k$ vs. Cutoff $X$")
ax1.set_xlabel("DAG Cutoff $X$ (thousands)")
ax1.set_ylabel("Certified Lower Bound on $k$")
ax1.set_ylim(2800, 9100)
ax1.grid(alpha=0.25)
ax1.legend(loc="lower right", frameon=True)

ax2.plot(dig_X, dig_v, marker="s", lw=2.0, color="#1b6b3a", label="Dual-LP Certified $\\log_{10}(n)$")
ax2.axhline(36067, color="#d97706", ls="-.", lw=1.3, label="Sorini (2001) heuristic: 36,067")
ax2.axhline(19907, color="#b22222", ls="--", lw=1.3, label="Borwein–Maitland–Skerritt (2013): 19,907")
ax2.annotate("40,731 digits", xy=(321, 40731), xytext=(205, 42200),
             arrowprops=dict(arrowstyle="->", color="#1b6b3a", lw=1.0),
             fontweight="bold", color="#1b6b3a")
ax2.set_title("(b) Certified Decimal Digits $\\log_{10}(n)$ vs. Cutoff $X$")
ax2.set_xlabel("DAG Cutoff $X$ (thousands)")
ax2.set_ylabel("Certified Lower Bound on $\\log_{10}(n)$")
ax2.set_ylim(16000, 45500)
ax2.grid(alpha=0.25)
ax2.legend(loc="lower right", frameon=True)

fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_lp_convergence.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig1_lp_convergence.png"), bbox_inches="tight")
plt.close(fig)

# =====================================================================
# Figure 2: Self-Throttling No-Divisor-Chain Sieve vs Unconstrained Sum
# =====================================================================
def primes_up_to(N):
    s = [True]*(N+1); s[0]=s[1]=False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = [False]*len(s[i*i::i])
    return [i for i, v in enumerate(s) if v]

ps_all = [p for p in primes_up_to(350000) if p >= 5]

# Unconstrained partial sums
u_k, u_s = [], []
acc = 0.0
for idx, p in enumerate(ps_all[:40], 1):
    acc += 1.0 / p
    u_k.append(idx); u_s.append(acc)

def greedy_curve(start_p):
    adm, ks, ss = [], [], []
    acc = 0.0
    for p in ps_all:
        if p < start_p: continue
        if all((p - 1) % q != 0 for q in adm):
            adm.append(p)
            acc += 1.0 / p
            if len(adm) <= 30 or len(adm) % 40 == 0 or acc >= 1.0:
                ks.append(len(adm)); ss.append(acc)
            if acc >= 1.0:
                break
    return ks, ss

k5, s5 = greedy_curve(5)
k7, s7 = greedy_curve(7)

fig, ax = plt.subplots(figsize=(7.6, 3.8))
ax.plot(u_k, u_s, lw=2.0, color="#7f7f7f", ls="--", label="Unconstrained primes $\\geq 5$ (crosses $1.0$ at $k=27$)")
ax.plot(k5, s5, lw=2.2, color="#1f4e79", label="No-Divisor-Chain, $p_1 = 5$ (crosses $1.0$ at $k=8{,}135$)")
ax.plot(k7, s7, lw=2.0, color="#b22222", ls="-.", label="No-Divisor-Chain, $p_1 = 7$ (still $<1.0$ at $k=9{,}500+$)")
ax.axhline(1.0, color="black", lw=1.0, ls=":")
ax.set_xscale("log")
ax.set_xlabel("Number of Selected Primes $k$ (log scale)")
ax.set_ylabel("Reciprocal Sum $\\sum_{i=1}^k 1/p_i$")
ax.set_title("Self-Throttling Effect of the Carmichael No-Divisor-Chain Constraint ($q$ does not divide $p_i - 1$)")
ax.set_ylim(0.15, 1.12)
ax.grid(alpha=0.25, which="both")
ax.legend(loc="lower right", frameon=True)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_dag_sieve.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig2_dag_sieve.png"), bbox_inches="tight")
plt.close(fig)

# =====================================================================
# Figure 3: Historical Progression of Giuga's Conjecture Digit Bounds
# =====================================================================
labels = [
    "Giuga\n(1950)",
    "Bedocchi\n(1985)",
    "Borwein et al.\n(1996)",
    "Borwein–Maitland–\nSkerritt (2013)",
    "Sorini (2001)\n[heuristic]",
    "Barmase & Neon\n(2026, certified)"
]
digits = [1000, 1700, 13887, 19908, 36067, 40731]
colors = ["#94a3b8", "#94a3b8", "#64748b", "#3b82f6", "#f59e0b", "#15803d"]

fig, ax = plt.subplots(figsize=(8.0, 3.8))
bars = ax.bar(labels, digits, color=colors, width=0.58, edgecolor="#1e293b", linewidth=0.7)
for bar, val in zip(bars, digits):
    ax.text(bar.get_x() + bar.get_width()/2, val + 900, f"{val:,}",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_ylabel("Minimum Decimal Digits of $n$")
ax.set_title("Historical Progression of Lower Bounds for Counterexamples to Giuga's Conjecture")
ax.set_ylim(0, 46500)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_historical_bounds.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig3_historical_bounds.png"), bbox_inches="tight")
plt.close(fig)
print("ALL 3 FIGURES GENERATED OK")
