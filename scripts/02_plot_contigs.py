"""Overview plots of the D. lutheri genome: how contiguous is it, and do any contigs look unusual?"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("results/contigs.tsv", sep="\t", header=None, usecols=[0, 1, 2])
df.columns = ["contig", "length", "gc"]
df = df.sort_values("length", ascending=False).reset_index(drop=True)

total = df["length"].sum()
cumulative = df["length"].cumsum()
n50 = df.loc[cumulative >= total / 2, "length"].iloc[0]
print(f"Contigs: {len(df)}   Total: {total/1e6:.2f} Mb   N50: {n50/1e3:.0f} kb")
print(f"GC: mean {df['gc'].mean():.2f}%, range {df['gc'].min():.2f}-{df['gc'].max():.2f}%")

mean, sd = df["gc"].mean(), df["gc"].std()
odd = df[(df["gc"] - mean).abs() > 3 * sd]
print(f"Unusual-GC contigs: {len(odd)}")
print(odd.to_string(index=False))

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))

left.step(range(1, len(df) + 1), cumulative / 1e6, where="post", color="#2E5395")
left.axhline(total / 2 / 1e6, ls="--", color="grey")
left.text(len(df) * 0.45, total / 2 / 1e6 + 1, f"half the genome (N50 = {n50/1e3:.0f} kb)", color="grey")
left.set_xlabel("Contigs, longest first")
left.set_ylabel("Cumulative length (Mb)")
left.set_title("How fast the genome adds up")

right.scatter(df["gc"], df["length"] / 1e6, s=18, alpha=0.7, color="#2E5395")
right.scatter(odd["gc"], odd["length"] / 1e6, s=40, color="#C00000", label="unusual GC")
right.set_xlabel("GC content (%)")
right.set_ylabel("Contig length (Mb)")
right.set_title("Does every contig look alike?")
if len(odd):
    right.legend()

fig.suptitle("Diacronema lutheri reference genome (GCA_019448385.1)")
fig.tight_layout()
fig.savefig("results/contig_overview.png", dpi=150)
print("Saved results/contig_overview.png")
