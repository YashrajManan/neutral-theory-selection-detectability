# Neutral Theory vs. Selection: Detectability Edge Cases

### Does statistical power to detect selection increase monotonically with selection strength — or is there a window where it doesn't?

## 1. Pitch

"We didn't detect selection" is often treated as evidence selection is absent, but a real, standing
critique of Kimura's Neutral Theory is that this conflates "absent" with "undetectable by this test at this
sample size." This project builds a real population-genetic forward simulation and a real neutrality test
from scratch, then runs an actual statistical power analysis to find out exactly where that boundary sits —
and finds it isn't where the textbook framing assumes.

## 2. Research question

At a given population size, sample size, and selection strength, what is the real statistical power of
Tajima's D — a standard, widely-used neutrality test — to correctly flag a selected site as non-neutral,
measured directly via Monte Carlo simulation rather than assumed from theory?

## 3. Biological background

Kimura's Neutral Theory (1968) proposes most molecular variation is selectively neutral. A persistent
methodological critique is that many standard neutrality tests have limited power against weak selection,
especially at small sample sizes — meaning a result "consistent with neutrality" may reflect low test power
rather than genuine absence of selection. This project builds the real machinery to quantify that boundary:
an infinite-sites Wright-Fisher forward simulation (new mutations enter a finite population at a constant
rate, each evolving under selection + drift until fixed or lost) and Tajima's D, derived directly from
Tajima's 1989 formula and cross-checked against an independent library implementation.

## 4. Data fetch

None — and deliberately so: the question here is a statistical test's own detection power, which has to be
established against a known, controlled ground truth (an exact selection coefficient and population size)
that no real dataset could provide, since the true presence or strength of selection in a real population
is exactly what's unknown and being tested for. The simulated Wright-Fisher populations, with their fully
known generating parameters, are the real object of study.

## 5. Data describe

Each simulation tracks a haploid population of size `N=150`. New mutations enter at a constant per-generation
rate (combined input `lam_per_genome x R` for `R` independent replicate "genomes," since unlinked sites never
interact and can be pooled into one simulation run and split afterward). Each currently-segregating site's
derived-allele count is updated every generation by selection-weighted reweighting followed by binomial
resampling (drift); a site that reaches 0 or `N` copies leaves the segregating pool for good. A finite sample
of `n` individuals is then drawn from the population to compute Tajima's D, mirroring how real sequencing
studies sample a subset of a much larger natural population.

## 6. Methods / Workflow — what we did

1. Built the infinite-sites Wright-Fisher forward simulation; confirmed under neutrality the model reaches
   a stable, non-trivial equilibrium (not collapsing to 0% or 100% derived-allele frequency), and that
   positive/negative selection visibly shifts that equilibrium in the expected direction.
2. Derived Tajima's D directly from Tajima's 1989 formula constants and cross-checked it against
   `scikit-allel`'s independent implementation on a hand-built example — confirmed an exact match.
3. Built an efficient "many independent replicate genomes from one simulation run" trick: since sites are
   unlinked, one simulation with a combined mutation rate can be split into many independent replicate
   samples afterward, avoiding hundreds of separate from-scratch burn-ins.
4. Built the neutral null distribution of D empirically and defined a two-sided 95% rejection region from
   its own percentiles — the correct way to calibrate a Monte Carlo power analysis.
5. Swept selection strength (`s`, equivalently `Ns`) and sample size (`n`) and computed real empirical
   power at every combination, without assuming the shape of the relationship in advance.

## 7. Results

Tajima's D matched `scikit-allel`'s independent implementation exactly (`-0.279223` both ways). The neutral
null calibrated correctly (mean D close to 0, as required by the test's own construction). Power as a
function of selection strength is **not monotonic**: it rises from the ~5% baseline, peaks around
`Ns≈1.5-3` (reaching `0.24-0.32` at the largest sample size tested), then falls back down as selection gets
stronger — in several cases below the nominal 5% false-positive rate by `Ns=60`. Power increases cleanly
with sample size at a fixed, moderately detectable selection strength (`s=0.05`): from `0.003-0.06` at
`n=10` up to `0.22-0.26` at `n=80`.

## 8. Biology interpretation of results

Tajima's D compares two estimators of the same mutation parameter that respond differently to *where* in
the frequency spectrum variation sits; under neutrality their expected difference is exactly zero by
construction, which is why the null calibration check matters so much. The non-monotonic power curve has a
real mechanistic explanation: a new, very strongly favored mutation sweeps to fixation (or is lost) so
quickly that very few such mutations are ever caught still segregating at an informative intermediate
frequency at any single census snapshot — the ones still segregating under very strong selection are
disproportionately very young and statistically similar to young neutral mutations. Moderate selection
gives mutations enough of an edge to reach observable frequencies while still lingering there long enough
to leave a detectable signature. This is a real, specific limitation of steady-state, snapshot-based
neutrality testing — not a claim that strong selection is generally undetectable by every method (sweep
scans built around linked neutral diversity are a different, complementary tool suited to exactly this
case).

## Learning through project

A test's own internal calibration checkpoint is worth more than eyeballing downstream numbers for
plausibility: before trusting any power estimate, this project's neutral null was first confirmed to
calibrate correctly on its own terms (mean Tajima's D close to 0, exactly as required by the test's
mathematical construction) — a check that catches a mis-specified simulation long before any selection
sweep would be worth interpreting. More broadly: "detectability increases with effect size" is an
intuitive but unexamined assumption for any detection method built around sampling a fast-moving process
at a single point in time — a very fast signal can be *harder* to catch than a moderate one, purely
because of how observation timing interacts with the underlying process's own speed.

## Reproducible code

`neutral_theory_selection_detectability.ipynb` — runs top to bottom locally (CPU only, no GPU or dataset
download needed; `pip install scikit-allel` if not already installed).

## Key tools

NumPy (forward simulation), a hand-derived Tajima's D cross-checked against `scikit-allel`, pandas (power
sweep tables), Matplotlib (hero figure).

## Skills demonstrated

Forward population-genetic simulation design, deriving a classical statistic from its own formula and
independently validating it, efficient Monte Carlo replicate generation via unlinked-site pooling, a real
statistical power analysis with empirically-calibrated null distributions, and catching a silent numerical
bug via a test's own mathematical self-consistency check rather than downstream plausibility judgment.

## Stretch goals

Testing whether the same non-monotonic detectability window appears for other standard neutrality tests
(Fu & Li's D, Fay & Wu's H); extending to a demographic model with population size changes; comparing
against a real published site-frequency spectrum to ask whether that study's sample size had power to
detect selection in a biologically plausible `Ns` range.
