# Compartments: the body-water lens on the energy of S at N = 10^6

> Exploration at finite X, not a proof. Every number below is read from `cache/compartments.pkl` (grids M = 2^17 … 2^24,
> X = 10^4 … 2·10^6); the physiology numbers come from the two sources listed at the end. Body-water physiology is used as
> a teacher for bookkeeping, not as a mechanism: nothing in S moves, mixes or equilibrates.

**What physiology does.** Total body water (TBW) is about 60% of an adult man's weight (42 L of 70 kg). Two thirds of it
is intracellular (28 L) and one third extracellular (14 L: 10.5 L interstitial and 3.5 L plasma) (StatPearls, *Physiology,
Water Balance*). Compartments are measured by **indicator dilution**: give a known amount I of an indicator that spreads
only through the compartment of interest, let it equilibrate, read its plasma concentration C, and set V = I/C. Heavy
water (D₂O, H₂¹⁸O) or antipyrine measures TBW; bromide, sucrose or inulin measures extracellular water; Evans blue
measures plasma. The rest follows by subtraction: intracellular = TBW − extracellular, interstitial = extracellular −
plasma, solids = body weight − TBW (Abdominal Key, ch. 2). The composition changes with growth. TBW falls from 94% of body
weight in the first trimester to about 86% at 24 weeks, 78% at term and about 70% at one year, as solids accrete.
Extracellular water falls (60% → 44% → ≈ 28%) and intracellular water rises (≈ 26% → 34% → 42%). It also depends on growth
conditions: infants small for gestational age have more water and less fat; infants large for gestational age have more
fat and less water.

**What the circle has.** One conserved total, by Parseval: (1/M)Σ_k |S(k/M)|² = ΣΛ² exactly. The windows
|N(α − a/q)| ≤ ½ around the reduced fractions a/q of a set of denominators play the indicators. A window reads the energy it
covers: a known share (its capture) of the Dirichlet spike centred in it, plus the sea it covers. Three compartments:

- **atoms**: the spikes at q = 1, 2 (α = 0 and ½);
- **forest**: the spikes with 3 ≤ q ≤ Q = ⌊q*(X)⌋, where q* is the denominator at which spike and sea are equal
  (fractals.md, M1);
- **sea**: the rest, i.e. the structureless part plus the spikes with q > Q that it drowns.

The no-fit spike model is S(a/q + β) ≈ μ(q)(ψ/φ(q))·v(β)/X, with v the Dirichlet kernel of energy X. So the φ(q) spikes of
denominator q carry ψ²μ²(q)/(Xφ(q)) together, and the shares are atoms = 2/L_N and forest = (Σ_{3≤q≤Q} μ²/φ)/L_N, with
L_N = XΣΛ²/ψ². L_N is the total energy divided by the energy of one reference spike (the one at 0, ψ²/X), a dilution
volume in the literal sense of a ratio. At X = 10^6, L_N = 12.81: the energy of S is spread as if over
12.8 copies of its main spike.

## 0. Certificates

| check | result | detail |
|:---|:---|:---|
| Parseval, primes, X = 10000 | pass |  |
| centre sea (reduced classes) = divisor sum of residue variances, X = 10000 | pass | max \|ratio - 1\| = 6.8e-14; prime powers p^k, p \| q, move it by up to 0.143 |
| Parseval, primes, X = 30000 | pass |  |
| centre sea (reduced classes) = divisor sum of residue variances, X = 30000 | pass | max \|ratio - 1\| = 2.9e-14; prime powers p^k, p \| q, move it by up to 0.259 |
| Parseval, primes, X = 100000 | pass |  |
| centre sea (reduced classes) = divisor sum of residue variances, X = 100000 | pass | max \|ratio - 1\| = 1.8e-13; prime powers p^k, p \| q, move it by up to 0.033 |
| Parseval, primes, X = 300000 | pass |  |
| centre sea (reduced classes) = divisor sum of residue variances, X = 300000 | pass | max \|ratio - 1\| = 7.0e-13; prime powers p^k, p \| q, move it by up to 0.148 |
| Parseval, primes, X = 1000000 | pass |  |
| centre sea (reduced classes) = divisor sum of residue variances, X = 1000000 | pass | max \|ratio - 1\| = 3.0e-12; prime powers p^k, p \| q, move it by up to 0.054 |
| window energy = capture x model, q <= 30 at X = 10^6 | pass | max \|ratio - 1\| = 0.0004 |
| Parseval, primes, X = 2000000 | pass |  |
| centre sea (reduced classes) = divisor sum of residue variances, X = 2000000 | pass | max \|ratio - 1\| = 1.3e-12; prime powers p^k, p \| q, move it by up to 0.072 |
| Parseval, Cramér | pass |  |
| Parseval, odd Cramér | pass |  |
| Parseval, primes | pass |  |
| Parseval, primes ≡ 1 (mod 4), weight 2Λ | pass |  |

The divisor identity is exact once the variances are taken over n coprime to q. With the prime powers p^k (p | q) that
S(a/q) actually contains, the centre level moves by up to 26% for a few q at small X; the local-sea reading
uses the actual centre level.

Each window's energy against its capture times the model spike, at X = 10^6 (small q, where the sea is negligible):

| q | window energy / (capture × model) | window capture, min–max |
|---:|---:|---:|
| 1 | 1.0002 | 0.8030–0.8030 |
| 2 | 1.0002 | 0.8030–0.8030 |
| 3 | 1.0000 | 0.7554–0.7554 |
| 5 | 1.0003 | 0.7539–0.7558 |
| 6 | 1.0000 | 0.7554–0.7554 |
| 7 | 1.0000 | 0.7549–0.8025 |
| 10 | 1.0002 | 0.7539–0.7558 |
| 11 | 1.0001 | 0.7548–0.8028 |
| 13 | 0.9998 | 0.7543–0.8028 |
| 14 | 0.9999 | 0.7549–0.8025 |
| 15 | 1.0004 | 0.7547–0.8029 |
| 30 | 1.0004 | 0.7547–0.8029 |

## B1. Compartments at X = 10^6

| indicator | what it reads | share of ΣΛ² |
|:---|:---|---:|
| whole circle (Parseval) | the conserved total ΣΛ² | 100% |
| windows at q = 1, 2 | raw reading | 12.54% |
| windows at 3 ≤ q ≤ Q = 512 | raw reading | 34.94% |

The readings become compartments by a three-window solution. Each window reads its capture of its own spikes plus the
sea it covers; the rest of the circle holds the spikes' leakage and the remaining sea:

| compartment | one average capture, uniform sea | exact capture, uniform sea | exact capture, sea from Σ_{d\|q} of X·FG(d) (no fit) | exact capture, sea measured at the centres | no-fit model |
|:---|---:|---:|---:|---:|---:|
| atoms (q = 1, 2) | 16.20% | 15.61% | 15.61% | 15.61% | 15.61% |
| forest (3 ≤ q ≤ Q) | 42.11% | 42.35% | 43.66% | 43.64% | 43.53% |
| sea (the rest) | 41.68% | 42.04% | 40.73% | 40.75% | 40.86% |

Two corrections matter, and each resembles one that tracer practice needs:

1. **Incomplete recovery.** A window catches only part of its spike, and how much depends on where the fraction falls
   between grid points. At X = 10^6 a window holds 8 or 9 grid points and catches 0.755 or 0.803 of its spike, not the
   average 0.7737. Using one average capture misreads the atoms by +3.8%. The error changes sign
   with the grid offset, so it is a recovery correction rather than a tracer leaking one way.
2. **The tracer outside its space, and uneven mixing.** A window also reads sea that is not its spike's: that is the
   counterpart of an indicator that reaches beyond its compartment (in standard tracer practice, not from the two
   sources: heavy water exchanges with non-aqueous hydrogen, and bromide enters cells a little). Plain dilution assumes
   that sea is evenly mixed. It is not (next table). With a uniform sea the forest reads -2.7%
   against the model. With the sea level measured at the spike centres it reads +0.3%; with the
   no-fit divisor form below it reads +0.3%.

**The sea is thinner near the forest.** Sea level per unit α in the windows of each band of denominators, at X = 10^6:

| windows of | # q | measured in the windows | measured at the centres | divisor sum of X·FG(d) (no fit) | uniform sea |
|:---|---:|:---|:---|:---|:---|
| q ∈ [31, 97] | 42 | 1.72 X | 1.69 X | 1.52 X | 5.38 X |
| q ∈ [101, 255] | 96 | 2.82 X | 2.66 X | 2.57 X | 5.38 X |
| q ∈ [257, 511] | 157 | 3.56 X | 3.36 X | 3.34 X | 5.38 X |

Near a reduced fraction a/q (q squarefree) the structureless part of S has mean square

    (1/φ(q)) Σ_{d|q} μ(q/d) · d · V(X; d),

where V(X; d) is the variance of the primes over the reduced residue classes mod d (V(X; 1) = V(X; 2) = 0). This follows
from Σ_{a reduced} e(a(r − s)/q) = c_q(r − s) = Σ_{d|q, d|r−s} d μ(q/d); the certificate above checks it to rounding. For prime
q it is V(X; q)·q/(q − 1), about X·FG(q) with FG(q) = ln q − γ − ln 2π − ln q/(q − 1) (Friedlander–Goldston form of the
Barban–Davenport–Halberstam / Montgomery–Hooley variance). For composite q the divisor sum is larger. Replacing each
V(X; d) by X·FG(d) gives the no-fit column, which falls 6–11% below the window values.
Every band sits well below the circle's average sea level, and the level climbs back toward it as q grows. (Bands with
q ≤ 30 are left out: there the sea is too small a fraction of the window to measure.)

## B2. Development: the composition of S as X grows

| X | L_N | q* (Q = ⌊q*⌋) | atoms: measured / model | forest: measured / model | sea: measured / model | forest with a uniform sea |
|---:|---:|---:|---:|---:|---:|---:|
| 1e4 | 8.15 | 74.8 (74) | 0.2447 / 0.2453 | 0.4440 / 0.4484 | 0.3114 / 0.3063 | 0.4184 |
| 3e4 | 9.29 | 115.3 (115) | 0.2154 / 0.2153 | 0.4412 / 0.4422 | 0.3434 / 0.3424 | 0.4197 |
| 1e5 | 10.49 | 189.9 (189) | 0.1906 / 0.1907 | 0.4355 / 0.4352 | 0.3739 / 0.3741 | 0.4177 |
| 3e5 | 11.61 | 303.1 (303) | 0.1723 / 0.1722 | 0.4359 / 0.4348 | 0.3919 / 0.3930 | 0.4208 |
| 1e6 | 12.81 | 512.9 (512) | 0.1561 / 0.1561 | 0.4364 / 0.4353 | 0.4075 / 0.4086 | 0.4235 |
| 2e6 | 13.50 | 697.5 (697) | 0.1481 / 0.1481 | 0.4367 / 0.4358 | 0.4152 / 0.4161 | 0.4250 |

With the centre-measured sea, the atoms and the forest match the no-fit model to within 0.3% at X ≥ 3·10^4
(1.0% at 10^4).

For comparison, the body (percent of body weight; fetal and infant values from the Abdominal Key chapter, adult from
StatPearls; "≈" marks a value obtained by subtraction):

| stage | total water | extracellular | intracellular | solids |
|:---|---:|:---|:---|---:|
| first trimester | 94% | — | — | ≈ 6% |
| 24 weeks | 86% | 60% | ≈ 26% | ≈ 14% |
| term | 78% | 44% | 34% | 22% |
| 1 year | 70% | ≈ 28% | 42% | 30% |
| adult man | 60% | 20% | 40% | ≈ 40% |

Both show one conserved total redistributing among compartments as a growth parameter increases. The **atoms** fall like
2/L_N ≈ 2/(ln X − 1), from 24.5% at X = 10^4 to 14.8% at
2·10^6: their energy 2ψ²/X grows only like X, while the total grows like X ln X. The **forest** stays near 0.44 and the **sea**
grows.

The leading-order law below is a model, not a theorem: it uses the Friedlander–Goldston asymptotic, the spike model up to
q* and q* from q = √(X/FG(q)). In that model the shares tend to **0, ½, ½**, because q* grows like √X up to logarithms and the
forest holds (ln q* + 1.33 − 2)/(ln X − 1). The atoms approach 0 like 2/ln X, but the forest and the sea approach ½ only like
ln ln X/(2 ln X), and the sea passes ½ near X ≈ 10^38 and returns to it from above:

| X | atoms | forest | sea |
|:---|---:|---:|---:|
| 10^6 | 0.156 | 0.435 | 0.409 |
| 10^12 | 0.075 | 0.450 | 0.475 |
| 10^20 | 0.044 | 0.464 | 0.492 |
| 10^40 | 0.022 | 0.478 | 0.500 |
| 10^100 | 0.009 | 0.489 | 0.502 |
| 10^1000 | 0.001 | 0.498 | 0.501 |

**Where the analogy stops.** The coherent share of S (atoms plus forest) *falls* with X, from 68.9% to 58.5%
(model limit ½), while a fetus *gains* solids. Which compartment of S "corresponds" to extracellular or intracellular water
is a choice with no physical basis. Only the bookkeeping is literal; the limit belongs to the leading-order model. X is a
parameter, not time, and there is no counterpart to birth or to the postnatal contraction of extracellular water.

## B3. "Growth aberration": the arithmetic of the set sets its composition (X = 10^6)

Four sets with the same total weight (≈ ψ) but different congruence structure. Each set's no-fit model spreads its weight
evenly over the residue classes it allows mod lcm(D, q) and sums them with e(ab/q). **Atoms** are the set's exact copies of
its main spike (model height equal to the height at q = 1); the forest windows sit where the set's model has spikes, up to
the primes' Q = 512.

| set (X = 10^6, same total weight ≈ ψ) | what it keeps | atom denominators (full-height spikes) | atoms: measured / model | forest: measured / model | sea | # forest denominators with model spikes |
|:---|:---|---:|---:|:---|---:|---:|
| Cramér | random n ≥ 3 with probability 1/ln n, weight ln n: no congruence structure | 1 | 0.0780 / 0.0779 | no forest windows (model has none); control in the primes' windows: +0.0001 | 0.9220 | 0 |
| odd Cramér | random odd n with probability 2/ln n, weight ln n: parity only | 1, 2 | 0.1555 / 0.1555 | no forest windows (model has none); control in the primes' windows: -0.0003 | 0.8445 | 0 |
| primes | Λ(n): every congruence condition of the primes | 1, 2 | 0.1561 / 0.1561 | 0.4235 / 0.4353 | 0.4204 | 312 |
| primes ≡ 1 (mod 4), weight 2Λ | one extra congruence condition | 1, 2, 4 | 0.1561 / 0.1560 | 0.3817 / 0.3941 | 0.4622 | 363 |

Removing congruence structure removes forest. Cramér's model has one atom (q = 1) and no forest. Restricting it to odd
numbers adds the atom at ½, and nothing more. The primes have the whole Farey forest. Read in the primes' forest windows,
the two Cramér sets give forests of +0.0001 and -0.0003: that is a
measurement, not just the model.

Adding a congruence condition does **not** add forest. The primes ≡ 1 (mod 4) live on one class mod 4, so
S(α + ¼) = i·S(α): the spikes at ¼ and ¾ are exact copies of the atoms at 0 and ½ and count as atoms. The set then has the
same atoms as the primes (15.6%) but a smaller forest (38.2% against
42.4%; model 39.4% against 43.5%) and a larger sea. Filing
the q = 4 copies under the forest instead would give 46.0%, which is the artefact a first draft of
this file reported. The forest readings here use the uniform sea, so they share the primes'
-2.7% bias; the comparison between sets is unaffected.

In body-composition terms (small for gestational age: more water, less fat) only the first direction carries over: less
structure, less forest, more sea. The "large for gestational age" direction has no counterpart here. The measurement is
literal; the nutrition language is not.

## Literal versus analogy

| body water | on the circle | status |
|:---|:---|:---|
| total body water, conserved and split into compartments | Parseval: ΣΛ² split into atoms, forest and sea | LITERAL bookkeeping |
| dilution volume V = I/C | L_N = total energy / energy of one spike | LITERAL ratio |
| V = I/C as the measuring rule | compartment = (window reading − sea covered)/capture, a recovery correction | ANALOGY (same arithmetic form; nothing is given or diluted) |
| indicators with different spaces (D₂O, bromide/inulin, Evans blue) | windows around different sets of denominators; nested subtraction | LITERAL nesting |
| incomplete recovery of an indicator | a window catches 0.755–0.803 of its spike; one average capture errs by +3.8% on the atoms, with a sign that depends on the grid offset | ANALOGY |
| an indicator reaching beyond its space | the sea a window covers | ANALOGY |
| equilibration (the indicator mixes evenly) | uniform sea: fails near rationals, where the sea is the divisor sum of the residue-class variances | MEASURED failure (with an exact identity) |
| development: TBW 94 → 86 → 78 → 70 → 60%, extracellular ↓, intracellular ↑ to one year | atoms 2/L_N ↓, forest ≈ 0.44, sea ↑; leading-order limit 0, ½, ½ | ANALOGY (X is not time; compartment pairing is a choice) |
| growth aberration: small for gestational age (more water, less fat) | less congruence structure: less forest, more sea | LITERAL measurement, ANALOGY wording |
| growth aberration: large for gestational age (more fat) | an extra congruence condition adds atoms (copies), not forest | NO COUNTERPART |
| osmotic equilibrium across membranes | equidistribution of primes over reduced classes (V(q) is the imbalance) | ANALOGY (no flow, no membrane) |
| postnatal contraction of extracellular water, diuresis | none | NO COUNTERPART |
| 60-40-20 rule | no fixed fractions at finite X; only the model limit (0, ½, ½) is constant | NO |

## What this does and does not show

It shows that the energy budget of S can be measured with body-water bookkeeping, and that two corrections physiology
knows are needed here too. Each window needs its exact capture. The sea is not evenly mixed: near a/q it sits at the
divisor sum of the variances of primes in progressions, well below its average. With both handled, the atoms and the
forest match the no-fit spike model to within 0.3% at X ≥ 3·10^4. Removing congruence structure from the set
removes forest; adding a condition adds copies of the main spike instead.

It does not show any mechanism shared with physiology, and it gives no new information about Goldbach's problem. The
forest is the same kind of quantity as the major-arc mass of results.md (larger Q, narrower arcs, without q = 1, 2), and
the sea is the minor-arc obstacle of escapes.md.

## References

- A. Tobias, B. D. Ballard, S. S. Mohiuddin, "Physiology, Water Balance", StatPearls (2022), NCBI Bookshelf NBK541059,
  https://www.ncbi.nlm.nih.gov/books/NBK541059/ (the full text could not be fetched from this environment; the figures used are the standard ones it
  states for a 70 kg man).
- "Body Composition in the Fetus and Newborn: Effects of Intrauterine Growth Aberration", ch. 2 of *Nephrology and
  Fluid/Electrolyte Physiology: Neonatology Questions and Controversies*, via https://abdominalkey.com/body-composition-in-the-fetus-and-newborn-effects-of-intrauterine-growth-aberration/ (open part only).
- Variance of primes in progressions: Barban (1964); Davenport & Halberstam (1966); Montgomery (1970); Hooley (1975);
  Friedlander & Goldston (1996), as in fractals.md.

Run: 18 s, peak RSS 760 MB.
