# Compartments: the body-water lens on the energy of S at N = 10^6

> Exploration at finite X, not a proof. Every number below is read from `cache/compartments.pkl` (grids M = 2^17 … 2^24,
> X = 10^4 … 2·10^6). Body-water physiology is used as a teacher for bookkeeping, not as a mechanism: nothing in S moves,
> mixes or equilibrates.

**What physiology does.** Total body water (TBW) is about 60% of an adult man's weight (42 L of 70 kg). Two thirds of it
is intracellular (28 L) and one third extracellular (14 L: 10.5 L interstitial and 3.5 L plasma) (StatPearls, *Physiology,
Water Balance*). Compartments are measured by **indicator dilution**: give a known amount I of an indicator that spreads
only through the compartment of interest, let it equilibrate, read its plasma concentration C, and set V = I/C. Heavy
water (D₂O, H₂¹⁸O) or antipyrine measures TBW; bromide, sucrose or inulin measures extracellular water; Evans blue
measures plasma. The rest follows by subtraction: intracellular = TBW − extracellular, interstitial = extracellular −
plasma, solids = body weight − TBW (Abdominal Key, ch. 2). The composition changes with growth. TBW falls from 94% of body
weight in the first trimester to 86% at 24 weeks, 78% at term and about 70% at one year, as solids accrete. Extracellular
water falls (60% → 44% → ≈ 28%) and intracellular water rises (26% → 34% → 42%). It also depends on growth conditions:
infants small for gestational age have more water and less fat; infants large for gestational age have more fat and
less water.

**What the circle has.** One conserved total, by Parseval: (1/M)Σ_k |S(k/M)|² = ΣΛ² exactly. The "indicators" are windows
|N(α − a/q)| ≤ ½ around the reduced fractions a/q of a set of denominators. A window "reads" the energy it covers; it
catches a known share (its capture) of the Dirichlet spike centred in it, and it also covers some sea. Three compartments:

- **atoms**: the spikes at q = 1, 2 (α = 0 and ½);
- **forest**: the spikes with 3 ≤ q ≤ q*(X), where q* is the denominator at which spike and sea are equal (fractals.md, M1);
- **sea**: the rest, i.e. the structureless part plus the spikes with q > q* that it drowns.

The no-fit spike model is S(a/q + β) ≈ μ(q)(ψ/φ(q))·v(β)/X, with v the Dirichlet kernel of energy X. So the φ(q) spikes of
denominator q carry ψ²μ²(q)/(Xφ(q)) together, and the shares are atoms = 2/L_N and forest = (Σ_{3≤q≤q*} μ²/φ)/L_N, with
L_N = XΣΛ²/ψ². **L_N is a dilution volume**: the total energy divided by the energy of one reference spike (the one at 0,
ψ²/X). At X = 10^6, L_N = 12.81: the energy of S is spread as if over 12.8 copies of its main spike.

## 0. Certificates

| check | result | detail |
|:---|:---|:---|
| Parseval, primes, X = 10000 | pass |  |
| Parseval, primes, X = 30000 | pass |  |
| Parseval, primes, X = 100000 | pass |  |
| Parseval, primes, X = 300000 | pass |  |
| Parseval, primes, X = 1000000 | pass |  |
| window energy = capture x model, q <= 30 at X = 10^6 | pass | max \|ratio - 1\| = 0.0004 |
| Parseval, primes, X = 2000000 | pass |  |
| Parseval, Cramér | pass |  |
| Parseval, odd Cramér | pass |  |
| Parseval, primes | pass |  |
| Parseval, primes ≡ 1 (mod 4), weight 2Λ | pass |  |

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

## B1. Compartments by indicator dilution (X = 10^6)

| indicator | what it reads | share of ΣΛ² |
|:---|:---|---:|
| whole circle (Parseval) | the conserved total ΣΛ² | 100% |
| windows at q = 1, 2 | raw reading | 12.54% |
| windows at 3 ≤ q ≤ q* | raw reading | 34.94% |

The readings are turned into compartments by the three-indicator solution. Each window reads its capture of its own
spikes plus the sea it covers; the rest of the circle holds the spikes' leakage plus the remaining sea:

| compartment | one average capture 0.7737, uniform sea | exact capture, uniform sea | exact capture, sea at X·FG(q) | exact capture, sea measured at the centres | no-fit model |
|:---|---:|---:|---:|---:|---:|
| atoms (q = 1, 2) | 16.20% | 15.61% | 15.61% | 15.61% | 15.61% |
| forest (3 ≤ q ≤ q*) | 42.11% | 42.35% | 43.78% | 43.64% | 43.53% |
| sea (the rest) | 41.69% | 42.04% | 40.61% | 40.75% | 40.86% |

Two corrections matter, and both have close counterparts in tracer practice:

1. **The indicator's distribution space.** A window's capture depends on where the fraction falls between grid points: at
   X = 10^6 a window holds 8 or 9 grid points and catches 0.755 or 0.803 of its spike, not the average 0.7737. Using one
   average capture misreads the atoms by +3.8%. The same kind of error happens when a
   tracer reaches more than its compartment (in standard tracer practice, not from the two sources above: heavy water
   exchanges with non-aqueous hydrogen, and bromide enters cells a little). This was also the A1 octave bug that the fractal-lens review caught.
2. **The equilibration assumption.** Plain dilution assumes the sea is uniformly mixed. It is not (next table). With a
   uniform sea the forest reads -2.7% against the model. With the sea level
   measured at the spike centres it reads +0.3%; with the no-fit X·FG(q) it reads
   +0.6%.

**The sea is thinner near the forest.** Sea level per unit α in the windows of each band of denominators, at X = 10^6:

| windows of | # q | measured in the windows | measured at the centres | X·FG(q) (no fit) | uniform sea |
|:---|---:|:---|:---|:---|:---|
| q ∈ [31, 100] | 42 | 1.72 X | 1.69 X | 1.38 X | 5.38 X |
| q ∈ [101, 255] | 96 | 2.82 X | 2.66 X | 2.37 X | 5.38 X |
| q ∈ [256, 512] | 157 | 3.56 X | 3.36 X | 3.14 X | 5.38 X |

Near a fraction with denominator q, the structureless part of S is set by the variance of the primes among the residue
classes mod q: about X·FG(q) = X(ln q − γ − ln 2π − Σ_{p|q} ln p/(p−1)) (Friedlander–Goldston form of the Barban–Davenport–Halberstam /
Montgomery–Hooley variance). That is well below the circle's average sea level. The level climbs back toward the
average as q grows, and X·FG(q) falls 12–20% below the window values. (Bands with q ≤ 30 are left out: there the sea is too
small a fraction of the window to measure.)

## B2. Development: the composition of S as X grows

| X | L_N | q* | atoms: measured / model | forest: measured / model | sea: measured / model | forest with a uniform sea |
|---:|---:|---:|---:|---:|---:|---:|
| 1e4 | 8.15 | 75 | 0.2447 / 0.2453 | 0.4440 / 0.4484 | 0.3114 / 0.3063 | 0.4184 |
| 3e4 | 9.29 | 115 | 0.2154 / 0.2153 | 0.4412 / 0.4422 | 0.3434 / 0.3424 | 0.4197 |
| 1e5 | 10.49 | 190 | 0.1906 / 0.1907 | 0.4355 / 0.4352 | 0.3739 / 0.3741 | 0.4177 |
| 3e5 | 11.61 | 303 | 0.1723 / 0.1722 | 0.4359 / 0.4348 | 0.3919 / 0.3930 | 0.4208 |
| 1e6 | 12.81 | 513 | 0.1561 / 0.1561 | 0.4364 / 0.4353 | 0.4075 / 0.4086 | 0.4235 |
| 2e6 | 13.50 | 698 | 0.1481 / 0.1481 | 0.4367 / 0.4358 | 0.4152 / 0.4161 | 0.4250 |

For comparison, the body (percent of body weight; fetal and infant values from the Abdominal Key chapter, adult from
StatPearls; "≈" marks a value obtained by subtraction):

| stage | total water | extracellular | intracellular | solids |
|:---|---:|:---|:---|---:|
| first trimester | 94% | — | — | ≈ 6% |
| 24 weeks | 86% | 60% | 26% | 14% |
| term | 78% | 44% | 34% | 22% |
| 1 year | 70% | ≈ 28% | 42% | 30% |
| adult man | 60% | 20% | 40% | ≈ 40% |

The bookkeeping matches closely: one conserved total redistributing among compartments as a growth parameter
increases, toward a limit. The **atoms** fall like 2/L_N ≈ 2/(ln X − 1), from 24.5%
at X = 10^4 to 14.8% at 2·10^6. Their energy 2ψ²/X grows only like X, while
the total grows like X ln X. The **forest** stays near 0.44. The **sea** grows. As X → ∞ the shares tend to **0, ½, ½**:
q* grows like √X up to logarithms, so the forest holds (ln q* + 1.33 − 2)/(ln X − 1) → ½ of the energy. The visible forest
ends up holding half the logarithm. The approach is slow (leading-order law):

| X | atoms | forest | sea |
|:---|---:|---:|---:|
| 10^6 | 0.156 | 0.435 | 0.409 |
| 10^12 | 0.075 | 0.450 | 0.475 |
| 10^20 | 0.044 | 0.464 | 0.492 |
| 10^40 | 0.022 | 0.478 | 0.500 |
| 10^100 | 0.009 | 0.489 | 0.502 |
| 10^1000 | 0.001 | 0.498 | 0.501 |

**Where the analogy stops.** The coherent share of S (atoms plus forest) *falls* with X, from 68.9% to 58.5%
(limit ½), while a fetus *gains* solids. Which compartment of S "corresponds" to extracellular or intracellular water is
a choice with no physical basis. Only the bookkeeping and the derived limit are literal. X is a parameter, not time, and
there is no counterpart to birth or to the postnatal contraction of extracellular water.

## B3. "Growth aberration": the arithmetic of the set sets its composition (X = 10^6)

Four sets with the same total weight (≈ ψ) but different congruence structure. The forest windows run to the primes' q*
(512), placed wherever each set's own no-fit model has spikes. That model is the set's weight spread
evenly over the residue classes it allows mod lcm(D, q), summed with e(ab/q).

| set (X = 10^6, same total weight ≈ ψ) | what it keeps | atoms: measured / model | forest: measured / model | sea | # denominators with spikes ≤ q* |
|:---|:---|---:|---:|---:|---:|
| Cramér | random n ≥ 3 with probability 1/ln n, weight ln n: no congruence structure | 0.0780 / 0.0779 | 0.0000 / 0.0000 | 0.9220 | 0 |
| odd Cramér | random odd n with probability 2/ln n, weight ln n: parity only | 0.1555 / 0.1555 | 0.0000 / 0.0000 | 0.8445 | 0 |
| primes | Λ(n): every congruence condition of the primes | 0.1561 / 0.1561 | 0.4235 / 0.4353 | 0.4204 | 312 |
| primes ≡ 1 (mod 4), weight 2Λ | one extra congruence condition | 0.0781 / 0.0780 | 0.4598 / 0.4721 | 0.4621 | 364 |

Less congruence structure means less forest. Cramér's model has one spike (q = 1). Restricting it to odd numbers adds the
spike at ½, and nothing more. The primes have the whole Farey forest. Keeping only the primes ≡ 1 (mod 4) adds spikes at
q = 4·(odd squarefree): 364 denominators instead of 312. The forest grows to
46.0% of a total energy that has doubled, so the atoms' share halves. These forest
readings use the uniform sea, so they share the primes' -2.7% bias; the
comparison between sets is unaffected. In body-composition terms (SGA: more water, less fat; LGA: more fat, less water)
the sea plays "water" and the forest plays "solids". That pairing is a reading device: the measurement is literal, the
nutrition language is not.

## Literal versus analogy

| body water | on the circle | status |
|:---|:---|:---|
| total body water, conserved and split into compartments | Parseval: ΣΛ² split into atoms, forest and sea | LITERAL bookkeeping |
| dilution V = I/C | compartment = (window reading − sea covered)/capture; L_N = total/(energy of one spike) | LITERAL algebra |
| indicators with different spaces (D₂O, bromide/inulin, Evans blue) | windows around different sets of denominators; nested subtraction | LITERAL nesting |
| indicator corrections (tracer outside its space) | per-window capture 0.755/0.803 instead of 0.7737; otherwise +3.8% on the atoms | LITERAL, same kind of correction |
| equilibration (the tracer mixes evenly) | uniform sea: fails near rationals, where the sea is X·FG(q) | MEASURED failure |
| development: TBW 94 → 86 → 78 → 70 → 60%, extracellular ↓, intracellular ↑ | atoms 2/L_N ↓, forest ≈ 0.44, sea ↑; limits 0, ½, ½ | ANALOGY (X is not time; compartment pairing is a choice) |
| growth aberration (SGA more water, LGA more fat) | less / more congruence structure, less / more forest | LITERAL measurement, ANALOGY wording |
| osmotic equilibrium across membranes | equidistribution of primes over reduced classes (V(q) is the imbalance) | ANALOGY (no flow, no membrane) |
| postnatal contraction of extracellular water, diuresis | none | NO COUNTERPART |
| 60-40-20 rule | no fixed fractions at finite X; only the limit (0, ½, ½) is constant | NO |

## What this does and does not show

It shows that the energy budget of S can be measured the way a physiologist measures body water, and that the
corrections physiology learned the hard way are needed here too. The windows need their exact capture. The sea is not
evenly mixed: it is thinner near the rationals, at the level of the variance of primes in progressions. Once both are
handled, the atoms and the forest match the no-fit spike model to within 0.3% at X ≥ 3·10^4 (1.0% at
10^4). As X grows the
composition drifts toward the limit (0, ½, ½) like 1/ln X.

It does not show any mechanism shared with physiology, and it gives no new information about Goldbach's problem: the
forest's energy is the major-arc energy of results.md, and the sea is the minor-arc obstacle of escapes.md.

## References

- A. Tobias, B. D. Ballard, S. S. Mohiuddin, "Physiology, Water Balance", StatPearls (2022), NCBI Bookshelf NBK541059,
  https://www.ncbi.nlm.nih.gov/books/NBK541059/ (the full text could not be fetched from this environment; the figures used are the standard ones it
  states for a 70 kg man).
- "Body Composition in the Fetus and Newborn: Effects of Intrauterine Growth Aberration", ch. 2 of *Nephrology and
  Fluid/Electrolyte Physiology: Neonatology Questions and Controversies*, via https://abdominalkey.com/body-composition-in-the-fetus-and-newborn-effects-of-intrauterine-growth-aberration/ (open part only).
- Variance of primes in progressions: Barban (1964); Davenport & Halberstam (1966); Montgomery (1970); Hooley (1975);
  Friedlander & Goldston (1996), as in fractals.md.

Run: 0.1 min, peak RSS 759 MB.
