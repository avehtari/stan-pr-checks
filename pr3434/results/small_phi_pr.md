## Reference certification

Worst relative disagreement per tag: route A (analytic, 90 digits)
against route B (140 digits, mp.diff), and A against C (second
formula, value only; count of points in parentheses).

| tag | value A-B | grad A-B | value A-C |
|---|---|---|---|
| NB | 0 | 0 | - |
| NB2 | 0 | 5.0e-147 | - |
| NB2L | 0 | 0 | - |
| GLM1 | 0 | 0 | - |

## 1. neg_binomial_lpmf(n | alpha, beta), beta = alpha / mu

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | 7.0e-16 | 1.8e-17 | 1.4e-16 | 1.4e-16 | - | 0 |
| phi:1e-8 | 1 | 0 | 5.2e-16 | 2.1e-17 | 8.2e-17 | 8.2e-17 | - | 0 |
| phi:1e-6 | 1 | 0 | 2.5e-16 | 1.4e-17 | 1.1e-16 | 1.1e-16 | - | 0 |

## 2. neg_binomial_2_lpmf(n | mu, phi)

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | 7.0e-16 | 1.8e-17 | 1.4e-16 | 1.4e-16 | 0 | 0 |
| phi:1e-8 | 1 | 0 | 5.3e-13 | 2.1e-14 | 8.2e-17 | 8.2e-17 | 0 | 0 |
| phi:1e-6 | 1 | 0 | 6.0e-14 | 3.3e-15 | 1.1e-16 | 1.1e-16 | 0 | 0 |

## 3. neg_binomial_2_log_lpmf(n | log(mu), phi)

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | 7.0e-16 | 1.8e-17 | 1.4e-16 | 1.4e-16 | 0 | 0 |
| phi:1e-8 | 1 | 0 | 3.8e-13 | 1.5e-14 | 8.2e-17 | 8.2e-17 | 2.4e-24 | 0 |
| phi:1e-6 | 1 | 0 | 2.8e-13 | 1.5e-14 | 1.1e-16 | 1.1e-16 | 4.3e-22 | 0 |

## 4a. neg_binomial_2_log_glm_lpmf, one observation

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | 7.0e-16 | 1.8e-17 | 1.4e-16 | 1.4e-16 | 0 | 0 |
| phi:1e-8 | 1 | 0 | 5.2e-16 | 2.1e-17 | 8.2e-17 | 8.2e-17 | 2.4e-24 | 0 |
| phi:1e-6 | 1 | 0 | 2.5e-16 | 1.4e-17 | 1.1e-16 | 1.1e-16 | 4.3e-22 | 0 |

## Summary (defect: abs and rel error both > 1e-8)

| tag | value defect from | grad defect from | worst val abs | worst lsg abs | thr/pts |
|---|---|---|---|---|---|
| NB | - | - | 7.0e-16 | 1.4e-16 | 0/3 |
| NB2 | - | - | 5.3e-13 | 1.4e-16 | 0/3 |
| NB2L | - | - | 3.8e-13 | 1.4e-16 | 0/3 |
| GLM1 | - | - | 7.0e-16 | 1.4e-16 | 0/3 |

## Worst points

- NB val 7.0e-16 at phi:1e-17 (1, 1e-17, 1e-17)
- NB lsg 1.4e-16 at phi:1e-17 (1, 1e-17, 1e-17)
- NB2 val 5.3e-13 at phi:1e-8 (1000, 1000, 1e-08)
- NB2 lsg 1.4e-16 at phi:1e-17 (1, 1, 1e-17)
- NB2L val 3.8e-13 at phi:1e-8 (1000, 6.908, 1e-08)
- NB2L lsg 1.4e-16 at phi:1e-17 (1, 0, 1e-17)
- GLM1 val 7.0e-16 at phi:1e-17 (1, 1, 0, 0, 1e-17)
- GLM1 lsg 1.4e-16 at phi:1e-17 (1, 1, 0, 0, 1e-17)
