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
| phi:1e-17 | 1 | 0 | inf | inf | 1.4e-16 | 1.4e-16 | - | 0 |
| phi:1e-8 | 1 | 0 | 5.0e-09 | 2.0e-10 | 8.2e-17 | 8.2e-17 | - | 0 |
| phi:1e-6 | 1 | 0 | 2.9e-11 | 1.6e-12 | 1.1e-16 | 1.1e-16 | - | 0 |

## 2. neg_binomial_2_lpmf(n | mu, phi)

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | inf | inf | 1.4e-16 | 1.4e-16 | 0 | 0 |
| phi:1e-8 | 1 | 0 | 7.9e-07 | 3.1e-08 | 8.2e-17 | 8.2e-17 | 0 | 0 |
| phi:1e-6 | 1 | 0 | 2.5e-09 | 1.4e-10 | 1.1e-16 | 1.1e-16 | 0 | 0 |

## 3. neg_binomial_2_log_lpmf(n | log(mu), phi)

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | inf | inf | 1.4e-16 | 1.4e-16 | 0 | 0 |
| phi:1e-8 | 1 | 0 | 7.9e-07 | 3.1e-08 | 8.2e-17 | 8.2e-17 | 2.4e-24 | 0 |
| phi:1e-6 | 1 | 0 | 2.5e-09 | 1.4e-10 | 1.1e-16 | 1.1e-16 | 4.3e-22 | 0 |

## 4a. neg_binomial_2_log_glm_lpmf, one observation

| config:decade | pts | thr | val abs | val rel | lsg abs | lsg rel | oth | ms max |
|---|---|---|---|---|---|---|---|---|
| phi:1e-17 | 1 | 0 | 7.0e-16 | 1.8e-17 | 1.7e-17 | 1.7e-17 | 0 | 0 |
| phi:1e-8 | 1 | 0 | 1.3e-12 | 5.1e-14 | 6.7e-17 | 6.7e-17 | 2.4e-24 | 0 |
| phi:1e-6 | 1 | 0 | 1.1e-13 | 6.0e-15 | 1.1e-16 | 1.1e-16 | 1.4e-14 | 0 |

## Summary (defect: abs and rel error both > 1e-8)

| tag | value defect from | grad defect from | worst val abs | worst lsg abs | thr/pts |
|---|---|---|---|---|---|
| NB | 1e-17 | - | inf | 1.4e-16 | 0/3 |
| NB2 | 1e-17 | - | inf | 1.4e-16 | 0/3 |
| NB2L | 1e-17 | - | inf | 1.4e-16 | 0/3 |
| GLM1 | - | - | 1.3e-12 | 1.1e-16 | 0/3 |

## Worst points

- NB val inf at phi:1e-17 (1, 1e-17, 1e-17)
- NB lsg 1.4e-16 at phi:1e-17 (1, 1e-17, 1e-17)
- NB2 val inf at phi:1e-17 (1, 1, 1e-17)
- NB2 lsg 1.4e-16 at phi:1e-17 (1, 1, 1e-17)
- NB2L val inf at phi:1e-17 (1, 0, 1e-17)
- NB2L lsg 1.4e-16 at phi:1e-17 (1, 0, 1e-17)
- GLM1 val 1.3e-12 at phi:1e-8 (1000, 1, 0, 6.908, 1e-08)
- GLM1 lsg 1.1e-16 at phi:1e-6 (100, 1, 0, 4.605, 1e-06)
