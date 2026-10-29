Columns: worst absolute error (relative error at that point);
max rel: the worst relative error over all points.

| function | tree | points | throw | NaN | inf | worst value (rel) | worst log-scale grad (rel) | max rel value | max rel grad |
|---|---|---|---|---|---|---|---|---|---|
| beta_binomial_lpmf | dev | 13177 | 0 | 0 | 0 | 2.2e+02 (4.1e-04) | 7.7e+03 (1.5e+06) | 3.1e+01 | 3.2e+06 |
| beta_binomial_lpmf | pr | 13177 | 0 | 0 | 0 | 7.1e-09 (2.8e-16) | 3.5e-10 (2.3e-14) | 7.0e+00 | 8.1e-05 |
| beta_neg_binomial_lpmf | dev | 159 | 0 | 0 | 0 | 1.0e+01 (3.8e+00) | 1.3e+01 (7.5e+00) | 3.8e+00 | 7.5e+00 |
| beta_neg_binomial_lpmf | pr | 159 | 0 | 0 | 0 | 2.6e-02 (4.9e-17) | 8.2e-02 (2.8e-16) | 3.9e-09 | 3.1e-08 |
| neg_binomial_lpmf | dev | 42 | 0 | 0 | 0 | 2.3e-13 (8.1e-14) | 3.0e+00 (0) | 8.1e-14 | 5.0e-01 |
| neg_binomial_lpmf | pr | 42 | 0 | 0 | 0 | 2.3e-13 (8.1e-14) | 1.0e-14 (1.9e-16) | 8.1e-14 | 2.7e-10 |
| neg_binomial_2_lpmf | dev | 310 | 0 | 0 | 0 | 5.6e-13 (1.6e-13) | 2.7e+01 (0) | 1.6e-13 | 5.7e-02 |
| neg_binomial_2_lpmf | pr | 310 | 0 | 0 | 0 | 5.6e-13 (1.6e-13) | 1.7e-14 (0) | 1.6e-13 | 1.0e-09 |
| neg_binomial_2_log_lpmf | dev | 310 | 0 | 0 | 0 | 9.6e-13 (2.8e-13) | 2.6e+01 (0) | 2.8e-13 | 4.8e-02 |
| neg_binomial_2_log_lpmf | pr | 310 | 0 | 0 | 0 | 9.6e-13 (2.8e-13) | 2.0e-13 (0) | 2.8e-13 | 5.9e-09 |
| neg_binomial_2_log_glm_lpmf | dev | 63 | 0 | 0 | 0 | 2.3e+01 (4.2e-01) | 7.1e+00 (0) | 1.0e+00 | 4.0e-02 |
| neg_binomial_2_log_glm_lpmf | pr | 63 | 0 | 0 | 0 | 4.5e-13 (2.7e-15) | 2.8e-14 (2.0e-09) | 1.1e-14 | 2.0e-09 |
| dirichlet_multinomial_lpmf | dev | 35 | 0 | 0 | 0 | 1.0e-12 (1.6e-14) | 2.8e+00 (4.2e-01) | 8.3e-14 | 1.0e+00 |
| dirichlet_multinomial_lpmf | pr | 35 | 0 | 0 | 0 | 1.0e-12 (1.6e-14) | 8.1e-15 (7.8e-12) | 8.3e-14 | 2.0e-10 |
| lkj_corr_lpdf | dev | 28 | 0 | 0 | 0 | 3.9e+00 (7.5e-02) | 8.2e+00 (1.5e-14) | 7.5e-02 | 1.0e+00 |
| lkj_corr_lpdf | pr | 28 | 0 | 0 | 0 | 1.5e-01 (2.7e-16) | 1.4e-01 (2.6e-16) | 7.3e-16 | 4.0e-16 |
| lkj_corr_cholesky_lpdf | dev | 28 | 0 | 0 | 0 | 3.9e+00 (7.5e-02) | 8.4e+00 (1.5e-14) | 7.5e-02 | 1.0e+00 |
| lkj_corr_cholesky_lpdf | pr | 28 | 0 | 0 | 0 | 2.4e-02 (4.5e-17) | 2.4e-02 (4.4e-17) | 2.8e-16 | 1.3e-16 |
| yule_simon_lccdf | dev | 28 | 0 | 0 | 0 | 7.6e+00 (1.1e-01) | 2.1e+00 (4.2e-01) | 1.1e-01 | 1.0e+00 |
| yule_simon_lccdf | pr | 28 | 0 | 0 | 0 | 1.7e-13 (1.5e-16) | 8.2e-15 (1.4e-16) | 2.2e-16 | 1.6e-16 |
| yule_simon_cdf | dev | 28 | 0 | 0 | 0 | 2.8e-15 (2.8e-15) | 5.1e-15 (0) | 2.8e-15 | 5.5e-12 |
| yule_simon_cdf | pr | 28 | 0 | 0 | 0 | 5.5e-17 (5.5e-17) | 7.5e-18 (7.6e-16) | 5.5e-17 | 1.5e-15 |
| yule_simon_lcdf | dev | 28 | 0 | 0 | 0 | 2.8e-15 (2.8e-05) | 5.1e-15 (0) | 1.0e+00 | 5.5e-12 |
| yule_simon_lcdf | pr | 28 | 0 | 0 | 0 | 7.4e-18 (7.5e-16) | 7.4e-18 (7.5e-16) | 1.0e+00 | 1.5e-15 |
| beta_binomial_cdf | dev | 81 | 4 | 0 | 0 | 2.1e-01 (4.0e-01) | 1.9e+01 (0) | 1.7e+34 | 7.6e+00 |
| beta_binomial_cdf | pr | 81 | 4 | 0 | 0 | 4.1e-08 (7.7e-08) | 1.1e-06 (6.5e-07) | 1.6e+22 | 1.5e-05 |
| beta_binomial_lcdf | dev | 81 | 4 | 16 | 0 | inf (inf) | 4.1e+04 (1.7e+03) | inf | 1.7e+03 |
| beta_binomial_lcdf | pr | 81 | 4 | 0 | 0 | 7.7e-08 (1.2e-07) | 2.3e-06 (7.2e-07) | 1.8e-07 | 7.7e-07 |
| beta_binomial_lccdf | dev | 81 | 4 | 0 | 0 | 3.7e-01 (4.9e-01) | 2.2e+01 (5.6e+00) | 1.6e+34 | 5.6e+00 |
| beta_binomial_lccdf | pr | 81 | 4 | 0 | 0 | 1.8e-07 (9.7e-09) | 2.1e-06 (7.3e-08) | 1.6e+22 | 1.5e-05 |
| beta_neg_binomial_cdf | dev | 104 | 11 | 0 | 0 | 3.3e-01 (8.1e+01) | 1.5e+01 (1.1e+03) | 2.2e+09 | 1.1e+03 |
| beta_neg_binomial_cdf | pr | 104 | 11 | 0 | 0 | 2.0e-09 (9.4e+00) | 2.0e-06 (1.4e-05) | 4.4e+08 | 5.8e-03 |
| beta_neg_binomial_lcdf | dev | 101 | 10 | 0 | 0 | 2.1e+01 (4.1e-01) | 4.3e+01 (1.3e+00) | 1.0e+00 | 1.0e+01 |
| beta_neg_binomial_lcdf | pr | 101 | 8 | 0 | 0 | 2.0e+01 (3.8e-01) | 6.9e+01 (2.4e+00) | 1.0e+00 | 2.4e+00 |
| beta_neg_binomial_lccdf | dev | 101 | 8 | 0 | 0 | 2.9e-01 (7.0e+01) | 1.3e+01 (2.2e+01) | 2.2e+09 | 8.2e+02 |
| beta_neg_binomial_lccdf | pr | 101 | 8 | 0 | 0 | 2.0e-09 (9.4e+00) | 3.4e-06 (5.1e-07) | 4.4e+08 | 5.8e-03 |
| student_t_lpdf | dev | 230 | 0 | 0 | 0 | 1.6e-01 (7.4e-05) | 1.4e-01 (0) | 1.5e-02 | 2.5e+02 |
| student_t_lpdf | pr | 230 | 0 | 0 | 0 | 1.3e-11 (1.2e-16) | 4.5e-12 (2.4e-12) | 3.2e-15 | 5.0e-08 |
| multi_student_t_cholesky_lpdf | dev | 1013 | 0 | 0 | 0 | 1.8e-03 (6.3e-04) | 8.9e-05 (0) | 6.3e-04 | 1.3e-01 |
| multi_student_t_cholesky_lpdf | pr | 1013 | 0 | 0 | 0 | 8.1e-15 (4.2e-16) | 4.5e-15 (7.4e-10) | 1.6e-15 | 1.5e-09 |
| multi_student_t_lpdf | dev | 1013 | 0 | 0 | 0 | 1.6e-03 (5.6e-04) | 2.2e-07 (2.0e-01) | 5.6e-04 | 2.0e-01 |
| multi_student_t_lpdf | pr | 1013 | 0 | 0 | 0 | 5.6e-15 (3.9e-16) | 6.3e-15 (5.0e-15) | 1.6e-15 | 2.4e-09 |
| yule_simon_lpmf | dev | 1013 | 0 | 0 | 0 | 3.9e-03 (3.9e+09) | 3.6e-03 (0) | 3.9e+09 | 3.2e-04 |
| yule_simon_lpmf | pr | 1013 | 0 | 0 | 0 | 6.9e-10 (2.1e-16) | 1.7e-10 (1.9e-16) | 1.9e-03 | 3.0e-11 |
| lbeta (var), d/dx | dev | 275 | 0 | 0 | 0 | - | 6.3e-01 (4.2e-01) | - | 1.3e+00 |
| lbeta (var), d/dx | pr | 275 | 0 | 0 | 0 | - | 1.0e-14 (1.2e-16) | - | 1.9e-13 |

Worst points (region, arguments):

- beta_binomial_lpmf, dev: value at MC_d (1e+06, 1e+06, 1.144e+18, 7.972e+17); grad at MC_d (9.161e+05, 1e+06, 1.079e+18, 9.891e+16)
- beta_binomial_lpmf, pr: value at SEAM (1e+06, 1e+06, 1.572, 4.66e+16); grad at UN_c (1e+06, 1e+06, 1.279e+11, 1.926e+09)
- beta_neg_binomial_lpmf, dev: value at ab:1e15 (5, 5, 1e+15, 2e+15); grad at ab:1e15 (5, 5, 1e+15, 2e+15)
- beta_neg_binomial_lpmf, pr: value at all:1e15 (57, 1e+15, 1e+15, 1e+15); grad at all:1e15 (0, 1e+15, 1e+15, 1e+15)
- neg_binomial_lpmf, dev: value at mu50:1e6 (50, 1e+06, 2e+04); grad at mu3:1e15 (3, 1e+15, 3.333e+14)
- neg_binomial_lpmf, pr: value at mu50:1e6 (50, 1e+06, 2e+04); grad at mu3:1e8 (57, 1e+08, 3.333e+07)
- neg_binomial_2_lpmf, dev: value at dense:1e12 (57, 50, 3.047e+12); grad at dense:1e15 (1, 3, 3.778e+15)
- neg_binomial_2_lpmf, pr: value at dense:1e12 (57, 50, 3.047e+12); grad at dense:1e15 (57, 50, 1e+15)
- neg_binomial_2_log_lpmf, dev: value at dense:1e15 (57, 3.912, 2.608e+15); grad at dense:1e15 (1, 1.099, 3.778e+15)
- neg_binomial_2_log_lpmf, pr: value at dense:1e15 (57, 3.912, 2.608e+15); grad at dense:1e15 (57, 3.912, 2.096e+15)
- neg_binomial_2_log_glm_lpmf, dev: value at five:1e15 (0, 1, 5, 3, 57, 0.5, -0.25, 1.25, 0, 2, 1, 0.75, 1e+15); grad at five:1e15 (2, 4, 3, 5, 6, 0.5, -0.25, 1.25, 0, 2, 1, 0.75, 1e+15)
- neg_binomial_2_log_glm_lpmf, pr: value at x1.37:1e10 (57, 0.5, -0.2, 0.75, 1.37e+10); grad at one:1e8 (57, 0.5, 1, 0.75, 1e+08)
- dirichlet_multinomial_lpmf, dev: value at k:1e15 (57, 3, 10, 2e+14, 3e+14, 5e+14); grad at k:1e15 (5, 0, 0, 2e+14, 3e+14, 5e+14)
- dirichlet_multinomial_lpmf, pr: value at k:1e15 (57, 3, 10, 2e+14, 3e+14, 5e+14); grad at k:1e8 (57, 3, 10, 2e+07, 3e+07, 5e+07)
- lkj_corr_lpdf, dev: value at I:1e15 (0, 0, 0, 3.3e+15); grad at K3:1e15 (0.3, -0.2, 0.1, 3.3e+15)
- lkj_corr_lpdf, pr: value at K3:1e15 (0.3, -0.2, 0.1, 3.3e+15); grad at K3:1e15 (0.3, -0.2, 0.1, 3.3e+15)
- lkj_corr_cholesky_lpdf, dev: value at I:1e15 (0, 1, 0, 0, 1, 3.3e+15); grad at K3:1e15 (0.3, 0.9539, -0.2, 0.1677, 0.9653, 3.3e+15)
- lkj_corr_cholesky_lpdf, pr: value at K3:1e15 (0.3, 0.9539, -0.2, 0.1677, 0.9653, 3.3e+15); grad at K3:1e15 (0.3, 0.9539, -0.2, 0.1677, 0.9653, 3.3e+15)
- yule_simon_lccdf, dev: value at a:1e15 (2, 1e+15); grad at a:1e15 (5, 1e+15)
- yule_simon_lccdf, pr: value at a:1e10 (57, 1e+10); grad at a:1e8 (57, 1e+08)
- yule_simon_cdf, dev: value at a:1e10 (1, 1e+10); grad at a:1e8 (1, 1e+08)
- yule_simon_cdf, pr: value at a:1e4 (2, 1e+04); grad at a:1e2 (1, 100)
- yule_simon_lcdf, dev: value at a:1e10 (1, 1e+10); grad at a:1e8 (1, 1e+08)
- yule_simon_lcdf, pr: value at a:1e2 (1, 100); grad at a:1e2 (1, 100)
- beta_binomial_cdf, dev: value at ab3:1e15 (29, 117, 1e+15, 3e+15); grad at ab3:1e15 (5, 117, 1e+15, 3e+15)
- beta_binomial_cdf, pr: value at ab3:1e2 (29, 117, 100, 300); grad at ab3:1e2 (29, 117, 100, 300)
- beta_binomial_lcdf, dev: value at tied:1e4 (0, 117, 1e+04, 1e+04); grad at ab3:1e15 (5, 117, 1e+15, 3e+15)
- beta_binomial_lcdf, pr: value at ab3:1e2 (29, 117, 100, 300); grad at ab3:1e2 (29, 117, 100, 300)
- beta_binomial_lccdf, dev: value at ab3:1e15 (29, 117, 1e+15, 3e+15); grad at ab3:1e15 (29, 117, 1e+15, 3e+15)
- beta_binomial_lccdf, pr: value at ab3:1e4 (57, 117, 1e+04, 3e+04); grad at ab3:1e4 (57, 117, 1e+04, 3e+04)
- beta_neg_binomial_cdf, dev: value at ab:1e15 (0, 5, 1e+15, 2e+15); grad at ab:1e15 (0, 5, 1e+15, 2e+15)
- beta_neg_binomial_cdf, pr: value at b:1e4 (0, 3, 5, 1e+04); grad at r:1e1 (0, 10, 5, 2)
- beta_neg_binomial_lcdf, dev: value at all:1e2 (0, 100, 100, 100); grad at all:1e2 (0, 100, 100, 100)
- beta_neg_binomial_lcdf, pr: value at all:1e2 (0, 100, 100, 100); grad at all:1e2 (0, 100, 100, 100)
- beta_neg_binomial_lccdf, dev: value at ab:1e15 (0, 5, 1e+15, 2e+15); grad at ab:1e15 (5, 5, 1e+15, 2e+15)
- beta_neg_binomial_lccdf, pr: value at b:1e4 (0, 3, 5, 1e+04); grad at b:1e1 (57, 3, 5, 10)
- student_t_lpdf, dev: value at B8 (4.708, 8.023e+13, -0.0893, 0.07273); grad at B8 (3, 1e+14, 0, 2)
- student_t_lpdf, pr: value at B8 (4.678, 1.905e+09, -0.6296, 0.01144); grad at B8 (3.8, 1.13e+08, 1.32, 0.0145)
- multi_student_t_cholesky_lpdf, dev: value at B6 (1, 1, 1e+12); grad at B6 (1, 1, 1e+12)
- multi_student_t_cholesky_lpdf, pr: value at B5 (-4.678, 3.569, 4.69e+05); grad at B5 (4.085, 4.879, 6.028e+07)
- multi_student_t_lpdf, dev: value at B6 (1, 1, 1e+12); grad at B5 (4.645, 1.15, 9.845e+07)
- multi_student_t_lpdf, pr: value at B3 (1.058, 4.934, 1.074e+07); grad at B4 (4.633, 4.778, 291.2)
- yule_simon_lpmf, dev: value at B6 (1, 1e+12); grad at B6 (1, 1e+12)
- yule_simon_lpmf, pr: value at B5 (7.061e+05, 2.667e+07); grad at B5 (8.737e+05, 4.19e+07)
- lbeta (var), d/dx, dev: value at -; grad at B7 (1.5, 1e+14)
- lbeta (var), d/dx, pr: value at -; grad at B7 (82.27, 3.994e+05)

Per region: worst value and log-scale gradient error

| function | region | n | value dev | value pr | grad dev | grad pr | thr/nonfin dev | thr/nonfin pr |
|---|---|---|---|---|---|---|---|---|
| beta_binomial_lpmf | EQ_a | 1000 | 2.1e-10 | 2.1e-10 | 2.3e-14 | 2.4e-14 | 0/0 | 0/0 |
| beta_binomial_lpmf | EQ_b | 1000 | 5.9e-10 | 3.4e-10 | 2.7e-09 | 1.1e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | EQ_c | 1000 | 3.2e-04 | 3.1e-10 | 3.9e-03 | 1.9e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | EQ_d | 1000 | 1.6e+02 | 2.2e-10 | 6.1e+03 | 2.3e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | UN_a | 1000 | 3.3e-10 | 3.3e-10 | 2.1e-14 | 2.1e-14 | 0/0 | 0/0 |
| beta_binomial_lpmf | UN_b | 1000 | 2.8e-10 | 3.0e-10 | 4.6e-09 | 1.7e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | UN_c | 1000 | 1.1e-04 | 2.3e-09 | 3.6e-03 | 3.5e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | UN_d | 1000 | 6.7e+01 | 2.1e-09 | 6.1e+03 | 2.9e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | MC_a | 1000 | 2.7e-10 | 2.7e-10 | 3.9e-14 | 4.0e-14 | 0/0 | 0/0 |
| beta_binomial_lpmf | MC_b | 1000 | 4.4e-10 | 3.0e-10 | 5.6e-09 | 9.2e-11 | 0/0 | 0/0 |
| beta_binomial_lpmf | MC_c | 1000 | 1.9e-04 | 5.5e-10 | 5.1e-03 | 2.9e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | MC_d | 1000 | 2.2e+02 | 3.6e-10 | 7.7e+03 | 2.6e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | SEAM | 1000 | 7.1e-09 | 7.1e-09 | 2.0e+03 | 3.0e-10 | 0/0 | 0/0 |
| beta_binomial_lpmf | REP | 177 | 8.1e+01 | 4.3e-14 | 2.4e+02 | 2.3e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e2 | 5 | 1.8e-13 | 3.3e-14 | 1.2e-13 | 1.1e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e2 | 5 | 6.8e-14 | 3.1e-15 | 1.9e-13 | 6.4e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e2 | 4 | 3.8e-14 | 2.3e-14 | 1.1e-13 | 1.5e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e2 | 5 | 1.2e-13 | 4.0e-15 | 1.5e-13 | 5.5e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e2 | 5 | 1.9e-13 | 1.4e-14 | 1.5e-13 | 1.2e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e4 | 5 | 3.0e-11 | 1.6e-12 | 3.4e-11 | 7.4e-13 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e4 | 5 | 6.0e-13 | 1.7e-15 | 2.9e-11 | 5.1e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e4 | 4 | 7.5e-14 | 1.8e-14 | 5.7e-11 | 3.2e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e4 | 5 | 3.0e-11 | 5.7e-13 | 2.1e-11 | 1.0e-13 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e4 | 5 | 4.5e-11 | 4.4e-14 | 9.7e-11 | 2.5e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e6 | 5 | 8.9e-10 | 4.1e-10 | 1.8e-09 | 7.0e-11 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e6 | 5 | 6.5e-10 | 2.9e-15 | 1.5e-09 | 3.9e-11 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e6 | 4 | 4.0e-14 | 4.0e-14 | 5.6e-10 | 1.3e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e6 | 5 | 3.5e-09 | 4.0e-11 | 3.2e-09 | 8.5e-11 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e6 | 5 | 5.5e-09 | 6.2e-14 | 5.1e-09 | 1.8e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e8 | 5 | 2.8e-07 | 4.3e-08 | 5.9e-07 | 6.0e-09 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e8 | 5 | 1.9e-07 | 4.9e-15 | 3.6e-07 | 4.4e-09 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e8 | 4 | 1.5e-13 | 1.5e-13 | 1.9e-07 | 1.2e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e8 | 5 | 1.5e-07 | 2.1e-08 | 4.7e-07 | 1.3e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e8 | 5 | 3.4e-07 | 5.5e-14 | 6.1e-07 | 1.3e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e10 | 4 | 1.8e-05 | 5.9e-07 | 5.2e-05 | 3.8e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e10 | 4 | 2.5e-14 | 6.8e-15 | 3.5e-05 | 1.6e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e10 | 4 | 2.9e-13 | 6.6e-14 | 1.8e-05 | 1.4e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e10 | 4 | 2.4e-05 | 5.6e-14 | 1.8e-05 | 6.0e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e10 | 5 | 6.0e-05 | 2.6e-13 | 7.1e-05 | 3.8e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e12 | 4 | 6.2e-03 | 4.8e-05 | 4.9e-03 | 3.0e-05 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e12 | 4 | 3.8e-14 | 3.9e-15 | 3.7e-03 | 1.6e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e12 | 4 | 1.7e-13 | 5.4e-14 | 2.2e-03 | 4.2e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e12 | 4 | 3.0e-03 | 5.3e-14 | 2.0e-03 | 4.0e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e12 | 5 | 7.8e-03 | 4.8e-13 | 4.0e-03 | 7.1e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | all:1e15 | 4 | 4.8e+00 | 2.6e-02 | 6.2e+00 | 8.2e-02 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | r:1e15 | 4 | 2.8e-14 | 1.0e-14 | 2.0e+00 | 1.2e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | a:1e15 | 4 | 2.0e-13 | 2.0e-13 | 7.1e+00 | 3.3e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | b:1e15 | 4 | 2.8e+00 | 4.8e-14 | 3.0e+00 | 2.0e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lpmf | ab:1e15 | 5 | 1.0e+01 | 6.9e-14 | 1.3e+01 | 2.9e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e2 | 4 | 1.6e-14 | 1.6e-14 | 8.9e-14 | 2.3e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e2 | 2 | 9.3e-15 | 9.3e-15 | 6.1e-14 | 3.6e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e4 | 4 | 2.0e-14 | 2.0e-14 | 1.7e-11 | 4.3e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e4 | 2 | 2.5e-14 | 2.5e-14 | 2.8e-11 | 1.7e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e6 | 4 | 2.2e-13 | 2.2e-13 | 1.0e-09 | 3.5e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e6 | 2 | 2.3e-13 | 2.3e-13 | 1.5e-09 | 5.3e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e8 | 4 | 3.8e-14 | 3.8e-14 | 2.8e-07 | 1.0e-14 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e8 | 2 | 1.3e-13 | 1.3e-13 | 3.1e-07 | 3.3e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e10 | 4 | 6.7e-14 | 6.7e-14 | 2.3e-05 | 1.8e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e10 | 2 | 4.1e-14 | 4.1e-14 | 1.8e-05 | 3.8e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e12 | 4 | 9.5e-14 | 9.5e-14 | 1.7e-03 | 1.8e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e12 | 2 | 2.3e-13 | 2.3e-13 | 1.3e-03 | 8.1e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu3:1e15 | 4 | 3.1e-14 | 3.1e-14 | 3.0e+00 | 7.9e-15 | 0/0 | 0/0 |
| neg_binomial_lpmf | mu50:1e15 | 2 | 4.0e-14 | 4.0e-14 | 2.1e+00 | 2.4e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e2 | 4 | 2.3e-14 | 2.3e-14 | 1.0e-13 | 1.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e2 | 2 | 9.3e-15 | 9.3e-15 | 7.4e-14 | 2.2e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e4 | 4 | 1.5e-13 | 1.5e-13 | 1.5e-11 | 6.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e4 | 2 | 8.1e-14 | 8.1e-14 | 3.3e-11 | 1.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e6 | 4 | 1.0e-13 | 1.0e-13 | 1.0e-09 | 6.7e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e6 | 2 | 1.2e-13 | 1.2e-13 | 1.7e-09 | 1.2e-14 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e8 | 4 | 7.4e-14 | 7.4e-14 | 2.7e-07 | 1.4e-14 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e8 | 2 | 2.5e-13 | 2.5e-13 | 2.5e-07 | 3.3e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e10 | 4 | 7.0e-14 | 7.0e-14 | 3.6e-05 | 4.0e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e10 | 2 | 4.1e-14 | 4.1e-14 | 3.6e-05 | 1.2e-14 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e12 | 4 | 1.3e-13 | 1.3e-13 | 1.4e-09 | 3.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e12 | 2 | 2.3e-13 | 2.3e-13 | 1.0e-09 | 9.3e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu3:1e15 | 4 | 4.2e-13 | 4.2e-13 | 1.4e-12 | 9.4e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | mu50:1e15 | 2 | 4.0e-14 | 4.0e-14 | 1.0e-12 | 1.3e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e2 | 4 | 3.8e-14 | 3.8e-14 | 7.6e-14 | 6.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e4 | 4 | 3.2e-14 | 3.2e-14 | 1.9e-11 | 3.6e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e6 | 4 | 4.6e-14 | 4.6e-14 | 4.4e-09 | 1.3e-16 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e8 | 4 | 1.7e-13 | 1.7e-13 | 6.5e-07 | 6.6e-16 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e10 | 4 | 4.5e-13 | 4.5e-13 | 1.1e-07 | 6.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e12 | 4 | 2.6e-13 | 2.6e-13 | 1.1e-09 | 8.4e-16 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | x1.37:1e15 | 4 | 3.0e-14 | 3.0e-14 | 1.1e-12 | 3.4e-15 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | dense:1e10 | 80 | 5.5e-13 | 5.5e-13 | 1.4e-04 | 1.5e-14 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | dense:1e12 | 80 | 5.6e-13 | 5.6e-13 | 1.3e-02 | 1.3e-14 | 0/0 | 0/0 |
| neg_binomial_2_lpmf | dense:1e15 | 80 | 5.1e-13 | 5.1e-13 | 2.7e+01 | 1.7e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e2 | 4 | 1.8e-14 | 1.8e-14 | 9.0e-14 | 1.7e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e2 | 2 | 1.1e-14 | 1.1e-14 | 7.0e-14 | 8.9e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e4 | 4 | 3.0e-14 | 3.0e-14 | 1.7e-11 | 4.9e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e4 | 2 | 1.4e-13 | 1.4e-13 | 2.8e-11 | 1.9e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e6 | 4 | 1.1e-13 | 1.1e-13 | 1.0e-09 | 1.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e6 | 2 | 2.3e-13 | 2.3e-13 | 1.5e-09 | 5.3e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e8 | 4 | 7.9e-14 | 7.9e-14 | 2.8e-07 | 3.7e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e8 | 2 | 1.3e-13 | 1.3e-13 | 3.1e-07 | 6.7e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e10 | 4 | 7.5e-14 | 7.5e-14 | 2.3e-05 | 4.0e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e10 | 2 | 4.1e-14 | 4.1e-14 | 1.8e-05 | 1.2e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e12 | 4 | 9.3e-14 | 9.3e-14 | 1.7e-03 | 6.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e12 | 2 | 5.0e-14 | 5.0e-14 | 1.3e-03 | 7.6e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu3:1e15 | 4 | 3.8e-14 | 3.8e-14 | 3.0e+00 | 9.4e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | mu50:1e15 | 2 | 5.9e-14 | 5.9e-14 | 2.1e+00 | 9.0e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e2 | 4 | 3.8e-14 | 3.8e-14 | 4.4e-14 | 6.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e4 | 4 | 3.2e-14 | 3.2e-14 | 2.5e-11 | 1.5e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e6 | 4 | 4.6e-14 | 4.6e-14 | 3.5e-09 | 9.2e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e8 | 4 | 5.5e-14 | 5.5e-14 | 5.5e-07 | 1.6e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e10 | 4 | 4.5e-13 | 4.5e-13 | 1.9e-05 | 6.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e12 | 4 | 2.6e-13 | 2.6e-13 | 2.2e-03 | 8.0e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | x1.37:1e15 | 4 | 3.0e-14 | 3.0e-14 | 2.0e+00 | 3.4e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | dense:1e10 | 80 | 6.5e-13 | 6.5e-13 | 1.2e-04 | 1.1e-13 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | dense:1e12 | 80 | 5.3e-13 | 5.3e-13 | 9.4e-03 | 1.2e-13 | 0/0 | 0/0 |
| neg_binomial_2_log_lpmf | dense:1e15 | 80 | 9.6e-13 | 9.6e-13 | 2.6e+01 | 2.0e-13 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e2 | 3 | 1.1e-13 | 1.0e-14 | 1.4e-13 | 6.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e2 | 2 | 7.1e-13 | 9.7e-15 | 1.6e-13 | 9.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e4 | 3 | 4.4e-11 | 2.4e-14 | 2.0e-11 | 7.0e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e4 | 2 | 1.4e-10 | 5.4e-14 | 7.0e-11 | 4.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e6 | 3 | 1.4e-09 | 1.6e-13 | 5.5e-10 | 7.0e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e6 | 2 | 7.1e-09 | 2.1e-13 | 6.2e-09 | 6.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e8 | 3 | 5.0e-08 | 1.5e-13 | 4.3e-07 | 2.8e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e8 | 2 | 1.2e-06 | 6.5e-14 | 2.3e-07 | 7.6e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e10 | 3 | 4.6e-05 | 3.5e-13 | 3.6e-05 | 7.2e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e10 | 2 | 1.0e-04 | 4.1e-13 | 3.6e-05 | 3.1e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e12 | 3 | 4.3e-03 | 1.8e-14 | 1.4e-09 | 1.8e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e12 | 2 | 2.6e-02 | 4.3e-13 | 3.6e-03 | 5.0e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | one:1e15 | 3 | 4.0e+00 | 1.8e-13 | 1.4e-12 | 1.0e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | five:1e15 | 2 | 2.3e+01 | 6.1e-14 | 7.1e+00 | 9.4e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e2 | 4 | 2.1e-13 | 2.1e-14 | 6.5e-14 | 6.9e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e4 | 4 | 2.5e-11 | 1.7e-13 | 3.1e-11 | 1.9e-14 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e6 | 4 | 4.1e-09 | 3.0e-14 | 4.0e-09 | 7.5e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e8 | 4 | 8.5e-07 | 2.6e-13 | 4.5e-07 | 4.6e-16 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e10 | 4 | 5.0e-05 | 4.5e-13 | 1.1e-07 | 1.5e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e12 | 4 | 7.2e-03 | 1.5e-13 | 4.9e-03 | 7.3e-15 | 0/0 | 0/0 |
| neg_binomial_2_log_glm_lpmf | x1.37:1e15 | 4 | 1.5e+00 | 1.6e-13 | 1.1e-12 | 6.0e-16 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e2 | 5 | 3.0e-14 | 3.0e-14 | 1.3e-13 | 4.1e-15 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e4 | 5 | 8.9e-14 | 8.9e-14 | 2.6e-11 | 5.7e-15 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e6 | 5 | 1.5e-13 | 1.5e-13 | 1.7e-09 | 4.5e-15 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e8 | 5 | 1.0e-13 | 1.0e-13 | 3.1e-07 | 8.1e-15 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e10 | 5 | 3.3e-13 | 3.3e-13 | 1.8e-05 | 5.4e-15 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e12 | 5 | 3.2e-13 | 3.2e-13 | 1.4e-03 | 3.4e-15 | 0/0 | 0/0 |
| dirichlet_multinomial_lpmf | k:1e15 | 5 | 1.0e-12 | 1.0e-12 | 2.8e+00 | 4.4e-15 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e2 | 2 | 1.2e-14 | 1.6e-14 | 6.5e-13 | 1.7e-14 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e4 | 2 | 2.0e-11 | 1.4e-12 | 5.9e-11 | 1.5e-12 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e6 | 2 | 4.2e-09 | 2.1e-10 | 3.4e-09 | 1.6e-10 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e8 | 2 | 7.8e-07 | 2.0e-08 | 1.3e-06 | 1.6e-08 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e10 | 2 | 1.1e-04 | 1.0e-06 | 1.7e-04 | 2.1e-06 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e12 | 2 | 1.1e-02 | 2.1e-04 | 1.0e-02 | 2.0e-04 | 0/0 | 0/0 |
| lkj_corr_lpdf | K3:1e15 | 2 | 3.7e+00 | 1.5e-01 | 8.2e+00 | 1.4e-01 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e2 | 2 | 2.7e-14 | 1.3e-15 | 2.8e-13 | 1.7e-16 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e4 | 2 | 1.8e-11 | 6.1e-16 | 4.0e-11 | 7.2e-17 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e6 | 2 | 3.9e-09 | 2.0e-15 | 4.2e-10 | 1.4e-16 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e8 | 2 | 8.0e-07 | 3.6e-15 | 4.8e-07 | 7.2e-17 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e10 | 2 | 1.1e-04 | 3.7e-15 | 3.7e-05 | 1.5e-17 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e12 | 2 | 1.1e-02 | 1.8e-15 | 2.4e-02 | 2.0e-16 | 0/0 | 0/0 |
| lkj_corr_lpdf | I:1e15 | 2 | 3.9e+00 | 2.9e-15 | 1.5e+00 | 1.4e-16 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e2 | 2 | 2.4e-14 | 4.0e-15 | 6.7e-13 | 1.2e-15 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e4 | 2 | 1.9e-11 | 7.7e-13 | 5.7e-11 | 1.2e-13 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e6 | 2 | 3.9e-09 | 5.9e-11 | 3.4e-09 | 7.5e-12 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e8 | 2 | 8.0e-07 | 7.3e-09 | 1.3e-06 | 7.8e-10 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e10 | 2 | 1.1e-04 | 2.0e-07 | 1.7e-04 | 4.8e-07 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e12 | 2 | 1.1e-02 | 4.2e-05 | 1.1e-02 | 3.6e-05 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | K3:1e15 | 2 | 3.9e+00 | 2.4e-02 | 8.4e+00 | 2.4e-02 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e2 | 2 | 2.7e-14 | 1.3e-15 | 2.8e-13 | 1.7e-16 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e4 | 2 | 1.8e-11 | 6.1e-16 | 4.0e-11 | 7.2e-17 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e6 | 2 | 3.9e-09 | 2.0e-15 | 4.2e-10 | 1.4e-16 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e8 | 2 | 8.0e-07 | 3.6e-15 | 4.8e-07 | 7.2e-17 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e10 | 2 | 1.1e-04 | 3.7e-15 | 3.7e-05 | 1.5e-17 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e12 | 2 | 1.1e-02 | 1.8e-15 | 2.4e-02 | 2.0e-16 | 0/0 | 0/0 |
| lkj_corr_cholesky_lpdf | I:1e15 | 2 | 3.9e+00 | 2.9e-15 | 1.5e+00 | 1.4e-16 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e2 | 4 | 7.9e-14 | 4.4e-15 | 1.3e-13 | 1.5e-15 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e4 | 4 | 1.4e-11 | 7.8e-14 | 9.0e-12 | 2.7e-15 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e6 | 4 | 2.0e-09 | 7.8e-14 | 1.0e-09 | 6.1e-15 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e8 | 4 | 3.6e-07 | 5.1e-14 | 2.6e-07 | 8.2e-15 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e10 | 4 | 2.8e-05 | 1.7e-13 | 1.8e-05 | 3.5e-15 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e12 | 4 | 2.1e-03 | 1.6e-14 | 2.2e-03 | 6.3e-16 | 0/0 | 0/0 |
| yule_simon_lccdf | a:1e15 | 4 | 7.6e+00 | 2.0e-14 | 2.1e+00 | 3.5e-15 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e2 | 4 | 4.5e-16 | 2.4e-17 | 4.5e-16 | 7.5e-18 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e4 | 4 | 3.2e-16 | 5.5e-17 | 5.5e-16 | 9.4e-20 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e6 | 4 | 2.8e-16 | 5.0e-17 | 1.1e-15 | 3.9e-22 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e8 | 4 | 2.4e-15 | 3.9e-17 | 5.1e-15 | 8.9e-24 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e10 | 4 | 2.8e-15 | 8.3e-18 | 4.6e-15 | 1.4e-25 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e12 | 4 | 2.1e-15 | 2.2e-17 | 4.0e-15 | 8.0e-28 | 0/0 | 0/0 |
| yule_simon_cdf | a:1e15 | 4 | 7.8e-16 | 8.0e-19 | 1.0e-15 | 2.3e-30 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e2 | 4 | 4.3e-16 | 7.4e-18 | 4.5e-16 | 7.4e-18 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e4 | 4 | 3.5e-16 | 9.1e-20 | 5.5e-16 | 9.1e-20 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e6 | 4 | 2.4e-16 | 2.3e-22 | 1.1e-15 | 4.6e-22 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e8 | 4 | 2.5e-15 | 9.0e-24 | 5.1e-15 | 9.2e-24 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e10 | 4 | 2.8e-15 | 1.4e-25 | 4.6e-15 | 1.4e-25 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e12 | 4 | 2.1e-15 | 8.5e-28 | 4.0e-15 | 9.0e-28 | 0/0 | 0/0 |
| yule_simon_lcdf | a:1e15 | 4 | 7.7e-16 | 2.4e-30 | 1.0e-15 | 2.3e-30 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e2 | 4 | 3.4e-08 | 3.4e-08 | 8.7e-07 | 8.7e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e2 | 5 | 4.1e-08 | 4.1e-08 | 1.1e-06 | 1.1e-06 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e4 | 4 | 2.3e-08 | 2.3e-08 | 7.0e-07 | 7.0e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e4 | 5 | 2.9e-08 | 2.9e-08 | 8.1e-07 | 8.1e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e6 | 4 | 2.2e-08 | 2.1e-08 | 6.4e-07 | 6.4e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e6 | 5 | 2.7e-08 | 2.7e-08 | 7.7e-07 | 7.7e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e8 | 4 | 3.4e-08 | 2.1e-08 | 7.6e-07 | 6.4e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e8 | 5 | 5.0e-08 | 2.7e-08 | 1.0e-06 | 7.7e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e10 | 4 | 2.8e-06 | 2.1e-08 | 5.0e-05 | 6.4e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e10 | 5 | 4.9e-06 | 2.7e-08 | 4.2e-05 | 7.7e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e12 | 4 | 7.6e-05 | 2.1e-08 | 2.0e-03 | 6.4e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e12 | 5 | 5.0e-04 | 2.7e-08 | 6.1e-03 | 7.7e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e15 | 4 | 1.4e-01 | 2.1e-08 | 3.7e+00 | 6.4e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e15 | 5 | 2.1e-01 | 2.7e-08 | 1.9e+01 | 7.7e-07 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e0 | 4 | - | - | - | - | 4/0 | 4/0 |
| beta_binomial_cdf | ab3:1e0 | 5 | 6.5e-15 | 6.5e-15 | 2.1e-14 | 2.1e-14 | 0/0 | 0/0 |
| beta_binomial_cdf | tied:1e1 | 4 | 2.2e-09 | 2.2e-09 | 4.1e-08 | 4.1e-08 | 0/0 | 0/0 |
| beta_binomial_cdf | ab3:1e1 | 5 | 3.3e-08 | 3.3e-08 | 6.7e-07 | 6.7e-07 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e2 | 4 | 2.8e+01 | 7.6e-08 | 4.1e+01 | 1.7e-06 | 0/0 | 0/0 |
| beta_binomial_lcdf | ab3:1e2 | 5 | 1.0e-01 | 7.7e-08 | 2.9e+00 | 2.3e-06 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e4 | 4 | inf | 5.5e-08 | 7.0e+01 | 1.4e-06 | 0/2 | 0/0 |
| beta_binomial_lcdf | ab3:1e4 | 5 | 7.7e+00 | 5.4e-08 | 2.6e+01 | 1.7e-06 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e6 | 4 | inf | 5.0e-08 | 8.2e+01 | 1.3e-06 | 0/1 | 0/0 |
| beta_binomial_lcdf | ab3:1e6 | 5 | 1.1e+01 | 5.2e-08 | 3.3e+01 | 1.6e-06 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e8 | 4 | inf | 5.0e-08 | 1.2e+02 | 1.3e-06 | 0/2 | 0/0 |
| beta_binomial_lcdf | ab3:1e8 | 5 | inf | 5.2e-08 | 4.4e+01 | 1.6e-06 | 0/1 | 0/0 |
| beta_binomial_lcdf | tied:1e10 | 4 | inf | 5.0e-08 | 1.0e+02 | 1.3e-06 | 0/3 | 0/0 |
| beta_binomial_lcdf | ab3:1e10 | 5 | 2.1e+01 | 5.2e-08 | 3.8e+01 | 1.6e-06 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e12 | 4 | inf | 5.0e-08 | 7.4e+01 | 1.3e-06 | 0/3 | 0/0 |
| beta_binomial_lcdf | ab3:1e12 | 5 | 2.6e+01 | 5.2e-08 | 4.1e+01 | 1.6e-06 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e15 | 4 | inf | 5.0e-08 | 8.0e+01 | 1.3e-06 | 0/1 | 0/0 |
| beta_binomial_lcdf | ab3:1e15 | 5 | inf | 5.2e-08 | 4.1e+04 | 1.6e-06 | 0/3 | 0/0 |
| beta_binomial_lcdf | tied:1e0 | 4 | - | - | - | - | 4/0 | 4/0 |
| beta_binomial_lcdf | ab3:1e0 | 5 | 2.6e-13 | 2.6e-13 | 2.2e-13 | 2.2e-13 | 0/0 | 0/0 |
| beta_binomial_lcdf | tied:1e1 | 4 | 7.8e-05 | 4.7e-09 | 1.3e-04 | 7.9e-08 | 0/0 | 0/0 |
| beta_binomial_lcdf | ab3:1e1 | 5 | 3.7e-06 | 6.2e-08 | 1.3e-04 | 1.4e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e2 | 4 | 6.0e-08 | 6.0e-08 | 1.4e-06 | 1.4e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e2 | 5 | 8.6e-08 | 8.6e-08 | 2.0e-06 | 2.0e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e4 | 4 | 4.1e-08 | 4.1e-08 | 1.1e-06 | 1.1e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e4 | 5 | 1.8e-07 | 1.8e-07 | 2.1e-06 | 2.1e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e6 | 4 | 3.8e-08 | 3.7e-08 | 9.8e-07 | 9.8e-07 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e6 | 5 | 1.7e-07 | 1.7e-07 | 2.0e-06 | 2.0e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e8 | 4 | 6.0e-08 | 3.7e-08 | 1.1e-06 | 9.8e-07 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e8 | 5 | 1.8e-07 | 1.7e-07 | 2.1e-06 | 2.0e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e10 | 4 | 2.8e-06 | 3.7e-08 | 5.0e-05 | 9.8e-07 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e10 | 5 | 4.9e-06 | 1.7e-07 | 9.3e-05 | 2.0e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e12 | 4 | 7.6e-05 | 3.7e-08 | 2.0e-03 | 9.8e-07 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e12 | 5 | 5.0e-04 | 1.7e-07 | 6.1e-03 | 2.0e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e15 | 4 | 1.5e-01 | 3.7e-08 | 4.3e+00 | 9.8e-07 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e15 | 5 | 3.7e-01 | 1.7e-07 | 2.2e+01 | 2.0e-06 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e0 | 4 | - | - | - | - | 4/0 | 4/0 |
| beta_binomial_lccdf | ab3:1e0 | 5 | 6.7e-15 | 6.7e-15 | 2.2e-14 | 2.2e-14 | 0/0 | 0/0 |
| beta_binomial_lccdf | tied:1e1 | 4 | 4.2e-09 | 4.2e-09 | 7.2e-08 | 7.2e-08 | 0/0 | 0/0 |
| beta_binomial_lccdf | ab3:1e1 | 5 | 1.3e-07 | 1.3e-07 | 1.9e-06 | 1.9e-06 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | all:1e2 | 2 | 1.4e-13 | 7.1e-15 | 2.2e-13 | 2.7e-13 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | r:1e2 | 2 | 7.8e-14 | 2.6e-15 | 1.7e-06 | 1.7e-06 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e2 | 4 | 6.0e-17 | 6.0e-17 | 9.5e-11 | 9.5e-11 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | b:1e2 | 2 | 3.0e-15 | 3.6e-15 | 7.5e-07 | 7.5e-07 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e2 | 5 | 4.2e-14 | 2.6e-15 | 5.5e-08 | 5.5e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | all:1e4 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_cdf | r:1e4 | 2 | 5.8e-10 | 5.7e-10 | 6.7e-08 | 6.7e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e4 | 4 | 2.6e-17 | 2.6e-17 | 1.5e-14 | 6.0e-16 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | b:1e4 | 2 | 2.0e-09 | 2.0e-09 | 4.0e-08 | 4.0e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e4 | 5 | 8.7e-12 | 3.7e-15 | 3.2e-08 | 3.2e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | all:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_cdf | r:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_cdf | a:1e6 | 4 | 3.2e-17 | 3.2e-17 | 1.5e-15 | 2.0e-20 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | b:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_cdf | ab:1e6 | 5 | 5.1e-10 | 1.3e-14 | 3.2e-08 | 3.1e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e8 | 4 | 4.7e-17 | 4.7e-17 | 4.7e-15 | 3.9e-23 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e8 | 5 | 2.1e-08 | 2.6e-14 | 6.0e-07 | 3.1e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e10 | 4 | 5.0e-17 | 5.0e-17 | 1.3e-14 | 5.4e-25 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e10 | 5 | 1.5e-06 | 6.1e-15 | 7.0e-05 | 3.1e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e12 | 4 | 2.2e-17 | 2.2e-17 | 2.3e-14 | 1.0e-26 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e12 | 5 | 2.7e-04 | 7.8e-15 | 4.9e-03 | 3.1e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e15 | 4 | 4.8e-18 | 4.8e-18 | 1.9e-14 | 8.2e-30 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e15 | 5 | 3.3e-01 | 2.8e-14 | 1.5e+01 | 3.1e-08 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | all:1e1 | 4 | 7.7e-15 | 1.6e-15 | 3.8e-07 | 3.8e-07 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | r:1e1 | 4 | 1.4e-15 | 1.1e-16 | 2.0e-06 | 2.0e-06 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | a:1e1 | 4 | 5.9e-16 | 1.5e-16 | 1.6e-07 | 1.6e-07 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | b:1e1 | 4 | 1.8e-15 | 1.4e-15 | 1.9e-06 | 1.9e-06 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | ab:1e1 | 4 | 4.5e-15 | 3.4e-15 | 4.0e-07 | 4.0e-07 | 0/0 | 0/0 |
| beta_neg_binomial_cdf | all:1e8 | 1 | - | - | - | - | 1/0 | 1/0 |
| beta_neg_binomial_cdf | r:1e8 | 1 | - | - | - | - | 1/0 | 1/0 |
| beta_neg_binomial_cdf | b:1e8 | 1 | - | - | - | - | 1/0 | 1/0 |
| beta_neg_binomial_lcdf | all:1e2 | 2 | 2.1e+01 | 2.0e+01 | 4.3e+01 | 6.9e+01 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | r:1e2 | 2 | 1.7e-12 | 5.1e-13 | 3.1e-04 | 3.1e-04 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e2 | 4 | 4.0e-17 | 1.2e-17 | 1.0e-10 | 1.0e-10 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | b:1e2 | 2 | 1.6e-11 | 1.6e-11 | 9.8e-04 | 9.8e-04 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e2 | 5 | 9.6e-12 | 1.3e-13 | 2.2e-06 | 2.2e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | all:1e4 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lcdf | r:1e4 | 2 | 1.9e-03 | 1.9e-03 | 8.3e-02 | 8.3e-02 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e4 | 4 | 7.6e-19 | 7.6e-19 | 1.5e-14 | 6.0e-16 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | b:1e4 | 2 | 2.3e+00 | 2.3e+00 | 1.2e+00 | 1.2e+00 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e4 | 5 | 2.1e-09 | 1.9e-13 | 1.5e-06 | 1.5e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | all:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lcdf | r:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lcdf | a:1e6 | 4 | 5.0e-21 | 5.0e-21 | 1.5e-15 | 1.5e-20 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | b:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lcdf | ab:1e6 | 5 | 8.0e-08 | 8.3e-13 | 1.8e-06 | 2.0e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e8 | 4 | 4.0e-22 | 2.4e-23 | 2.6e-14 | 1.2e-23 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e8 | 5 | 2.9e-06 | 3.7e-13 | 1.5e-04 | 2.0e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e10 | 4 | 1.6e-24 | 4.9e-25 | 1.3e-14 | 5.5e-25 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e10 | 5 | 1.4e-04 | 6.1e-13 | 1.4e-02 | 2.0e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e12 | 4 | 3.4e-26 | 1.2e-26 | 2.3e-14 | 1.2e-26 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e12 | 5 | 6.8e-02 | 9.4e-13 | 6.3e-01 | 2.0e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e15 | 4 | 1.0e-28 | 1.4e-29 | 4.4e-14 | 1.5e-29 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e15 | 5 | 6.0e-01 | 3.4e-13 | 2.2e+01 | 2.0e-06 | 2/0 | 0/0 |
| beta_neg_binomial_lcdf | all:1e1 | 4 | 1.4e-12 | 2.6e-13 | 2.4e-05 | 2.4e-05 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | r:1e1 | 4 | 1.1e-14 | 2.2e-16 | 1.6e-05 | 1.6e-05 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | a:1e1 | 4 | 9.0e-16 | 2.3e-16 | 2.7e-07 | 2.7e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | b:1e1 | 4 | 3.2e-14 | 2.9e-15 | 3.3e-05 | 3.3e-05 | 0/0 | 0/0 |
| beta_neg_binomial_lcdf | ab:1e1 | 4 | 1.5e-13 | 9.3e-14 | 2.0e-05 | 2.0e-05 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | all:1e2 | 2 | 1.4e-13 | 7.1e-15 | 2.7e-13 | 2.7e-13 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | r:1e2 | 2 | 8.2e-14 | 2.7e-15 | 1.8e-06 | 1.8e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e2 | 4 | 3.3e-15 | 1.2e-14 | 6.5e-08 | 6.5e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | b:1e2 | 2 | 3.0e-15 | 3.6e-15 | 7.5e-07 | 7.5e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e2 | 5 | 4.2e-14 | 2.0e-14 | 1.5e-07 | 1.5e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | all:1e4 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lccdf | r:1e4 | 2 | 5.8e-10 | 5.7e-10 | 6.7e-08 | 6.7e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e4 | 4 | 4.7e-15 | 6.1e-15 | 8.9e-11 | 5.2e-11 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | b:1e4 | 2 | 2.0e-09 | 2.0e-09 | 4.0e-08 | 4.0e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e4 | 5 | 9.3e-12 | 8.1e-15 | 1.1e-07 | 1.1e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | all:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lccdf | r:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lccdf | a:1e6 | 4 | 9.4e-15 | 9.4e-15 | 1.1e-09 | 7.6e-13 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | b:1e6 | 2 | - | - | - | - | 2/0 | 2/0 |
| beta_neg_binomial_lccdf | ab:1e6 | 5 | 6.5e-10 | 4.3e-14 | 1.0e-07 | 9.9e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e8 | 4 | 9.6e-14 | 1.8e-14 | 2.7e-07 | 1.1e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e8 | 5 | 2.7e-08 | 1.4e-13 | 7.4e-07 | 9.9e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e10 | 4 | 1.2e-13 | 1.2e-13 | 1.8e-05 | 7.4e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e10 | 5 | 3.6e-06 | 3.6e-14 | 7.2e-05 | 9.9e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e12 | 4 | 3.3e-14 | 4.3e-13 | 1.9e-03 | 5.6e-15 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e12 | 5 | 2.7e-04 | 2.4e-13 | 5.0e-03 | 9.9e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e15 | 4 | 1.1e-13 | 1.2e-13 | 8.2e+00 | 1.2e-14 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e15 | 5 | 2.9e-01 | 7.1e-14 | 1.3e+01 | 9.9e-08 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | all:1e1 | 4 | 9.0e-15 | 1.7e-15 | 8.7e-07 | 8.7e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | r:1e1 | 4 | 1.8e-14 | 2.3e-15 | 3.4e-06 | 3.4e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | a:1e1 | 4 | 1.7e-14 | 5.1e-15 | 8.2e-07 | 8.2e-07 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | b:1e1 | 4 | 4.6e-15 | 2.7e-15 | 3.4e-06 | 3.4e-06 | 0/0 | 0/0 |
| beta_neg_binomial_lccdf | ab:1e1 | 4 | 1.2e-14 | 4.7e-15 | 8.7e-07 | 8.7e-07 | 0/0 | 0/0 |
| student_t_lpdf | B8 | 230 | 1.6e-01 | 1.3e-11 | 1.4e-01 | 4.5e-12 | 0/0 | 0/0 |
| multi_student_t_cholesky_lpdf | B1 | 200 | 2.1e-15 | 2.1e-15 | 2.1e-16 | 1.8e-16 | 0/0 | 0/0 |
| multi_student_t_cholesky_lpdf | B2 | 200 | 3.7e-14 | 4.6e-15 | 1.5e-14 | 2.5e-15 | 0/0 | 0/0 |
| multi_student_t_cholesky_lpdf | B3 | 200 | 1.7e-07 | 6.4e-15 | 1.5e-07 | 4.3e-15 | 0/0 | 0/0 |
| multi_student_t_cholesky_lpdf | B4 | 200 | 3.1e-13 | 4.5e-15 | 4.5e-13 | 2.3e-15 | 0/0 | 0/0 |
| multi_student_t_cholesky_lpdf | B5 | 200 | 1.4e-07 | 8.1e-15 | 1.5e-07 | 4.5e-15 | 0/0 | 0/0 |
| multi_student_t_cholesky_lpdf | B6 | 13 | 1.8e-03 | 2.5e-15 | 8.9e-05 | 8.7e-17 | 0/0 | 0/0 |
| multi_student_t_lpdf | B1 | 200 | 1.9e-15 | 1.9e-15 | 3.7e-16 | 3.7e-16 | 0/0 | 0/0 |
| multi_student_t_lpdf | B2 | 200 | 2.3e-14 | 4.6e-15 | 2.3e-14 | 4.7e-15 | 0/0 | 0/0 |
| multi_student_t_lpdf | B3 | 200 | 2.0e-07 | 5.6e-15 | 1.4e-07 | 5.0e-15 | 0/0 | 0/0 |
| multi_student_t_lpdf | B4 | 200 | 2.1e-13 | 4.3e-15 | 5.2e-13 | 6.3e-15 | 0/0 | 0/0 |
| multi_student_t_lpdf | B5 | 200 | 1.8e-07 | 5.6e-15 | 2.2e-07 | 4.8e-15 | 0/0 | 0/0 |
| multi_student_t_lpdf | B6 | 13 | 1.6e-03 | 2.2e-15 | 1.7e-07 | 9.2e-16 | 0/0 | 0/0 |
| yule_simon_lpmf | B1 | 200 | 6.1e-15 | 4.3e-15 | 3.9e-16 | 3.9e-16 | 0/0 | 0/0 |
| yule_simon_lpmf | B2 | 200 | 3.1e-13 | 2.7e-14 | 1.1e-13 | 7.0e-15 | 0/0 | 0/0 |
| yule_simon_lpmf | B3 | 200 | 2.4e-07 | 3.3e-14 | 2.7e-07 | 1.9e-15 | 0/0 | 0/0 |
| yule_simon_lpmf | B4 | 200 | 3.4e-09 | 1.5e-12 | 1.5e-12 | 9.2e-13 | 0/0 | 0/0 |
| yule_simon_lpmf | B5 | 200 | 2.2e-07 | 6.9e-10 | 2.5e-07 | 1.7e-10 | 0/0 | 0/0 |
| yule_simon_lpmf | B6 | 13 | 3.9e-03 | 2.7e-15 | 3.6e-03 | 1.5e-16 | 0/0 | 0/0 |
| lbeta (var), d/dx | B7 | 275 | - | - | 6.3e-01 | 1.0e-14 | 0/0 | 0/0 |
