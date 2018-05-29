If you use our work in your research, please cite as:

B. Dumitrescu and P. Irofti, Dictionary Learning Algorithms and Applications, Springer, 2018
```
@book{DL_book,
author    = {Dumitrescu, B. and Irofti, P.},
title     = {Dictionary Learning Algorithms and Applications},
year      = {2018},
publisher = {Springer},
}
```

Run setupDL.m before executing any of the scripts described below.

================================================================================
Chapter 2

  Examples and figures:
    tab_2_1_overlap:  Denoising with overlapping patches

================================================================================
Chapter 3

  Basic functions:
    run_denoising:    Run denoising tests with varied images, noises and
                      algorithms. SSIM and PSNR results displayed on screen.
		      Data set saved in the data directory.
		      Needs plenty of disk space and execution time.

    run_repl:         Run replacement tests for fig_3_3
    run_replc:        Run replacement tests for fig_3_4

    run_DL:           Run in-depth DL tests on images for 500 iterations.
                      Data set used by most figures from this chapter.
                      Needs plenty of disk space and execution time.

    run_reg_fig:      Run regularized DL tests on synthetic data.
    run_reg_table:    Run regularized DL tests on images.
    run_reg_table_itN:Extend the above tests with regularized SimCO results
                      executed with multiple internal iterations.

    run_recov:        Run dictionary recovery tests on synthetic data.
		      Needs disk space and execution time.

  Examples and figures:
    fig_3_2_init:     Initialization with random vectors and signals
    fig_3_3_repl:     AK-SVD with several atom replacement strategies
    fig_3_4_replc:    AK-SVD with several atom replacement strategies (centered)
    fig_3_5_6_iter:   RMSE evolution over 50 and 500 iterations
    fig_3_7_8_sn:     Final RMSE with various dictionary sizes and sparsity
    tab_3_1_lena:     PSNR and SSIM when denoising Lena (see run_denoising)
    fig_3_9_recov:    Evolution of atom recovery percentage (see run_recov)
    fig_3_10_denoise: AK-SVD denoising for various images and noises
    fig_3_11_laols:   RMSE evolution for DL with OMP and LAOLS

================================================================================
Chapter 4

  Basic functions:
    gen_rand_untf:  generate random UNTF, Algorithm 4.1

  Examples and figures:
    fig_4_1_reg: RMSE for standard and regularized algorithms on synth data
    tab_4_1_reg: RMSE for standard and regularized algorithms on images
    ex_4_12_coh: distribution of UNTF atom products
    ex_4_13_coh: distribution of atom products for learned dictionaries
    ex_4_16_coh: bound coherence with IPR

================================================================================
Chapter 5
  Basic functions:
    apru:        DL with a 1-norm penalty, Algorithm 5.1
    dl_select:   dictionary selection, Algorithm 5.2
    online_dl:   online DL, Algorithm 5.3
    rls_dl:      RLS DL, Algorithm 5.4

  Examples and figures:
    ex_5_2_apru:    APrU vs AK-SVD
    ex_5_9_select:  dictionary selection vs AK-SVD
    ex_5_12_online: comparison of online DL algorithms
    ex_5_13_miss:   AK-SVD with incomplete data

================================================================================
Chapter 6

  Examples and figures:
    ex_6_1_size:    RMSE for different dictionary sizes and sparsity levels
    ex_6_2:         ??? Stagewise K-SVD

================================================================================
Chapter 8

Classification:
    clas_src_dl:     Design SRC classifier with dictionaries trained for each class, section 8.3
    clas_discrim_dl: Design discriminative DL classifier, section 8.5.2
    fig_8_6_class:   Discriminative and LC-KSVD succesful classification
    fig_8_7_coeff:   Sparse representation coefficients used for classification

================================================================================
Chapter 9
  Basic functions:
    omp_ker:         kernel OMP, batch version equivalent to Algorithm 9.1
    aksvd_ker:       kernel AK-SVD Algorithm 9.2
    aksvd_ker_discr: AK-SVD for kernel discriminative DL, section 9.6.3

  Classification:
    clas_src_dl_ker: SRC classifier with dictionaries trained for each class,
                     using kernel AK-SVD
 
  Examples and figures:
    ex_sec_9_6_kerdl: classification with kernel methods for a simple
                      two-class problem
    ex_sec_9_6_2_nystrom: Nystrom method with simple two-class problem
