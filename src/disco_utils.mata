********************************************************************************
*** BEGIN MATA CODE
********************************************************************************

mata:

// -----------------------------------------------------------------------------
// Output structure for confidence intervals
struct CI_out {
    real matrix qdiff_lower,    // Lower bound for quantile differences
              qdiff_upper,      // Upper bound for quantile differences
              cdiff_lower,      // Lower bound for CDF differences
              cdiff_upper       // Upper bound for CDF differences
}

// Final returned structure for disco
struct disco_out {
    real matrix weights, quantile_diff, cdf_diff, quantile_synth, cdf_synth,
                quantile_t, cdf_t, cids
}

// -----------------------------------------------------------------------------
// Returns unique values from a vector
real vector get_unique(real vector x)
{
    real scalar N
    real vector sortedx, y

    sortedx = sort(x,1)
    y       = uniqrows(sortedx)
    return(y)
}

// -----------------------------------------------------------------------------
// Evenly spaced probability grid on [q_min, q_max]. Built with the same
// scalar recursion as the original code so the grid values are bit-identical.
real vector disco_prob_grid(real scalar Npoints, real scalar q_min,
                            real scalar q_max)
{
    real scalar j2
    real vector p

    p = J(Npoints,1,.)
    for (j2=1; j2<=Npoints; j2++) {
        p[j2,1] = q_min + (q_max - q_min)*(j2-1)/(Npoints-1)
    }
    return(p)
}

// -----------------------------------------------------------------------------
// Type=7 quantiles of PRE-SORTED data at probabilities p, vectorized. The
// clamping reproduces the p<=0, p>=1, and boundary branches of the original
// per-point loop; missing p follows the p>=1 branch.
real vector disco_quantile_sorted(real vector Xs, real vector p)
{
    real scalar N
    real vector pp, ones, alpha, fl, gam

    N = rows(Xs)
    if (N==0) return(J(rows(p),1,.))
    if (N==1) return(J(rows(p),1,Xs[1]))

    pp    = editmissing(p, 1)
    ones  = J(rows(pp),1,1)
    alpha = (N-1):*pp :+ 1
    fl    = floor(alpha)
    fl    = rowmax((fl, ones))
    fl    = rowmin((fl, (N-1):*ones))
    gam   = alpha :- fl
    gam   = rowmax((gam, 0:*ones))
    gam   = rowmin((gam, ones))
    return(Xs[fl]:*(1:-gam) :+ Xs[fl:+1]:*gam)
}

// -----------------------------------------------------------------------------
// Compute quantiles at arbitrary probabilities p using type=7 interpolation
real vector disco_quantile_points(real vector X, real vector p)
{
    return(disco_quantile_sorted((cols(X)>1 ? sort(X',1) : sort(X,1)), p))
}

// -----------------------------------------------------------------------------
// Compute quantiles at a user-chosen number of points (Npoints)
real vector disco_quantile(real vector X, real scalar Npoints,
                           real scalar q_min, real scalar q_max)
{
    return(disco_quantile_points(X, disco_prob_grid(Npoints, q_min, q_max)))
}

// -----------------------------------------------------------------------------
// Empirical CDF of PRE-SORTED data at grid points: #(X <= g)/N, exact under
// ties. Binary search per grid point when the grid is small relative to the
// data; otherwise one merge pass, where data rows (flag 0) sort before equal
// grid rows (flag 1) so the running data count at a grid row is #(X <= g).
real vector disco_cdf_sorted(real vector Xs, real vector grid)
{
    real scalar N, G, g5, lo, hi, mid
    real vector gc, ord, isX, cum, sel, out
    real matrix Z

    gc = (cols(grid)>1 ? grid' : grid)
    N  = rows(Xs)
    G  = rows(gc)

    if (8*G <= N) {
        out = J(G,1,0)
        for (g5=1; g5<=G; g5++) {
            lo = 0
            hi = N + 1
            while (hi - lo > 1) {
                mid = floor((lo+hi)/2)
                if (Xs[mid] <= gc[g5]) lo = mid
                else                   hi = mid
            }
            out[g5] = lo/N
        }
        return(out)
    }

    Z   = (Xs, J(N,1,0)) \ (gc, J(G,1,1))
    ord = order(Z, (1,2))
    isX = (ord :<= N)
    cum = runningsum(isX)
    sel = selectindex(1 :- isX)
    out = J(G,1,.)
    out[ord[sel] :- N] = cum[sel] :/ N
    return(out)
}

// -----------------------------------------------------------------------------
// Compute CDF values at specified grid points
real vector cdf_builder(real vector x, real vector grid)
{
    return(disco_cdf_sorted((cols(x)>1 ? sort(x',1) : sort(x,1)), grid))
}

// -----------------------------------------------------------------------------
// Helper to compute CDF at given grid points
real vector cdf_at_points(real vector x, real vector grid)
{
    return(cdf_builder(x, grid))
}

// -----------------------------------------------------------------------------
// Sorted with-replacement resample read directly off a PRE-SORTED cell:
// exponential spacings generate the order statistics of n iid U(0,1) draws,
// so Xs[ceil(n*U_(i))] is the sorted version of a bootstrap resample without
// an O(n log n) sort per draw.
real vector disco_boot_sorted(real vector Xs)
{
    real scalar n
    real vector S

    n = rows(Xs)
    if (n==0) return(Xs)
    S = runningsum(-ln(runiform(n+1,1)))
    return(Xs[ceil((S[|1 \ n|] :/ S[n+1]) :* n)])
}

// -----------------------------------------------------------------------------
// Split the panel once into per-(period, unit) outcome vectors, sorted
// ascending. Column 1 is the target, column 1+j is cids[j]. The estimator,
// bootstrap, and permutation test only touch these cells afterwards, which
// avoids rescanning the full dataset inside hot loops.
pointer(real colvector) matrix disco_build_cells(real vector y, real vector id,
    real vector tt, real scalar target_id, real scalar T_max, real vector cids)
{
    real scalar t6, j6, Jn
    real vector yt, idt, cd6
    pointer(real colvector) matrix P

    Jn = length(cids)
    P  = J(T_max, Jn+1, NULL)
    for (t6=1; t6<=T_max; t6++) {
        yt  = select(y,  tt:==t6)
        idt = select(id, tt:==t6)
        cd6 = select(yt, idt:==target_id)
        P[t6,1] = (rows(cd6)>0 ? &(sort(cd6,1)) : &J(0,1,.))
        for (j6=1; j6<=Jn; j6++) {
            cd6 = select(yt, idt:==cids[j6])
            P[t6,1+j6] = (rows(cd6)>0 ? &(sort(cd6,1)) : &J(0,1,.))
        }
    }
    return(P)
}

// -----------------------------------------------------------------------------
// Quantile/CDF stores for all units and periods, computed once from the
// sorted cells. QG and CG stack G rows per period for every period; WM holds
// the M-grid weighting objects for the pre-treatment periods only: quantile
// functions when mixture==0, CDFs when mixture==1. M-grid objects are never
// needed for post-treatment periods or for the other mode.
void disco_precompute(pointer(real colvector) matrix P, real scalar T0,
    real scalar T_max, real scalar M, real scalar G,
    real scalar q_min, real scalar q_max, real scalar mixture,
    real vector gridG, real vector gridM,
    real matrix QG, real matrix CG, real matrix WM)
{
    real scalar nu, t7, u7, npre
    real vector pG, pM, cell

    nu   = cols(P)
    npre = max((T0-1, 0))
    pG   = disco_prob_grid(G, q_min, q_max)
    pM   = disco_prob_grid(M, q_min, q_max)

    QG = J(G*T_max, nu, .)
    CG = J(G*T_max, nu, .)
    WM = J(M*npre,  nu, .)

    for (t7=1; t7<=T_max; t7++) {
        for (u7=1; u7<=nu; u7++) {
            cell = *(P[t7,u7])
            QG[|(t7-1)*G+1, u7 \ t7*G, u7|] = disco_quantile_sorted(cell, pG)
            CG[|(t7-1)*G+1, u7 \ t7*G, u7|] = disco_cdf_sorted(cell, gridG)
            if (t7<=T0-1) {
                WM[|(t7-1)*M+1, u7 \ t7*M, u7|] = (mixture==1 ?
                    disco_cdf_sorted(cell, gridM) : disco_quantile_sorted(cell, pM))
            }
        }
    }
}

// -----------------------------------------------------------------------------
// 1. Quadratic-Programming Weights from PRE-COMPUTED Quantiles
real vector disco_solve_weights(real matrix control_quantiles, real vector target_quantiles, real scalar simplex)
{
    real scalar J
    real matrix Gm, CE, CI, C_in
    real vector g0, ce0, ci0, w, res

    J    = cols(control_quantiles)
    C_in = control_quantiles

    Gm = 2 :* (C_in' * C_in)
    g0 = -2 :* (C_in' * target_quantiles)

    CE   = J(J,1,1)
    ce0  = -1

    if (simplex==1) {
        CI  = I(J)
        ci0 = J(J,1,0)
    }
    else {
        CI  = J(J,1,1)
        ci0 = -1e20
    }

    res = solve_quadprog(Gm, g0, CE, ce0, CI, ci0)
    w   = res[1..J]

    return(w)
}

//------------------------------------------
// Mixture (CDF-based) weights via linear programming
//------------------------------------------
real vector disco_mixture_weights(real matrix control_cdf, real vector target_cdf, real scalar simplex)
{
    // --- Declarations ---
    real scalar J, Grows, total_vars, val, ss
    real matrix c, ecmat, lowerbd, upperbd
    real vector bec, sol, w, tcdf
    class LinearProgram scalar q

    // "Grows" = number of grid points, same as rows in control_cdf
    // J       = number of controls (columns)
    Grows      = rows(control_cdf)
    J          = cols(control_cdf)
    total_vars = J + 2*Grows
    tcdf       = (cols(target_cdf)>1 ? target_cdf' : target_cdf)

    // objective: minimize the sum of the 2*Grows slack variables
    c = (J(1, J, 0), J(1, 2*Grows, 1))

    // equality constraints: sum of weights = 1 (top row), then for each grid
    // point sum_j w_j F_j(y_g) - s+_g + s-_g = F_target(y_g)
    ecmat = (J(1, J, 1), J(1, 2*Grows, 0)) \ (control_cdf, -I(Grows), I(Grows))
    bec   = 1 \ tcdf

    // Bounds
    lowerbd = J(1, total_vars, .)
    upperbd = J(1, total_vars, .)

    // If simplex=1, force weights >= 0
    if (simplex == 1) {
        lowerbd[|1, 1 \ 1, J|] = J(1, J, 0)
    }
    // Slacks must be >= 0
    lowerbd[|1, J+1 \ 1, total_vars|] = J(1, 2*Grows, 0)

    // Solve linear program
    q = LinearProgram()
    q.setCoefficients(c)
    q.setMaxOrMin("min")
    q.setEquality(ecmat, bec)
    q.setBounds(lowerbd, upperbd)

    val = q.optimize()
    if (q.errorcode() != 0) {
        // fallback uniform
        w = J(J, 1, 1/J)
        return(w)
    }
    sol = q.parameters()
    w   = sol[1..J]' // first J parameters = the weights

    // Enforce sum(w)=1 if minor numeric drift
    ss = sum(w)
    if (abs(ss - 1) > 1e-8) w = w / ss

    return(w)
}

// -----------------------------------------------------------------------------
// Run the DiSCo estimator from the precomputed stores, treating column tc as
// the target and the remaining columns as controls. Pre-treatment periods use
// their own weights and post-treatment periods the averaged weights, exactly
// like the original two-loop implementation.
struct disco_out scalar disco_run_from_store(real matrix QG, real matrix CG,
    real matrix WM, real scalar tc, real scalar T0, real scalar T_max,
    real scalar M, real scalar G, real vector gridG,
    real scalar simplex, real scalar mixture)
{
    struct disco_out scalar r
    real scalar nu, Jn, t8, m8, p8, gg8
    real vector idx, w, W_avg, Tq, Tc, Q_synth, C_synth, wfinal
    real matrix weights_store, period_weights, WMt, CqG, CcG
    real matrix quantile_diff, cdf_diff, Q_synth_all, C_synth_all
    real matrix Q_target_all, C_target_all

    nu  = cols(QG)
    Jn  = nu - 1
    idx = select((1::nu), (1::nu):!=tc)

    //-------------------------------------------------------
    // 1. Weights for each pre-treatment period
    //-------------------------------------------------------
    weights_store  = J(max((T0-1,0)), Jn, .)
    period_weights = J(T_max, Jn, .)

    for (t8=1; t8<=T0-1; t8++) {
        WMt = WM[|(t8-1)*M+1, 1 \ t8*M, nu|]
        if (mixture==0) {
            w = disco_solve_weights(WMt[., idx], WMt[., tc], simplex)
        }
        else {
            w = disco_mixture_weights(WMt[., idx], WMt[., tc], simplex)'
        }
        weights_store[t8, .]  = w
        period_weights[t8, .] = w
    }

    //-------------------------------------------------------
    // 2. Average weights for the post period
    //-------------------------------------------------------
    if (T0-1 < 1) {
        W_avg = J(Jn,1,1/Jn)
    }
    else {
        W_avg = (colsum(weights_store) / (T0-1))'
    }

    for (t8=T0; t8<=T_max; t8++) {
        period_weights[t8, .] = W_avg'
    }

    //-------------------------------------------------------
    // 3. Build final synthetic distributions per period
    //-------------------------------------------------------
    quantile_diff = J(G, T_max, .)
    cdf_diff      = J(G, T_max, .)
    Q_target_all  = J(G, T_max, .)
    C_target_all  = J(G, T_max, .)
    Q_synth_all   = J(G, T_max, .)
    C_synth_all   = J(G, T_max, .)

    for (t8=1; t8<=T_max; t8++) {
        Tq  = QG[|(t8-1)*G+1, tc \ t8*G, tc|]
        Tc  = CG[|(t8-1)*G+1, tc \ t8*G, tc|]
        CqG = QG[|(t8-1)*G+1, 1 \ t8*G, nu|][., idx]
        CcG = CG[|(t8-1)*G+1, 1 \ t8*G, nu|][., idx]

        wfinal = period_weights[t8, .]'

        if (mixture == 0) {
            Q_synth = CqG * wfinal
            // empirical CDF of the G synthetic quantiles = counts/G
            C_synth = disco_cdf_sorted(sort(Q_synth,1), gridG)
        }
        else {
            C_synth = CcG * wfinal
            Q_synth = J(G,1,.)
            for (m8=1; m8<=G; m8++) {
                p8  = (m8-1)/(G-1)
                gg8 = 1
                while (gg8 < length(C_synth) & C_synth[gg8] < p8) gg8++
                if (gg8>G) gg8 = G
                Q_synth[m8] = gridG[gg8]
            }
        }

        Q_target_all[, t8] = Tq
        C_target_all[, t8] = Tc
        Q_synth_all[, t8]  = Q_synth
        C_synth_all[, t8]  = C_synth

        quantile_diff[, t8] = Tq - Q_synth
        cdf_diff[, t8]      = Tc - C_synth
    }

    r.weights         = W_avg
    r.quantile_diff   = quantile_diff
    r.cdf_diff        = cdf_diff
    r.quantile_synth  = Q_synth_all
    r.cdf_synth       = C_synth_all
    r.quantile_t      = Q_target_all
    r.cdf_t           = C_target_all
    r.cids            = J(0,0,.)
    return(r)
}

// -----------------------------------------------------------------------------
// Main DISCO function with both M & G
struct disco_out disco_full_run(real vector y, real vector id, real vector tt,
                              real scalar target_id, real scalar T0, real scalar T_max,
                              real scalar M, real scalar G,
                              real scalar q_min, real scalar q_max,
                              real scalar simplex, real scalar mixture)
{
    struct disco_out scalar r
    real scalar amin, amax
    real vector uid, cids, gridG, gridM
    real matrix QG, CG, WM
    pointer(real colvector) matrix P

    uid  = get_unique(id)
    cids = select(uid, uid:!=target_id)

    amin = min(y)
    amax = max(y)

    st_numscalar("amin", amin)
    st_numscalar("amax", amax)

    gridG = range(amin, amax, (amax - amin)/(G-1))
    gridM = range(amin, amax, (amax - amin)/(M-1))

    P = disco_build_cells(y, id, tt, target_id, T_max, cids)

    QG = .
    CG = .
    WM = .
    disco_precompute(P, T0, T_max, M, G, q_min, q_max, mixture, gridG, gridM,
                     QG, CG, WM)

    r = disco_run_from_store(QG, CG, WM, 1, T0, T_max, M, G, gridG,
                             simplex, mixture)
    r.cids = cids
    return(r)
}

// -----------------------------------------------------------------------------
// Pre/post RMSE ratio of the quantile treatment effects
real scalar disco_ratio_from_diff(real matrix quantile_diff,
                                  real scalar T0, real scalar T_max)
{
    real scalar pre_dist, pre_count, post_dist, post_count, dist_t, t_loop

    pre_dist   = 0
    pre_count  = 0
    post_dist  = 0
    post_count = 0

    for (t_loop=1; t_loop<=T_max; t_loop++) {
        dist_t = mean((quantile_diff[,t_loop]:^2))
        if (t_loop<T0) {
            pre_dist   = pre_dist + dist_t
            pre_count  = pre_count + 1
        }
        else {
            post_dist  = post_dist + dist_t
            post_count = post_count + 1
        }
    }

    if (pre_count==0 | post_count==0) return(.)

    return(sqrt((post_dist/post_count))/sqrt((pre_dist/pre_count)))
}

// -----------------------------------------------------------------------------
// Compute ratio, updated to pass M (and G) into disco_full_run
real scalar disco_compute_ratio(real vector y, real vector id, real vector tt,
                              real scalar target_id, real scalar T0, real scalar T_max,
                              real scalar M, real scalar G,
                              real scalar q_min, real scalar q_max,
                              real scalar simplex, real scalar mixture)
{
    struct disco_out scalar rr

    rr = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    return(disco_ratio_from_diff(rr.quantile_diff, T0, T_max))
}

// -----------------------------------------------------------------------------
// Permutation test. The per-unit quantile/CDF stores are computed once and
// reused for every placebo run: only the weight solves and the synthetic
// distributions are redone with each unit treated as the target in turn.
real scalar disco_permutation_test(real vector y, real vector id, real vector tt,
                                 real scalar target_id, real scalar T0, real scalar T_max,
                                 real scalar M, real scalar G,
                                 real scalar q_min, real scalar q_max,
                                 real scalar simplex, real scalar mixture)
{
    struct disco_out scalar rr
    real scalar amin, amax, actual_ratio, rj, J, count, j_loop
    real vector uid, cids, gridG, gridM
    real matrix QG, CG, WM
    pointer(real colvector) matrix P

    uid  = get_unique(id)
    cids = select(uid, uid:!=target_id)
    J    = length(cids)

    amin = min(y)
    amax = max(y)
    gridG = range(amin, amax, (amax - amin)/(G-1))
    gridM = range(amin, amax, (amax - amin)/(M-1))

    P = disco_build_cells(y, id, tt, target_id, T_max, cids)

    QG = .
    CG = .
    WM = .
    disco_precompute(P, T0, T_max, M, G, q_min, q_max, mixture, gridG, gridM,
                     QG, CG, WM)

    rr = disco_run_from_store(QG, CG, WM, 1, T0, T_max, M, G, gridG,
                              simplex, mixture)
    actual_ratio = disco_ratio_from_diff(rr.quantile_diff, T0, T_max)

    count = 0
    for (j_loop=1; j_loop<=J; j_loop++) {
        rr = disco_run_from_store(QG, CG, WM, 1+j_loop, T0, T_max, M, G, gridG,
                                  simplex, mixture)
        rj = disco_ratio_from_diff(rr.quantile_diff, T0, T_max)
        if (rj>=actual_ratio) count = count + 1
    }

    return((count+1)/(J+1))
}

//--------------------------------------------------------------
// disco_bootstrap_CI
//  * Resamples each (period, unit) cell, re-estimates the weights on the
//    resampled pre-treatment data, and rebuilds the synthetic distributions
//  * Data cells are pre-split and pre-sorted once; resamples are drawn
//    already sorted; only the M-grid objects needed for the weight step of
//    the requested mode are computed
//  * Builds pointwise or uniform intervals for both quantile & cdf
//--------------------------------------------------------------
struct CI_out scalar disco_bootstrap_CI(
    real vector y,
    real vector id,
    real vector tt,
    real scalar target_id,
    real scalar T0,
    real scalar T_max,
    real scalar M,
    real scalar G,
    real scalar q_min,
    real scalar q_max,
    real scalar simplex,
    real scalar mixture,
    real scalar boots,
    real scalar cl,
    real scalar uniform,
    real matrix quantile_diff,
    real matrix cdf_diff
)
{
    //----------------------------------------------------------
    // Declarations
    //----------------------------------------------------------
    struct CI_out scalar co

    real scalar amin, amax, Jn, nu
    real scalar b_loop, t_loop, j_loop, m8, p8, gg8
    real scalar t_loop2, idx_loop3, idx_loop4
    real scalar alpha, lower_idx, upper_idx
    real scalar t_loop3, index_q, m_i2, c_i2
    real scalar alpha_q, cval_q, base_val2, cval_c, base_val_c

    real vector uid, cids, gridG, gridM, pG, pM, vals, tmp
    real vector bs, target_w, w, W_avg, Q_synth, C_synth
    real matrix quantile_diff_boot, cdf_diff_boot
    real matrix qdiff_lower, qdiff_upper, cdiff_lower, cdiff_upper
    real matrix qdiff_mat, qdiff_err, sorted_diffs_q
    real matrix cdf_mat, cdf_err, sorted_diffs_c
    real vector qmax_abs, cmax_abs
    real matrix tQG, tCG, cQG, cCG, Wmat, weights_all
    pointer(real colvector) matrix P

    //----------------------------------------------------------
    // 1. Build grids and pre-split the data into sorted cells
    //----------------------------------------------------------
    amin = min(y)
    amax = max(y)

    // final G-based grid for the final distribution
    if (G>1 & amax>amin) {
        gridG = range(amin, amax, (amax - amin)/(G-1))
    } else {
        gridG = J(G,1,amin)
    }

    // M-based grid for weighting
    if (M>1 & amax>amin) {
        gridM = range(amin, amax, (amax - amin)/(M-1))
    } else {
        // corner case
        gridM = J(M,1,amin)
    }

    pG = disco_prob_grid(G, q_min, q_max)
    pM = disco_prob_grid(M, q_min, q_max)

    uid  = get_unique(id)
    cids = select(uid, uid:!=target_id)
    Jn   = length(cids)
    nu   = Jn + 1

    P = disco_build_cells(y, id, tt, target_id, T_max, cids)

    //----------------------------------------------------------
    // 2. Storage for bootstrap draws and per-replication work
    //----------------------------------------------------------
    quantile_diff_boot = J(G*T_max, boots, .)
    cdf_diff_boot      = J(G*T_max, boots, .)

    tQG         = J(G, T_max, .)
    tCG         = J(G, T_max, .)
    cQG         = J(G*T_max, Jn, .)
    cCG         = J(G*T_max, Jn, .)
    Wmat        = J(M, Jn, .)
    weights_all = J(max((T0-1,0)), Jn, .)
    target_w    = J(M, 1, .)

    //----------------------------------------------------------
    // 3. Main bootstrap loop
    //----------------------------------------------------------
    for (b_loop=1; b_loop<=boots; b_loop++) {

        // resample every cell and rebuild quantiles/CDFs
        for (t_loop=1; t_loop<=T_max; t_loop++) {
            bs = disco_boot_sorted(*(P[t_loop,1]))
            tQG[, t_loop] = disco_quantile_sorted(bs, pG)
            tCG[, t_loop] = disco_cdf_sorted(bs, gridG)
            if (t_loop<=T0-1) {
                target_w = (mixture==1 ? disco_cdf_sorted(bs, gridM) :
                                         disco_quantile_sorted(bs, pM))
            }

            for (j_loop=1; j_loop<=Jn; j_loop++) {
                bs = disco_boot_sorted(*(P[t_loop,1+j_loop]))
                cQG[|(t_loop-1)*G+1, j_loop \ t_loop*G, j_loop|] = disco_quantile_sorted(bs, pG)
                cCG[|(t_loop-1)*G+1, j_loop \ t_loop*G, j_loop|] = disco_cdf_sorted(bs, gridG)
                if (t_loop<=T0-1) {
                    Wmat[, j_loop] = (mixture==1 ? disco_cdf_sorted(bs, gridM) :
                                                   disco_quantile_sorted(bs, pM))
                }
            }

            // weights from the resampled pre-treatment data
            if (t_loop<=T0-1) {
                if (mixture==0) {
                    weights_all[t_loop, .] = disco_solve_weights(Wmat, target_w, simplex)
                }
                else {
                    weights_all[t_loop, .] = disco_mixture_weights(Wmat, target_w, simplex)'
                }
            }
        }

        // average pre-treatment weights
        if (T0-1>0) {
            W_avg = (colsum(weights_all)/(T0-1))'
        } else {
            // corner case T0=1 => fallback uniform
            W_avg = J(Jn,1,1/Jn)
        }

        // final synthetic distribution per period
        for (t_loop=1; t_loop<=T_max; t_loop++) {
            if (mixture==0) {
                Q_synth = cQG[|(t_loop-1)*G+1, 1 \ t_loop*G, Jn|] * W_avg
                C_synth = disco_cdf_sorted(sort(Q_synth,1), gridG)
            }
            else {
                C_synth = cCG[|(t_loop-1)*G+1, 1 \ t_loop*G, Jn|] * W_avg
                Q_synth = J(G,1,.)
                for (m8=1; m8<=G; m8++) {
                    p8  = (m8-1)/(G-1)
                    gg8 = 1
                    while (gg8<length(C_synth) & C_synth[gg8]<p8) gg8++
                    if (gg8>G) gg8 = G
                    Q_synth[m8] = gridG[gg8]
                }
            }

            quantile_diff_boot[|(t_loop-1)*G+1, b_loop \ t_loop*G, b_loop|] = tQG[, t_loop] :- Q_synth
            cdf_diff_boot[|(t_loop-1)*G+1, b_loop \ t_loop*G, b_loop|] = tCG[, t_loop] :- C_synth
        }
    }

    //----------------------------------------------------------
    // 4. Allocate final CI matrices
    //----------------------------------------------------------
    qdiff_lower = J(G,T_max,.)
    qdiff_upper = J(G,T_max,.)
    cdiff_lower = J(G,T_max,.)
    cdiff_upper = J(G,T_max,.)

    //----------------------------------------------------------
    // 5. Build intervals
    //----------------------------------------------------------
    if (uniform==0) {
        // Pointwise intervals
        for (t_loop2=1; t_loop2<=T_max; t_loop2++) {
            // quantile
            for (idx_loop3=1; idx_loop3<=G; idx_loop3++) {
                vals       = quantile_diff_boot[(t_loop2-1)*G + idx_loop3,.]'
                tmp        = sort(vals, 1)

                alpha      = (1-cl)/2
                lower_idx  = max((ceil(alpha*boots), 1))
                upper_idx  = min((ceil((1-alpha)*boots), boots))

                qdiff_lower[idx_loop3,t_loop2] = tmp[lower_idx]
                qdiff_upper[idx_loop3,t_loop2] = tmp[upper_idx]
            }

            // cdf
            for (idx_loop4=1; idx_loop4<=G; idx_loop4++) {
                vals       = cdf_diff_boot[(t_loop2-1)*G + idx_loop4,.]'
                tmp        = sort(vals, 1)

                alpha      = (1-cl)/2
                lower_idx  = max((ceil(alpha*boots), 1))
                upper_idx  = min((ceil((1-alpha)*boots), boots))

                cdiff_lower[idx_loop4,t_loop2] = tmp[lower_idx]
                cdiff_upper[idx_loop4,t_loop2] = tmp[upper_idx]
            }
        }
    }
    else {
        // Uniform intervals
        for (t_loop3=1; t_loop3<=T_max; t_loop3++) {
            // quantile
            qdiff_mat = rowshape(
                quantile_diff_boot[(t_loop3-1)*G+1 .. t_loop3*G, .], G
            )
            qdiff_err = abs(qdiff_mat :- quantile_diff[., t_loop3])
            qmax_abs  = colmax(qdiff_err)
            sorted_diffs_q = sort(qmax_abs', 1)

            alpha_q = cl
            index_q = floor(alpha_q * boots)
            if (index_q<1)   index_q=1
            if (index_q>boots) index_q=boots

            cval_q = sorted_diffs_q[index_q]
            for (m_i2=1; m_i2<=G; m_i2++) {
                base_val2 = quantile_diff[m_i2,t_loop3]
                qdiff_lower[m_i2,t_loop3] = base_val2 - cval_q
                qdiff_upper[m_i2,t_loop3] = base_val2 + cval_q
            }

            // cdf
            cdf_mat = rowshape(
                cdf_diff_boot[(t_loop3-1)*G+1 .. t_loop3*G, .], G
            )
            cdf_err = abs(cdf_mat :- cdf_diff[., t_loop3])
            cmax_abs = colmax(cdf_err)
            sorted_diffs_c = sort(cmax_abs', 1)

            alpha_q = cl
            index_q = floor(alpha_q * boots)
            if (index_q<1)   index_q=1
            if (index_q>boots) index_q=boots

            cval_c = sorted_diffs_c[index_q]
            for (c_i2=1; c_i2<=G; c_i2++) {
                base_val_c = cdf_diff[c_i2,t_loop3]
                cdiff_lower[c_i2,t_loop3] = base_val_c - cval_c
                cdiff_upper[c_i2,t_loop3] = base_val_c + cval_c
            }
        }
    }

    //----------------------------------------------------------
    // 6. Return final CI matrices
    //----------------------------------------------------------
    co.qdiff_lower = qdiff_lower
    co.qdiff_upper = qdiff_upper
    co.cdiff_lower = cdiff_lower
    co.cdiff_upper = cdiff_upper
    return(co)
}



// -----------------------------------------------------------------------------
// Wrapper for main DISCO (kept same but with M & G as arguments).
// -----------------------------------------------------------------------------
real scalar disco_wrapper(real vector y, id, tt,
                         real scalar target_id,
                         real scalar T0, real scalar T_max,
                         real scalar M, real scalar G,
                         real scalar q_min, real scalar q_max,
                         real scalar simplex, real scalar mixture)
{
    struct disco_out scalar results

    results = disco_full_run(y, id, tt, target_id, T0, T_max, M, G,
                             q_min, q_max, simplex, mixture)

    st_matrix("weights", results.weights')
    st_matrix("quantile_diff", results.quantile_diff)
    st_matrix("cdf_diff", results.cdf_diff)
    st_matrix("quantile_synth", results.quantile_synth)
    st_matrix("quantile_t", results.quantile_t)
    st_matrix("cdf_synth", results.cdf_synth)
    st_matrix("cdf_t", results.cdf_t)
    st_matrix("cids", results.cids')

    return(0)
}

// -----------------------------------------------------------------------------
// Wrapper for CI now passing M & G to disco_bootstrap_CI
// -----------------------------------------------------------------------------
real scalar disco_ci_wrapper(real vector y, id, tt,
                           real scalar target_id,
                           real scalar T0, real scalar T_max,
                           real scalar M, real scalar G,
                           real scalar q_min, real scalar q_max,
                           real scalar simplex, real scalar mixture,
                           real scalar boots, real scalar cl, real scalar uniform,
                           real matrix quantile_diff, real matrix cdf_diff)
{
    struct CI_out scalar results

    results = disco_bootstrap_CI(y, id, tt, target_id, T0, T_max, M, G,
                                 q_min, q_max, simplex, mixture,
                                 boots, cl, uniform,
                                 quantile_diff, cdf_diff)

    st_matrix("qdiff_lower", results.qdiff_lower)
    st_matrix("qdiff_upper", results.qdiff_upper)
    st_matrix("cdiff_lower", results.cdiff_lower)
    st_matrix("cdiff_upper", results.cdiff_upper)

    return(0)
}
// Compute summary stats
real scalar compute_summary_stats(string scalar agg, real vector sample_points,
                                real scalar T0, real scalar T_max, real matrix quantile_diff,
                                real matrix cdf_diff, real scalar CI, real scalar cl)
{
    if (!anyof(("quantile", "cdf", "quantileDiff", "cdfDiff"), agg)) {
        errprintf("Invalid aggregation type\n")
        return(1)
    }

    real scalar is_cdf
    is_cdf = (agg == "cdf" | agg == "cdfDiff")

    real vector grid_points
    if (is_cdf) {
        real scalar amin, amax
        amin = st_numscalar("amin")
        amax = st_numscalar("amax")

        real vector grid_points_temp
        if (max(sample_points) <= 1 & min(sample_points) >= 0) {
            grid_points_temp = amin :+ sample_points :* (amax - amin)
            if (min(grid_points_temp) > amin) grid_points = amin \ grid_points_temp'
            else grid_points = grid_points_temp'
            if (max(grid_points_temp) < amax) grid_points = grid_points \ amax
        } else {
            grid_points = sort(sample_points', 1)
        }
    } else {
        grid_points = sort(sample_points', 1)
        if (min(grid_points) > 0) grid_points = 0 \ grid_points
        if (max(grid_points) < 1) grid_points = grid_points \ 1
    }

    real scalar n_intervals
    n_intervals = length(grid_points) - 1

    real matrix summary_stats
    summary_stats = J(n_intervals * (T_max - T0 + 1), 7, .)

    real scalar row, G, M, t, i
    row = 1

    real vector grid, idx, prob_grid

    for(t = T0; t <= T_max; t++) {
        for(i = 1; i <= n_intervals; i++) {
            summary_stats[row, 1] = t
            summary_stats[row, 2] = grid_points[i]
            summary_stats[row, 3] = grid_points[i + 1]

            if (is_cdf) {
                amin = st_numscalar("amin")
                amax = st_numscalar("amax")
                G = rows(cdf_diff)
                grid = range(amin, amax, (amax - amin)/(G-1))'
                idx = selectindex(grid :>= grid_points[i] :& grid :<= grid_points[i + 1])
            } else {
                M = rows(quantile_diff)
                prob_grid = range(0, 1, 1/(M-1))'
                idx = selectindex(prob_grid :>= grid_points[i] :& prob_grid :<= grid_points[i + 1])
            }

            if (length(idx) > 0) {
                if (is_cdf) {
                    summary_stats[row, 4] = mean(cdf_diff[idx, t])
                } else {
                    summary_stats[row, 4] = mean(quantile_diff[idx, t])
                }
            }

            if (CI) {
                real matrix diff_lower, diff_upper
                if (is_cdf) {
                    diff_lower = st_matrix("cdiff_lower")
                    diff_upper = st_matrix("cdiff_upper")
                } else {
                    diff_lower = st_matrix("qdiff_lower")
                    diff_upper = st_matrix("qdiff_upper")
                }

                if (length(idx) > 0) {
                    if (is_cdf) {
                        summary_stats[row, 4] = mean(cdf_diff[idx, t])
                        summary_stats[row, 6] = mean(diff_lower[idx, t])
                        summary_stats[row, 7] = mean(diff_upper[idx, t])
                        summary_stats[row, 5] = (summary_stats[row, 7] - summary_stats[row, 6])/(2*1.96)
                    } else {
                        summary_stats[row, 4] = mean(quantile_diff[idx, t])
                        summary_stats[row, 6] = mean(diff_lower[idx, t])
                        summary_stats[row, 7] = mean(diff_upper[idx, t])
                        summary_stats[row, 5] = (summary_stats[row, 7] - summary_stats[row, 6])/(2*1.96)
                    }
                }
            }

            row++
        }
    }

    st_matrix("summary_stats", summary_stats)

    return(0)
}

end
********************************************************************************
*** END MATA CODE
********************************************************************************
