fast_confSAM_1 = function(p, PM, includes.id = TRUE,
                           cutoff = 0.01, reject = "small", alpha = 0.05,
                           method="simple",  ncombs=1000, verbose_from = 1) {

  if (ncol(PM) != length(p)){
    stop("invalid permutation matrix")
  }

  if (includes.id & !(min(PM[1,] == p))) {
    stop("first row of matrix provided does not equal vector p provided.")
  }

  if (length(cutoff) != 1 & length(cutoff) != length(p)) {
    stop("length of cutoff should be 1 or length(p)")
  }

  # if identity permutation (p) not included in PM: add it to PM
  if (!includes.id) {
    PM = cbind(p, PM)
  }

  # Dimensions of PM
  w = nrow(PM)    # each row corresponds to a permutation
  m = ncol(PM)    # each column corresponds to a test

  # rejection function
  if (reject== "small") {
    reject_func = function(x) x < cutoff
  }

  if (reject== "large"){
    reject_func = function(x) x > cutoff
  }

  if (reject== "absolute"){
    reject_func = function(x) abs(x) > cutoff
  }

  reject_num  = function(x) sum(reject_func(x))

  # Which of the identity permutations exceed the cutoff?
  Rset = reject_func(p)

  # Partition tests based on whether identity permutation rejects:
  # number of rejections for each permutation in both partitions
  nrej_R = apply(PM[, Rset], 1, reject_num)
  nrej_Rc  = apply(PM[, !Rset],  1, reject_num)

  # total number of rejections for each permutation
  nrej = nrej_R + nrej_Rc

  if (nrej[1] == 0) {
    stop("Error: No rejections.")
  }

  # Critical value
  k = ceiling((1 - alpha) * w)


  ## Simple
  simple = min(sort(nrej, partial = k)[k], nrej[1])
  est = min(sort(nrej, partial = floor(w / 2))[floor(w / 2)], nrej[1])

  if (method == "simple") {
    out = c(nrej[1], est, simple)
    names(out) = c("#rejections:", "Simple estimate of #fp:",
                    "Simple conf. bound for #fp:")
    return(out)
  }


  ## Approx

  if (method == "approx") {

    # Initializing upper bound
    bound = nrej[1]
    last_i = 1

    # initialize l as follows
      # explanation: sum(l > nrejs) <= sum(l > nrej_Rc), for all l, as
      # nrejs = nrej_Rc + nrejs_rcombs >= nrej_Rc (element-wise)
      # so, sum(l > nrej_Rc) < k implies sum(l > nrejs) < k
      # hence, we only need to check l >= start_l
    start_l = sort(nrej_Rc, partial = k)[k]

    if (start_l >= nrej[1]) {
      out <- c(nrej[1], est, min(bound, simple), ncombs)
      names(out) <- c("#rejections:", "Simple estimate of #fp:",
                      "Appr. cl.testing-based bound for #fp:", "ncombs")
      return(out)
    }

    # Generate the random combs from rejection set
    indR = which(Rset)
    rcombs = replicate(ncombs, sample(indR, size = nrej[1], replace = FALSE))
    if (is.null(dim(rcombs))) {
      rcombs = matrix(rcombs, nrow = 1)
    }

    # Find the smallest value of l so that for all larger values of l,
    # no comb leads to a rejection
    for (l in start_l:nrej[1]) {
      any_reject = FALSE

      # Is any of the rcombs rejected?
      for (i in last_i:ncombs) { # can start at last_i due to explanation below
        elements = rcombs[1:l, i]
        PM_temp = PM[, elements, drop = F]

        nrejs_rcombs = apply(PM_temp, 1, reject_num)

        # use the precomputation of nrej_Rc
        nrejs = nrej_Rc + nrejs_rcombs

        # Reject this rcomb?
        if (sum(l > nrejs) < k) {
          any_reject = TRUE
          break
        }
      }

      # NOTE: for fixed i, sum(l > nrejs) weakly increases in l
      # explanation: l increases by 1 and each element of nrejs by at most 1
      # if they all increase by 1, sum(l > nrejs) remains unchanged
      # if some do not increase by 1, sum(l > nrejs) may increase or not change
      # so, if l increases, we don't have to re-check i's we already checked
      last_i = i

      # Report progress
      if (i > verbose_from * ncombs) {
        progress = paste0("l: ", l, ". i: ", i)
        print(progress)
      }

      if (!any_reject) {
        bound = l - 1
        break
      }
    }

    out <- c(nrej[1], est, min(bound, simple), ncombs)
    names(out) <- c("#rejections:", "Simple estimate of #fp:",
                  "Appr. cl.testing-based bound for #fp:", "ncombs")
    return(out)
  }
}
