  # Dimensions of PM
  w = nrow(PM)    # each row corresponds to a permutation
  m = ncol(PM)    # each column corresponds to a test

  # rejection function
  reject_func = function(x) abs(x) > cutoff
  reject_num  = function(x) sum(reject_func(x))

  # Which of the identity permutation exceed the cutoff?
  Rset = reject_func(p)

  # Partition tests based on whether identity permutation rejects
  indR = which(Rset)
  indRc = which(!Rset)

  # number of rejections for each permutation in both partitions
  nrej_R = apply(PM[, Rset], 1, reject_num)
  nrej_Rc  = apply(PM[, !Rset],  1, reject_num)

  # total number of rejections for each permutation
  nrej = nrej_R + nrej_Rc


  # Critical value
  k = ceiling((1 - alpha) * w)


  ## Simple
  simple <- min( sort(nrej, partial = k)[k] , nrej[1] )
  est <- min( sort(nrej, partial = floor(0.5*w))[floor(0.5*w)] , nrej[1] )

  tic()
  ## Approx

  # Initializing upper bound
  upper_bound = 0
  last_i = 1

  # Generate the random combs from rejection set
  rcombs = gen_rcombs(ncombs, indR, nrej[1])


  # Find the smallest value of l so that for all larger values of l,
  # no comb leads to a rejection
  for (l in 1:nrej[1]) {
    any_reject = FALSE
    for (i in last_i:ncombs) { # can start at last_i due to explanation below
      elements = rcombs[1:l, i]
      PM_temp = PM[, elements, drop = F]

      nrejs_rcombs = apply(PM_temp, 1, reject_num)
      nrejs = nrej_Rc + nrejs_rcombs

      # Reject?
        # NOTE: for fixed i, the below sum weakly increases in l
        # explanation: l increases by 1 and each element of nrejs by at most 1
        # if they all increase by 1, the sum remains unchanged
        # if some do not increase by 1, the sum may increase or remain unchanged
      if (sum(l > nrejs) < k) {
        any_reject = TRUE
        break
      }
    }

    last_i = i

    if (i > 10) {
      print(c(l, i, sum(l > nrej_Rc)))
    }

    if (!any_reject) {
      upper_bound = l - 1
      break
    }
  }

  appctbound = upper_bound

  out <- c(nrej[1],est,min(appctbound, simple), cutoff, ncombs)
  names(out) <- c("#rejections:", "Simple estimate of #fp:",
                "Appr. cl.testing-based bound for #fp:", "cutoff", "ncombs")
  out
