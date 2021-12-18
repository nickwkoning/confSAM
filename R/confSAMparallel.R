confSAMparallel <- function(p, PM, includes.id=TRUE,
                            cutoff=0.01, reject="small", alpha=0.05,
                            method="simple",  ncombs=1000, ncores = 1) {

  if (ncol(PM)!=length(p) & nrow(PM)!=length(p)){
    stop("invalid permutation matrix")
  }

  if (ncol(PM)!=length(p) & nrow(PM)==length(p)){
    PM<-t(PM)
  }

  w <- nrow(PM)    #each row corresponds to a perm
  m <- ncol(PM)

  if(includes.id & !(min(PM[1,]==p))) {
    stop("first row/column of matrix provided does not equal vector p provided.")
  }

  if( length(cutoff)!=1 & length(cutoff)!=length(p) ) {
    stop("length of cutoff should be 1 or length(p)")
  }

  if(!includes.id){
    PMid <- matrix(nrow=w+1,ncol=m)
    PMid[2:(w+1),] <- PM
    PMid[1,] <- p
    PM <- PMid
    w<-nrow(PM)  # i.e. w <- w+1
  }




  k <- ceiling((1-alpha)*w)

  if(reject== "small"){
    nrej <- apply( PM, 1, function(x) {sum(x<cutoff)} )
  }
  if(reject== "large"){
    nrej <- apply( PM, 1, function(x) {sum(x>cutoff)} )
  }
  if(reject== "absolute"){
    nrej <- apply( PM, 1, function(x) {sum(x>cutoff)+sum(-x>cutoff) } )
  }




  simple <- min( sort(nrej, partial = k)[k] , nrej[1] )
  est <- min( sort(nrej, partial = floor(0.5*w))[floor(0.5*w)] , nrej[1] )



  if(method=="simple"){
    out <- c(nrej[1], est, simple)
    names(out) <- c("#rejections:", "Simple estimate of #fp:",
                    "Simple conf. bound for #fp:")
    return(out)
  }


  if(method=="full" | method=="approx"| method=="csc"){
    #make vector indR with indices corresponding to rejected set:
    if(reject== "small"){
      Rset <- (p < cutoff)
    }
    if(reject== "large"){
      Rset <- (p > cutoff)
    }
    if(reject== "absolute"){
      Rset <- (abs(p) > cutoff)
    }


    indR = which(Rset)
    indRc = which(!Rset)
  }


  if(reject== "small"){
    reject_fun = function(x) sum(x<cutoff)
  }
  if(reject== "large"){
    reject_fun = function(x) sum(x>cutoff)
  }
  if(reject== "absolute"){
    reject_fun = function(x) sum(abs(x)>cutoff)
  }


  if(method=="approx"){
    if(ncombs > 1e6){ message("The procedure might be computationally
    infeasible since ncombs is very large.")}

    if (file.exists("upper_bounds.txt")) {
      file.remove("upper_bounds.txt")
    }

    if (file.exists("stop.txt")) {
      file.remove("stop.txt")
    }

    write(0, "upper_bounds.txt", append = TRUE)

    zap <- foreach(
      j = 1:ncores,
      .combine= 'c'
    ) %dopar% {

      # spread values of l evenly over the cores
      for (l in which(1:nrej[1] %% ncores == j)) {

        # Check if lowest upper bound has been exceeded
        if (file.exists("stop.txt")) {
          stopping_l = unlist(read.delim("stop.txt", header = FALSE))

          # If l is smaller than the lowest upper bound found so far,
          # then the current l may still yield a lower upper bound
          if (l >= stopping_l) {
            return(NULL)
          }
        }

        # Generate the random combs from rejection set
        rcombs = replicate(ncombs, sample(indR, size=l, replace=FALSE))
        if (is.null(dim(rcombs))) {
          rcombs = matrix(rcombs, nrow = 1)
        }

        for (i in 1:ncombs) {
          nrejs = apply( PM[, c(indRc, rcombs[,i])], 1, reject_fun)

          if (sum(l > nrejs) < k) {
            content = unlist(read.delim("upper_bounds.txt", header = FALSE))

            if (l >= max(content)) {
              write(l, "upper_bounds.txt")
              break
            }
          }
        }

        # lowest upper bound has been exceeded
        if (i == ncombs) {
            write(l, "stop.txt")
        }
      }
      return(l)
    } #loop l

    appctbound = unlist(read.delim("upper_bounds.txt", header = FALSE))
    file.remove("upper_bounds.txt")

    if (file.exists("stop.txt")) {
      file.remove("stop.txt")
    }

    out <- c(nrej[1],est,min(appctbound, simple))
    names(out) <- c("#rejections:", "Simple estimate of #fp:",
                    "Appr. cl.testing-based bound for #fp:")
    return(out)

  }

  if(method=="csc"){  #conservative shortcut

    if(reject != "small"){ stop("The conservative shortcut is only useful if the test
    statistics are p-values and the smallest p-values
    are rejected") }

    ord <- order(nrej)

    S <- (PM[ord,]<cutoff)%*% Rset
    U <- numeric(nrej[1])

    Vsc <- 0
    M <- simple    #start checking for M=simple and then lower M
    found <- FALSE

    while(found==FALSE){
      .s <- 0
      if(sort(p)[1]<cutoff){
        if(sum(Rset)>1){
          .s <- apply( (PM[,Rset==TRUE] <cutoff) ,2,sum)
        }
      }

      SIGMA <- 0
      if(nrej[1]-simple>0){
        SIGMA  <- sum(  (sort(.s))[1:(nrej[1]-M)]  )
      }

      #calculation of maxA:

      maxAfound <- FALSE
      s <- 0
      nrej.sorted <- sort(nrej)

      while(maxAfound==FALSE){
        N_s <- sum(nrej.sorted < nrej.sorted[k]-s)
        M_s <- k-1-N_s
        Ks <- pmax(0, (S- nrej.sorted + nrej.sorted[k]-s)[(N_s+1):w] )
        Ks.sorted <- sort(Ks, decreasing=TRUE)
        sum1<- sum(S[1:N_s])
        sum2 <- sum( pmin(S[(N_s+1):w], nrej.sorted[(N_s+1):w]-nrej.sorted[k]+s)    )
        sum3 <- sum(Ks.sorted[1:M_s])
        Maxeraf <- sum1+sum2+sum3
        if(SIGMA <= Maxeraf){
          maxAfound <- TRUE; maxA <- s-1
        }
        s <- s+1
      }

      U[M] <- min( nrej[1], nrej.sorted[k]-1-maxA )

      if(M<=U[M]){
        found <- TRUE
        Vsc <- min(M,simple,nrej[1])
      }
      M <- M-1
    }

    out <- c(nrej[1],est,Vsc)
    names(out) <- c("#rejections:", "Simple estimate of #fp:",
                    "Bound #fp based on shortcut:")
    return(out)

  }

} # end function confSAM

