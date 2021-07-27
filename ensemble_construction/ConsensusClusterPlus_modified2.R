################################################################################
# This script is from the ‘ConsensusClusterPlus’ bioconductor package.
# maxK parameter is modified so that the function only takes one k as input.

ConsensusClusterPlus.M <- function( d=NULL,
                                    oneK = (2:3),
                                    reps=10,
                                    pItem=0.8,
                                    pFeature=1,
                                    clusterAlg="hc",
                                    title="untitled_consensus_cluster",
                                    innerLinkage="average",
                                    finalLinkage="average",
                                    distance="pearson",
                                    ml=NULL,
                                    seed=NULL,
                                    weightsItem=NULL,
                                    weightsFeature=NULL,
                                    verbose=F,
                                    corUse="everything" ) {
  ##description: runs consensus subsamples
  if(is.null(seed)==TRUE){
    seed=timeSeed = as.numeric(Sys.time())
  }
  set.seed(seed)
  message(paste("cluster alg: ",clusterAlg))
  #distance=ifelse( inherits(d,"dist"), attr( d, "method" ), "pearson" )


  if(is.null(ml)==TRUE){

    if ( ! class( d ) %in% c( "dist", "matrix", "ExpressionSet" ) ) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }

    if ( inherits( d, "dist" ) ) {
      ## if d is a distance matrix, fix a few things so that they don't cause
      ## problems with the analysis
      ## Note, assumption is that if d is a distance matrix, the user doesn't
      ## want to sample over the row features
      if ( is.null( attr( d, "method" ) ) ) {
        attr( d, "method" ) <- distance <- "unknown - user-specified"
      }
      if ( is.null( distance ) || ( distance != attr( d, "method" ) ) ) {
        distance <- attr( d, "method" )
      }

      if ( ( ! is.null( pFeature ) ) && ( pFeature < 1 ) ) {
        message("Cannot use the pFeatures parameter when specifying a",
                "distance matrix as the data object\n")
        pFeature <- 1
      }
      if ( ! is.null( weightsFeature ) ) {
        message("Cannot use the weightsFeature parameter when specifying",
                "a distance matrix as the data object\n")
        weightsFeature <- NULL
      }
      if ( clusterAlg == "km" ) {
        message("Note: k-means will cluster the distance matrix you",
                "provided.  This is similar to kmdist option when",
                "suppling a data matrix")
        ##d <- as.matrix( d )  #this is now done w/in ccRun
      }
    } else {
      if ( is.null( distance ) ) {
        ## we should never get here, but just in case
        distance <- "pearson"
      }
    }

    if ( ( clusterAlg == "km" ) && inherits( distance, "character" ) &&
         ( distance != "euclidean" ) ) {
      message("Note: The km (kmeans) option only supports a euclidean distance",
              "metric when supplying a data matrix.  If you want to cluster a",
              "distance matrix using k-means use the 'kmdist' option, or use a",
              "different algorithm such as 'hc' or 'pam'.  Changing distance",
              "to euclidean")
      distance <- 'euclidean'
    }


    if ( inherits( d,"ExpressionSet" ) ) {
      d <- exprs(d)
    }

    ml <- ccRun( d=d,
                 oneK=oneK,
                 repCount=reps,
                 diss=inherits(d,"dist"),
                 pItem=pItem,
                 pFeature=pFeature,
                 innerLinkage=innerLinkage,
                 clusterAlg=clusterAlg,
                 weightsFeature=weightsFeature,
                 weightsItem=weightsItem,
                 distance=distance,
                 verbose=verbose,
                 corUse=corUse)
  }
  res=list();

  for (tk in 2:(length(oneK)+1)){
    if(verbose){
      message(paste("consensus ",oneK[tk-1]))
    }
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=finalLinkage);
    message("clustered")
    ct = cutree(hc,oneK[tk-1])
    names(ct) = colnames(d)
    if(class(d)=="dist"){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm

    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,
                     ml=ml[[tk]])
  }

  return(res)
}


ccRun <- function( d=d,
                   oneK=NULL,
                   repCount=NULL,
                   diss=inherits( d, "dist" ),
                   pItem=NULL,
                   pFeature=NULL,
                   innerLinkage=NULL,
                   distance=NULL,
                   clusterAlg=NULL,
                   weightsItem=NULL,
                   weightsFeature=NULL,
                   verbose=NULL,
                   corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",(length(oneK)+1))
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);

  ## necessary if d is a dist object and attr( d, "method" ) == NULL
  if (is.null( distance ) ) distance <- 'euclidean'
  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra",
                            "binary", "minkowski", "pearson", "spearman" )

  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set
    ## (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
           ( pFeature < 1 ) ) {
      message("user-supplied data is a distance matrix; ignoring",
              "user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring",
               "user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else if(FALSE){ ## d is a data matrix
    #message("d is a data matrix")
    #exit()
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
           ( is.null( pFeature ) ||
               ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not
      #  sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &
             ( class(try(get(distance),silent=T))!="function") )
          stop("unsupported distance.")

        if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
        }else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene
      ## features, so main.dist.obj is left as NULL
    }
  }


  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )
    message(paste("length subcols: ", length(sample_x$subcols)))
    this_dist = NA
    message(!is.null( main.dist.obj ))
    if (  is.null( main.dist.obj ) ) {
      message("boot cols: ")
      boot.cols <- sample_x$subcols
      message(length(boot.cols))
      this_dist <- as.matrix(d)[ boot.cols, boot.cols ]
      message(dim(this_dist))
      this_dist <- this_dist * (diag(rep(-1,length(boot.cols)))+1)
      #this_dist <- apply(this_dist, c(1,2), FUN = function(x) round(x, digits=10))
      message(paste("this distt: ", this_dist[1,1]))
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        #this_dist <- as.dist( this_dist )
        this_dist <- this_dist
        #message(paste("this distt: ", this_dist[1,1]))
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {
      ## if main.dist.obj is NULL, then d is a data matrix, and either:
      ##   1) clusterAlg is 'km'
      ##   2) pFeatures < 1 or weightsFeatures have been specified, or
      ##   3) both
      ## so we can't use a main distance object and for every iteration, we will
      ## have to re-calculate either
      ##   1) the distance matrix (because we're also sampling the features as
      ##      well), or
      ##   2) the submat (if using km)

      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &
             ( class(try(get(distance),silent=T))!="function")  )
          stop("unsupported distance.")
        if( ( class(try(get(distance),silent=T))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
        }else{
          if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,
                                        method=distance) )
          }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
          }
        }
        attr( this_dist, "method" ) <- distance
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
               ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }

          this_dist <- sample_x$submat
        }
      }
    }

    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ## mCount is possible number of times that two sample occur in same random
    ## sample, independent of k
    ## mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] )

    ##use samples for each k
    for (k in 2:(length(oneK)+1)){
      if(verbose){
        message(paste("  k =",k))
        #message(paste("  oneK =",length(oneK)))
        #message(paste("  K =",oneK[k-1]))
      }
      if (i==1){
        #message("pre init")
        ml[[k]] = mConsist #initialize
        #message("post init")
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        message("hc")
        ##prune to k for hc
        this_assignment = cutree(this_cluster,oneK[k-1])

      }else if(clusterAlg=="kmdist"){
      message("kmdist")
        this_assignment = kmeans(this_dist, oneK[k-1], iter.max = 10,
                                 nstart = 1,
                                 algorithm = c("Hartigan-Wong") )$cluster

      }else if(clusterAlg=="km"){
        message("km")
        ## this_dist should now be a matrix corresponding to the result from
        ## sampleCols
        this_assignment <- kmeans( t( this_dist ),
                                   oneK[k-1],
                                   iter.max = 10,
                                   nstart = 1,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if (clusterAlg=="pam"){
        message("pam")
        this_assignment <- pam( x=this_dist,
                                oneK[k-1],
                                diss=TRUE,
                                metric=distance,
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        message("running hook")
        #message(dim(this_dist))
        this_assignment <- get(clusterAlg)(this_dist, oneK[k-1])
        message("ran hook")
      }
      ##add to tally
      message("begin connect assign")
      ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     sample_x[[3]] )
      message("end connect assign")
    }
  }


  ##consensus fraction
  res = vector(mode="list",(length(oneK)+1))
  for (k in 2:(length(oneK)+1)){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  message("end fraction")
  return(res)
}


connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  ## input: named vector of cluster assignments, matrix to add connectivities
  ## output: connectivity matrix
  names( clusterAssignments ) <- sampleKey
  # list samples by clusterId
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric(
    names( clusterAssignments[ clusterAssignments %in% i ] ) ) )

  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    # product of arrays with * function; with above indicator (1/0) statement
    # updates all cells to indicate the sample pair was observed int the same
    # cluster;
    updt <- outer( cl, cl )
    m <- m + updt
  }
  return(m)
}



sampleCols <- function( d,
                        pSamp=NULL,
                        pRow=NULL,
                        weightsItem=NULL,
                        weightsFeature=NULL ){
  ## returns a list with the sample columns, as well as the sub-matrix & sample
  ## features (if necessary)
  ## if no sampling over the features is performed, the submatrix & sample
  ## features are returned as NAs
  ## to reduce memory overhead

  print(paste("d dim: ",dim(d)))
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  message(paste("space: ",space))
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem))

  this_sample <- sampRows <- NA


  if ( inherits( d, "matrix" ) ) {
    message("matrixxxxxxxxx")
    weightsFeature=NULL
    message(paste("if weights: ",(! is.null( weightsFeature ) )))
    if ( (! is.null( pRow ) ) &&
           ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over
      ## the row/gene/features
      space = nrow(d)
      message(paste("space2: ",space))
      sampleN = floor(space*pRow)
      print(sampleN)
      sampRows = sort( sample(space, sampleN, replace = FALSE,
                              prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}



triangle <- function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m


  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half

  fm = t(nm)+nm
  diag(fm) = diag(m)

  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]

  if(mode==1){
    return(vm) #vector
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }

}
