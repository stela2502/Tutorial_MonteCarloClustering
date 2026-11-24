
MonteCarloData <- R6::R6Class("MonteCarloData",
  public = list(
    data = NULL,
    vars_old = NULL,
    vars_new= NULL,
    changed_clusters=NULL,
    drc=NULL,
    Mgroups = NULL,
    T=NULL,
    c=NULL,
    Eold = NULL,
    failed=0,
    good=0,
    initialize = function(T=1, c= .999, drc=NULL) {
      self$data = Matrix::t(cellexalvrR::cellexalObj@data)
      self$Mgroups = cbind( colnames(cellexalvrR::cellexalObj@data), sample(as.character(1:10), ncol(cellexalvrR::cellexalObj@data), replace=TRUE) )
      self$vars_old = lapply( split(self$Mgroups[,1], self$Mgroups[,2]), function(group){ 
        var( matrix(self$data[group,]))
        } )
      self$changed_clusters = c()
      self$c=c
      self$T = T
      self$Eold = +Inf
      self$failed=0
      self$drc = drc
    },
    print = function(){
      cat("This is an example class R6::MonteCarloData\n")
      cat(paste("the data table has", ncol(self$data),"columns and", nrow(self$data), "rows\n" ) )
      cat( paste( "and the Mgrouping has", length( unique(self$Mgroups[,2])), "groups\n"  ) )
    },
    cluster =function( startT=NULL,  ){

      i = 0;
      pb <- progress::progress_bar$new( total=steps )
      for ( i in 1:steps ){
        pb$tick()
        self$r_shift()
      }
      f()
      print ("Finished")
    },
    r_shift= function( ){
      #set.seed( as.numeric(paste( sep="","1", stringr::str_sub(as.character(round(as.numeric(Sys.time())*10000)), -5,-1)) ))
      old = self$Mgroups
      A = sample( rownames(self$data), 1 )
      B = sample( rownames(self$data), 1 )
      
      if ( self$Mgroups[match(A,self$Mgroups[,1]),2] != self$Mgroups[match(B,self$Mgroups[,1]),2] ){
        ## do the shift
        #print( paste("change",A,"to", B, "from", self$Mgroups[match(A,self$Mgroups[,1]),2], "to",  self$Mgroups[match(B,self$Mgroups[,1]),2]) )
        tmp = self$Mgroups[match(A,self$Mgroups[,1]),2]
        self$Mgroups[match(A,self$Mgroups[,1]),2] = self$Mgroups[match(B,self$Mgroups[,1]),2]
        self$Mgroups[match(B,self$Mgroups[,1]),2] = tmp

        self$changed_clusters = c( self$Mgroups[match(A,self$Mgroups[,1]),2], 
          self$Mgroups[match(B,self$Mgroups[,1]),2] ) 

        Enew = self$getEnergy( )
        ## OK this is far too simple here. The temperature cooling algorithm needs to be added.
        if ( Enew < self$Eold ) {
          self$Eold = Enew
          self$good = self$good + 1
        }else if (  exp( - ( Enew- data$Eold) / self$T) > runif(1) ) {
          self$Eold = Enew
          self$good = self$good + 1
        }else{
          self$failed =  self$failed +1 
        }
      }
      self$T = self$T * self$c
    },
    getEnergy  = function(  ){
      self$vars_new = self$vars_old

      if ( length( self$changed_clusters) != 2) {
       stop("I need exactly two changed clusters here!")
      }
      self$vars_new[[ self$changed_clusters[1] ]] = var( matrix(self$data[ match( self$changed_clusters[1], self$Mgroups[,2]), ]))
      self$vars_new[[ self$changed_clusters[2] ]] = var( matrix(self$data[ match( self$changed_clusters[2], self$Mgroups[,2]), ]))
      return( mean(unlist(self$vars_new)) )

    },
    plot = function(main =""){
        plot(drc, col=rainbow(length(unique(self$Mgroups[,2])))[as.numeric(self$Mgroups[order(self$Mgroups[,1]),2])], main=main )
    }
  )
)

