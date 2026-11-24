# Tutorial_MonteCarloClustering

## How to get there

Create a new git repo.
Cretae a new R package using the ``` Rcpp.package.skeleton('YourPackageName') ```.

Commit this to github.


### Get yourself a test dataset

Small enough to not be of great (computational) problem, but complicated enough to be a little bit chalenging.

```
library(Matrix)
set.seed(2785695)
len = 5000

dims = round(runif(100) *10000)
dims = dims[ which(dims > 500 )]

dat =  unlist( lapply( 1:5, function(i) { sample(x = 1:dims[i], size = len/5, replace = TRUE) } ))

#dat [sort(sample(x = 1:len, size = round( len * .7) ))] = 0
dat = matrix(dat, nrow=100) ## we will cluster these 100 rows
dim(dat)

data = Matrix( dat, sparse=TRUE)

dim(data)
length(data@x)


set.seed(5)
len = 5000

dims = round(runif(100) *100000)
dims = dims[ which(dims > 5000 )]

dat =  unlist( lapply( 1:5, function(i) { sample(x = 1:dims[i], size = len/5, replace = TRUE) } ))

#dat [sort(sample(x = 1:len, size = round( len * .7) ))] = 0
dat = matrix(dat, nrow=100) ## we will cluster these 100 rows
dim(dat)

data2 = Matrix( dat, sparse=TRUE)

dim(data2)
length(data2@x)

data = rbind(data, data2)

```

Now check what you got - it is not impressive:

```
pcr = prcomp( t(data) )
plot(pcr$x)
```

Possibly 5 groups? Lets check if we can use a MonteCarloc clustering on that data.

## Implement monte Carlo

We need a method to get the amount of order in the data.
Easiest might actually be a mean of col vars.

```
mean(apply( data[group,], 1, var))
```

And this summed up for all groups:

```
order = sample( 1:nrow(data), nrow(data) )
groups = split( order, rep( 1:2,  nrow(data)/2 ))

sum(unlist (lapply( groups, function(group){ mean(apply( data[group,], 1, var)) } )))
```

This all seams quite straight forward. Now check how we could implement the Monte Carlo part.
I am a 100% newbee on that. So pain is to be expected.

### The system should implement a simulated annealing process

Annealing is dependant on the temperature. So the higher the temperature the easier it is for the variables to change cluster.
So we need to store the "amount of order" we calulate each time in a histogram and then get from this histogram whether a move has been beneficial or bad.
This value better be between 0 and 1 so that we can have a simple overall temperature to check if it was a beneficial move.

I first need to get an impression of how this data might look:


random shift:

```
r_shift <- function( Mgroups ){
	set.seed(-1) ## unset the seed
	A = sample( 1:100, 1 )
	B = sample( 1:100, 1 )
	
	if ( A != B ){
		## do sthe shift
		tmp = Mgroups[match(A,Mgroup[,1]),2]
		Mgroups[match(A,Mgroup[,1]),2] = Mgroups[match(B,Mgroup[,1]),2]
		Mgroups[match(B,Mgroup[,1]),2] = tmp
	}
	return ( Mgroups )
}
```

And while we are at it we need a function to get from this Mgroups to the energy value:

```
getEnergy = function( Mgroups, data ){
	return sum(unlist (lapply( split(Mgroups[,1], Mgroups[,2]), function(group){ mean(apply( data[group,], 1, var)) } )))
}
```

Now get us a list of 10000 values?

```
dist = lapply(1:10000, function() { } )
```

Ha and that is where I definitely want this to become an R6 class.
This does not work if I have my data in an S4 class as I would need to change the Mgroups while running this function. And that is not possible with an R S4 class.


```
library(R6)
library(reshape2)

Data <- R6Class("MonteCarloData",
  public = list(
    data = data,
    Mgroups = melt(groups),
    initialize = function() {
    }
  )
)

Data = Data$new()
```

And a short change to to the two functions we had:

```
r_shift <- function( MCdata ){
	set.seed( as.numeric(stringr::str_sub(as.character(round(as.numeric(Sys.time())*10000)), -5,-1)) )
	A = sample( 1:nrow(MCdata$data), 1 )
	B = sample( 1:nrow(MCdata$data), 1 )
	
	if ( A != B ){
		## do sthe shift
		#print( paste("change",A,"to", B, "from", MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2], "to",  MCdata$Mgroups[match(B,Mgroups[,1]),2]) )
		tmp = MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2]
		MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2] = MCdata$Mgroups[match(B,Mgroups[,1]),2]
		MCdata$Mgroups[match(B,MCdata$Mgroups[,1]),2] = tmp
	}
	return ( MCdata )
}

getEnergy <- function( MCdata ){
	return( mean(unlist (lapply( split(MCdata$Mgroups[,1], MCdata$Mgroups[,2]), function(group){ 
		#warning( paste( group ) )
		apply( MCdata$data[group,], 2, var)
	} ))))
}
```

Now we get our dist:

```
dist = unlist( lapply(1:1000, function(id) { getEnergy( r_shift(Data) ) } ) )

hist(dist)
```

## OK this was naive

But at least I cleared out a problem in my getEnergy function ( I had the energy calculates over 1 not over 2 as it should be...).

So reading Shamits Way to run the algo (Fig 4.3 on page 106 in his Thesis.pdf) the r_shift needs a work over! As well as out class.

We need Eold and T.

```
Data <- R6Class("MonteCarloData",
  public = list(
    data = t(cellexalObj@data),
    Mgroups = cbind( colnames(cellexalObj@data), sample(1:10, ncol(cellexalObj@data), replace=TRUE) ),
    T=1,
    c=0.999,
    Eold = +Inf,
    failed=0,
    initialize = function() {
    }
  )
)

Data = Data$new()

drc = cellexalObj@drc[['DDRtree']]
plot(drc, col=rainbow(length(unique(Data$Mgroups[,2])))[as.numeric(Data$Mgroups[order(Data$Mgroups[,1]),2])], main="After the clustering" )

```

And the R_shift needs to implement the check if better.

```
r_shift <- function( MCdata ){
	set.seed( as.numeric(paste( sep="","1", stringr::str_sub(as.character(round(as.numeric(Sys.time())*10000)), -5,-1)) ))
	old = MCdata$Mgroups
	A = sample( rownames(MCdata$data), 1 )
	B = sample( rownames(MCdata$data), 1 )
	
	if ( MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2] != MCdata$Mgroups[match(B,MCdata$Mgroups[,1]),2] ){
		## do the shift
		#print( paste("change",A,"to", B, "from", MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2], "to",  MCdata$Mgroups[match(B,MCdata$Mgroups[,1]),2]) )
		tmp = MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2]
		MCdata$Mgroups[match(A,MCdata$Mgroups[,1]),2] = MCdata$Mgroups[match(B,MCdata$Mgroups[,1]),2]
		MCdata$Mgroups[match(B,MCdata$Mgroups[,1]),2] = tmp
		Enew = getEnergy( MCdata )
		if ( Enew < MCdata$Eold ) {
			MCdata$Eold = Enew
			MCdata$failed = 0
			return ( MCdata )
		}
		else {
			if ( exp( - ( Enew - MCdata$Eold ) / MCdata$T ) > runif(1)  ){
				print ( paste( exp( - ( Enew - MCdata$Eold) / MCdata$T ) ,"at temp", MCdata$T ) )
				#browser()
				#print("changed grouping")
				MCdata$Eold = Enew
				MCdata$failed = 0
			}
			else{
				#print( "failed change")
				MCdata$failed = MCdata$failed +1
				if ( MCdata$failed > 5000){
					stop( "5000 failed tries in a row - optimal?" )
				}
				MCdata$Mgroups = old ## rewind
			}
		}
	}
	MCdata$T = MCdata$T * MCdata$c
	invisible ( MCdata )
}
```

And the updated Function to run this:

```

png( file="before.png", width=800, height=800)
plot(drc, col=rainbow(length(unique(Data$Mgroups[,2])))[as.numeric(Data$Mgroups[order(Data$Mgroups[,1]),2])], main="After the clustering" )
dev.off()

startTime <- Sys.time()

while ( Data$T > .0 ) {
	r_shift( Data )
}

endTime <- Sys.time()

print(endTime - startTime)

png( file="after.png", width=800, height=800)
plot(drc, col=rainbow(length(unique(Data$Mgroups[,2])))[as.numeric(Data$Mgroups[order(Data$Mgroups[,1]),2])], main="After the clustering" )
dev.off()
```

### This takes forever to finish.

This code is horrible slow. OK I also have a print hidden in the function.

Let's tackle the first idea of improvement and calculate the varianze only for the changed clusters.

Change our object to store these values:

```
MonteCarloData <- R6Class("MonteCarloData",
  public = list(
    data = NULL,
    vars_old = NULL,
    vars_new= NULL,
    changed_clusters=NULL,
    Mgroups = NULL,
    T=1,
    c=0.999,
    Eold = +Inf,
    failed=0,
    initialize = function() {
    	self$data = t(cellexalObj@data)
    	self$Mgroups = cbind( colnames(cellexalObj@data), sample(1:10, ncol(cellexalObj@data), replace=TRUE) )
    	self$vars_old = lapply( split(self$Mgroups[,1], self$Mgroups[,2]), function(group){ 
			var( matrix(self$data[group,]))
		} )
		self$changed_clusters = c()
		self$c=0.999
		self$T = 1
		self$Eold = +Inf,
		self$failed=0
    }
  )
)

Data = MonteCarloData$new()

drc = cellexalObj@drc[['DDRtree']]
plot(drc, col=rainbow(length(unique(Data$Mgroups[,2])))[as.numeric(Data$Mgroups[order(Data$Mgroups[,1]),2])], main="After the clustering" )

```

And the getEnergy needs to be changed, too:

So now we are taking it to a next level - look into the R subdirectory. From here on the project will be a real R class.

```
devtools::load_all()
data = MonteCarloData$new()

drc = cellexalvrR::cellexalObj@drc[['DDRtree']]

png( file="before.png", width=800, height=800)
plot(drc, col=rainbow(length(unique(data$Mgroups[,2])))[as.numeric(data$Mgroups[order(data$Mgroups[,1]),2])], main="After the clustering" )
dev.off()

startTime <- Sys.time()

while ( data$T > .0 ) {
	r_shift( data )
}

endTime <- Sys.time()

print(endTime - startTime)

png( file="after.png", width=800, height=800)
plot(drc, col=rainbow(length(unique(Data$Mgroups[,2])))[as.numeric(Data$Mgroups[order(Data$Mgroups[,1]),2])], main="After the clustering" )
dev.off()

```



