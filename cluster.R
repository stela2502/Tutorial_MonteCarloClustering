

cluster <- function(data, steps = 1e+6 ){

	i = 0;
	pb <- progressprogress_bar$new( steps )
	for ( i in 1:steps ){
		r_shift( data )
	}

	print "Finished"
	
	invisible ( data )
}