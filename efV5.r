efV5 = function(formula, input_data, model = c("Sph","Exp", "Gau"),
     kappa = 0, fix.values = c(NA,NA,NA),								verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA), 
                                miscFitOptions = list(),...)
# This function automatically fits a variogram to input_data
{

    # Check for anisotropy parameters
 #   if('alpha' %in% names(list(...))) warning('Anisotropic variogram model fitting not supported, see the documentation of efV2 for more details.')
     print("in efV5")  
     m1 <- c("Sph","Exp","Gau")

     n2 <- length(m1)
     yfin <- 9*(10^10)


    # Take the misc fit options and overwrite the defaults by the user specified ones
    miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5)
  
     miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
     

    # Create boundaries
    longlat = !is.projected(input_data)
    if(is.na(longlat)) longlat = FALSE
    diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1,2]                # 0.35 times the length of the central axis through the area
    boundaries = c(2,4,6,9,12,15,25,35,50,65,80,100) * diagonal * 0.35/100         # Boundaries for the bins in km
     

  	# If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
  	if(!is(GLS.model, "variogramModel")) {
  		experimental_variogram = variogram(formula, input_data,boundaries = boundaries, ...)
  	} else {
  		if(verbose) cat("Calculating GLS sample variogram\n")
  		g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
  		experimental_variogram = variogram(g, boundaries = boundaries, ...)
  	}


	
    # request by Jon Skoien
    if(miscFitOptions[["merge.small.bins"]]) {
      if(verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
      while(TRUE) {
          if(length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
          boundaries = boundaries[2:length(boundaries)]			
          if(!is(GLS.model, "variogramModel")) {
              experimental_variogram = variogram(formula, input_data,boundaries = boundaries, ...)
          } else {
              experimental_variogram = variogram(g, boundaries = boundaries, ...)
          }
      }	
}

    

    # set initial values
    if(is.na(start_vals[1])) {  # Nugget
        initial_nugget = min(experimental_variogram$gamma)
    } else {
        initial_nugget = start_vals[1]
    }
    if(is.na(start_vals[2])) { # Range
        initial_range = 0.1 * diagonal   # 0.10 times the length of the central axis through the area
    } else {
        initial_range = start_vals[2]
    }
    if(is.na(start_vals[3])) { # Sill
        initial_sill = mean(c(max(experimental_variogram$gamma), median(experimental_variogram$gamma)))
    } else {
        initial_sill = start_vals[3]
    }
    
    # Determine what should be automatically fitted and what should be fixed
    # Nugget
    if(!is.na(fix.values[1]))
    {
        fit_nugget = FALSE
        initial_nugget = fix.values[1]
    } else
        fit_nugget = TRUE

    # Range
    if(!is.na(fix.values[2]))
    {
        fit_range = FALSE
        initial_range = fix.values[2]
    } else
        fit_range = TRUE

    # Partial sill
    if(!is.na(fix.values[3]))
    {
        fit_sill = FALSE
        initial_sill = fix.values[3]
    } else
        fit_sill = TRUE

     




    size <- 0
    id <- 0
    comm <- 1
    id <- .Fortran("fmpib",as.integer(comm),as.integer(id))[[2]]
    size <- .Fortran("fmpia",as.integer(comm),as.integer(size))[[2]]
    lastnode <- size
    big.one <- vector("list",length=size)
     
     tmp <- floor(n2/lastnode)
     cat("stuff",id,size,"\n")

     mystart <- id * tmp + 1
     myend <- (id+1) * tmp 
     if(id==(lastnode-1))myend <- n2 
     vgm_list <- NULL
     SSerr_list <- NULL
   
   print("setup")
   cat("mo",mystart,myend,"\n")

    # Automatically testing different models, the one with the smallest sums-of-squares is chosen

   n4 <- (myend - mystart + 1)
   mod.list <- as.list(m1[mystart:myend])
   kap.list <- as.list(kappa[mystart:myend])
   model_fit <- vector("list",length=n4)
   
   


for(j in 1:n4) {
      print("loop")
      cat(j,initial_nugget,"\n")
        obj <-  try(fit.variogram(experimental_variogram,
                        model = vgm(psill=initial_sill-initial_nugget,
 model=mod.list[[j]], range=initial_range,
        nugget=initial_nugget,kappa = kap.list[[j]]),
        fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
	debug.level = 0), TRUE)


		if("try-error" %in% class(obj)) {
			#print(traceback())
			warning("An error has occured during variogram fitting. Used:\n", 
					"\tnugget:\t", nugget, 
					"\n\tmodel:\t", model, 
					"\n\tpsill:\t", psill,
					"\n\trange:\t", range,
					"\n\tkappa:\t",ifelse(kappa == 0, NA, kappa),
					"\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj)
		}
		
			if(!is.null(obj)) {	# skip models that failed
				vgm_list[[j]] = obj
				SSerr_list  = c(SSerr_list, attr(obj, "SSErr"))

				}

}


     
# Check for negative values in sill or range coming from fit.variogram
    # and NULL values in vgm_list, and remove those with a warning
    strange_entries <-  sapply(vgm_list, function(v) any(c(v$psill, v$range) < 0) | is.null(v)) 

        if(any(strange_entries)) {
        SSerr_list[strange_entries] = (9*10^10)
	}


	      
	fake <- data.frame(model=factor("Nug","Exp"),psill=c(0,0.2),range=rep(0,2),
	kappa=rep(0,2),ang1=rep(0,2),ang2=rep(0,2),ang3=rep(0,2),
	anis1=rep(1,2),anis2=rep(1,2))
	attr(fake, "SSErr") <- (9*10^10)
	class(fake) <- c("variogramModel","data.frame")


	n3 <- n4 - length(SSerr_list)
	if(length(SSerr_list)!=n4) {
	SSerr_list <- c(SSerr_list,
		rep(9*(10^10),n3))
	}
        

	etest <- min(SSerr_list)
	big.two <- which.min(SSerr_list)
        fmod <- vgm_list[[big.two]]

	if(etest != (9*10^10)) {
	      big.list <- list(vgm=fmod,SSErr=etest)

	      }
	      else {
	      big.list <- list(vgm=fake,SSErr=etest)
	      }
	      
	      print("list")
	      print(big.list)
    big.one[[id+1]] <- big.list
    yfina <- .Fortran("freda",as.double(big.list$SSErr),as.integer(id),as.double(0.00),as.integer(0))

	


	if(id==0) {
    sse <- yfina[[3]]
    irank <- as.integer(yfina[[4]])
    print("final sse")
    print(sse)
    print(irank)
    print(big.one[[irank+1]])
    print(big.one[irank])
	big.fin <- big.one[[irank+1]]
	assign("sse.fin",sse,.GlobalEnv)
	       
	return(big.fin)       
	}


}
