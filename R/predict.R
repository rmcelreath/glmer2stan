########################################################################
# functions to compute posterior predictions from glmer2stan fit models

# to do
# (*) replace missing variables with means from data used in fitting
# (*) when data missing entirely, use data used in fitting
# (*) structure result as data frame?

stanpredict <- function( stanfit , data , vary_prefix="vary_" , fixed_prefix="beta_" , probs=c(0.025,0.975) , nsims=1e4 ) {
    
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring <- gsub( ":" , "_X_" , astring , fixed=TRUE )
        astring <- gsub( "(" , "" , astring , fixed=TRUE )
        astring <- gsub( ")" , "" , astring , fixed=TRUE )
        astring
    }
    
    logistic <- function (x) 
    {
        p <- 1/(1 + exp(-x))
        p <- ifelse(x == Inf, 1, p)
        p
    }
    
    rgamma2 <- function (n, mu, scale) 
    {
        rgamma(n, shape = mu/scale, scale = scale)
    }
    
    rordlogit <- function (n, a, phi = 0) 
    {
        a <- c(as.numeric(a), Inf)
        k <- 1:length(a)
        p <- dordlogit(k, a = a, phi = phi, log = FALSE)
        y <- sample(k, size = n, replace = TRUE, prob = p)
        y
    }

    
    # check params
    
    if ( missing(data) ) {
        stop("no data")
    }
    if ( class(stanfit) != "stanfit" ) {
        stop("requires stanfit object")
    }
    if ( is.null( attr(stanfit,"formulas") ) ) {
        stop("stanfit object not fit by glmer2stan (cannot parse formulas)")
    }
    
    # compute vary and fixef components of linear model
    # need:
    # number of formulas
    # number of cases for each formula
    # ranef structure from parsed formula
    # fixef structure
    # cluster_vars
    
    fp <- attr( stanfit , "formulas" )$fp
    num_formulas <- length(fp)
    cluster_vars <- attr( stanfit , "formulas" )$cluster_vars
    var_suffix <- attr( stanfit , "formulas" )$var_suffix
    family <- attr( stanfit , "formulas" )$family
    formula <- attr( stanfit , "formulas" )$formula
    
    # check cluster vars are class integer
    for ( i in 1:length(cluster_vars) ) {
        varname <- names(cluster_vars)[i]
        if ( !is.null( data[[ varname ]] ) ) {
            data[[ varname ]] <- as.integer( data[[ varname ]] )
        } else {
            # no cluster values in data!
            # so just substitute with 0 => not use random effects in prediction
            data[[ varname ]] <- as.integer(0)
        }
    }
    
    # check for missing outcome values in data
    # don't need outcomes, but design.matrix wants them for its input
    # so can just add placeholders
    for ( f in 1:num_formulas ) {
        if ( is.null(data[[ fp[[f]]$yname ]]) ) 
            data[[ fp[[f]]$yname ]] <- 0
        # for binomial, also need failure count
        # so stop and warn user
        if ( family[[f]]=="binomial" ) {
            if ( length( formula[[f]][[2]] )==3 ) {
                # cbind construction
                bintotname <- deparse( formula[[f]][[2]][[3]] )
                # check for presence
                if ( is.null( data[[ bintotname ]] ) )
                    stop( paste( "Failure count '" , bintotname , "' required as element of data list." , sep="" ) )
            } else {
                # bernoulli construction
                # do nothing for now
            }
        }
    }
    
    # build data (design) matrix by re-parsing formulas with new data frame
    # this is mainly needed so user doesn't have to manually construct interactions
    
    data_list <- as.data.frame(data) # uses repetition to fill all cases
    
    data_temp <- list()
    for ( f in 1:num_formulas ) {
        temp_parse <- parse_formula( formula[[f]] , data_list )
        data_temp[[f]] <- temp_parse$dat
        
        # prep binomial size variable
        if ( family[[f]]=="binomial" ) {
            varname <- paste( "bin_total" , var_suffix[f] , sep="" )
            if ( class( temp_parse$y )=="matrix" ) {
                # cbind() form
                # split matrix into two distinct columns
                bin_tot <- as.integer( temp_parse$y[,1] + temp_parse$y[,2] )
                temp_parse$y <- as.integer( temp_parse$y[,1] ) # extract first column (success count)
                data_list[[ varname ]] <- bin_tot
            } else {
                # bernoulli trials
                data_list[[ varname ]] <- rep(1,nrow(data_list))
            }
        }
    }
    
    # merge dat tables
    for ( f in 1:num_formulas ) {
        for ( j in 1:ncol(data_temp[[f]]) ) {
            acolname <- colnames( data_temp[[f]] )[j]
            data_list[[ acolname ]] <- data_temp[[f]][j][[1]]
        }
    }
    
    # assign( "zz" , data_list , envir=.GlobalEnv )
    
    # compute number of cases from longest vector in data
    N_all <- nrow(data_list)
    # convert column names
    new_names <- sapply( colnames(data_list) , undot )
    colnames(data_list) <- new_names
    
    # extract samples
    post <- extract( stanfit )
    
    # prep empty result
    result <- list()
    
    # loop over formulas and compute
    for ( f in 1:num_formulas ) {
        
        # compute varying component of GLM
        # cases in rows, samples in cols
        vary <- matrix( 0 , nrow=N_all , ncol=length(post$lp__) )
        if ( length( fp[[f]]$ranef ) > 0 ) {
            for ( i in 1:length( fp[[f]]$ranef ) ) {
                cluster_name <- undot( names(fp[[f]]$ranef)[i] )
                parname <- paste( vary_prefix , cluster_name , sep="" )
                nterms <- length( fp[[f]]$ranef[[i]] )
                total_terms <- sum( cluster_vars[[ cluster_name ]] ) # across all formulas
                fnames <- fp[[f]]$ranef[[i]]
                for ( j in 1:nterms ) {
                    for ( k in 1:N_all ) { # loop over cases
                        jstart <- 0
                        if ( f > 1 ) {
                            jstart <- sum( cluster_vars[[ cluster_name ]][ 1:(f-1) ] )
                        }
                        dox <- 1
                        if ( undot(fnames[j]) != "Intercept" ) {
                            xname <- paste( undot(fnames[j]) , var_suffix[f] , sep="" )
                            dox <- data_list[[ xname ]][k]
                        }
                        cname <- paste( cluster_name , var_suffix[f] , sep="" )
                        if ( data_list[[cname]][k] > 0 ) { 
                        # index 0 indicates not to use varying effect in prediction
                            if ( total_terms > 1 ) {
                                vary[k,] <- vary[k,] + post[[ parname ]][ , data_list[[cname]][k] , j+jstart ] * dox
                            } else {
                                vary[k,] <- vary[k,] + post[[ parname ]][ , data_list[[cname]][k] ] * dox
                            }
                        }
                    } #k
                } #j
            } #i
        }
        
        # compute fixed component of GLM
        nterms <- length( fp[[f]]$fixef )
        fnames <- fp[[f]]$fixef
        glm <- matrix( 0 , nrow=N_all , ncol=length(post$lp__) )
        for ( j in 1:nterms ) {
            prefix <- fixed_prefix
            for ( k in 1:N_all ) {
                if ( undot(fnames[j])=="Intercept" ) {
                    if ( family[[f]]!="ordered" ) {
                        parname <- paste( undot(fnames[j]) , var_suffix[f] , sep="" )
                        glm[k,] <- glm[k,] + post[[ parname ]]
                    }
                } else {
                    parname <- paste( prefix , undot(fnames[j]) , var_suffix[f] , sep="" )
                    xname <- paste( undot(fnames[j]) , var_suffix[f] , sep="" )
                    dox <- data_list[[ xname ]][k]
                    glm[k,] <- glm[k,] + post[[ parname ]] * dox
                }
            } #k
        } #j
        
        # join
        glm2 <- vary + glm
        
        # simulate observations
        # uses outcome density
        if ( family[[f]] == "gaussian" ) {
            # gaussian noise around mean
            parname <- paste( "sigma" , var_suffix[f] , sep="" )
            outsim <- sapply( 1:nrow(glm2) , function(i) quantile( 
                rnorm( nsims , glm2[i,] , post[[parname]] ) , probs=probs ) )
        }
        if ( family[[f]] == "binomial" ) {
            # apply inverse link
            glm2 <- logistic( glm2 )
            # simulate
            bintotname <- paste( "bin_total" , var_suffix[f] , sep="" )
            outsim <- sapply( 1:nrow(glm2) , function(i) quantile( 
                rbinom( nsims , prob=glm2[i,] , size=data_list[[bintotname]][i] )/data_list[[bintotname]][i] , probs=probs ) )
        }
        if ( family[[f]] == "ordered" ) {
            message( "Outcome simulations for family 'ordered' not yet implemented" )
            # NYI - rordlogit doesn't vectorize over posterior cutpoints
            cutsname <- paste( "cutpoints" , var_suffix[f] , sep="" )
            outsim <- sapply( 1:nrow(glm2) , function(i) quantile( 
                rordlogit( nsims , phi=glm2[i,] , a=post[[cutsname]] ) , probs=probs ) )
        }
        if ( family[[f]] == "gamma" ) {
            glm2 <- exp( glm2 ) # log link
            parname <- paste( "theta" , var_suffix[f] , sep="" )
            outsim <- sapply( 1:nrow(glm2) , function(i) quantile( 
                rgamma2( nsims , mu=glm2[i,] , scale=1/post[[parname]] ) , probs=probs ) )
        }
        if ( family[[f]] == "poisson" ) {
            glm2 <- exp( glm2 ) # log link
            outsim <- sapply( 1:nrow(glm2) , function(i) quantile( 
                rpois( nsims , lambda=glm2[i,] ) , probs=probs ) )
        }
        
        # store results
        # compute summary stats
        mu <- apply( glm2 , 1 , mean )
        ci <- apply( glm2 , 1 , quantile , probs=probs )
        obs <- outsim
        
        # result
        # name by outcome variable
        result[[ fp[[f]]$yname ]] <- list( mu=mu , mu.ci=ci , obs.ci=obs )
        
    } #f

    result
    
}

# test

# data(UCBadmit)
# m <- glmer2stan( cbind(admit,reject) ~ (1|dept) + male , data=UCBadmit , family="binomial" )

# z <- stanpredict( m , data=UCBadmit )[[1]]
# UCBadmit$p.admit <- UCBadmit$admit / UCBadmit$applications
# plot( UCBadmit$p.admit , ylim=c(0,1) , pch=ifelse(UCBadmit$male==1,16,1) )
# lines( 1:nrow(UCBadmit) , z$mu )
# shade( z$mu.ci , 1:nrow(UCBadmit) )
# shade( z$obs.ci , 1:nrow(UCBadmit) )

# z <- stanpredict( m , data=list( reject=100 , male=c(0,1) , dept=0 ) )[[1]]
# plot( 0:1 , type="n" , xlim=c(0,1) , ylim=c(0,1) , xlab="female/male" , ylab="prob admit" )
# lines( 0:1 , z$mu )
# shade( z$mu.ci , 0:1 )
# shade( z$obs.ci , 0:1 )

# data(chimpanzees)
# m2 <- glmer2stan( pulled.left ~ (1+prosoc.left|actor) + prosoc.left * condition - condition , data=chimpanzees , family="binomial" , varpriors="weak" , sample=TRUE )
# z <- stanpredict( m2 , data=chimpanzees )[[1]]