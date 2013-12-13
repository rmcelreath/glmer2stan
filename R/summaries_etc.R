
# pretty print estimates
stanmer <- function( fit , digits=2 , probs=c(0.025,0.975) , fixed_prefix="beta_" , vary_prefix="vary_" ) {
    spot <- function(pattern,target) length( grep( pattern , target , fixed=TRUE ) )>0
    unprefix <- function(astring) {
        astring <- gsub( fixed_prefix , "" , astring , fixed=TRUE )
        astring <- gsub( vary_prefix , "" , astring , fixed=TRUE )
        astring <- gsub( "sigma_" , "" , astring , fixed=TRUE )
        astring <- gsub( "Sigma_" , "" , astring , fixed=TRUE )
        astring <- gsub( "Intercept" , "(Intercept)" , astring , fixed=TRUE )
        astring
    }
    unprefix2 <- function(astring) {
        astring <- gsub( fixed_prefix , "" , astring , fixed=TRUE )
        astring <- gsub( "Intercept" , "(Intercept)" , astring , fixed=TRUE )
        astring
    }
    post <- extract( fit , permuted=TRUE )
    post.mu <- list()
    post.se <- list()
    for ( i in 1:length(post) ) {
        dims <- length( dim( post[[i]] ) )
        name <- names(post)[i]
        if ( name!="lp__" ) {
            if ( dims==1 ) {
                post.mu[[ name ]] <- mean( post[[i]] )
                post.se[[ name ]] <- sd( post[[i]] )
            } else {
                post.mu[[ name ]] <- apply( post[[i]] , 2:dims , mean )
                post.se[[ name ]] <- apply( post[[i]] , 2:dims , sd )
            }
        }
    }
    
    fixlist <- c()
    ranlist <- c()
    
    for ( i in 1:length(post.mu) ) {
        name <- names(post.mu)[i]
        
        if ( spot( vary_prefix , name ) ) {
        
            ranlist <- c( ranlist , name )
            
        } else {
        
            vname <- paste( vary_prefix , unprefix(name) , sep="" )
            if ( is.null( post.mu[[vname]] ) & !spot("cutpoints",name) & name!="dev" ) {
                fixlist <- c( fixlist , name )
            }
        
        }
        
    }#i
    
    if ( length(fixlist)>0 ) {
    
        fix <- data.frame( Expectation=as.numeric(post.mu[fixlist]) , StdDev=as.numeric(post.se[fixlist]) )
        rownames(fix) <- sapply( fixlist , unprefix2 )
        
        if ( !is.null(probs) ) {
            q <- matrix( 0 , nrow=nrow(fix) , ncol=length(probs) )
            for ( k in 1:length(fixlist) ) {
                q[k,] <- quantile( post[[ fixlist[k] ]] , probs )
            }
            colnames(q) <- paste( round(probs*100,1) , "%" , sep="" )
            fix <- cbind( fix , q )
        }
        
        cat( paste( "glmer2stan model: " , fit@stanmodel@model_name , "\n\n" , sep="" ) )
        cat( "Level 1 estimates:\n" )
        print( round( fix , digits ) )
    
    } # fixlist
    
    cat( "\nLevel 2 estimates:\n" )
    cat( "(Std.dev. and correlations)\n" )
    
    sigresid <- "NULL"
    
    for ( varyname in ranlist ) {
    
        clustername <- unprefix(varyname)
        
        gini <- attr( fit , "ranef" )[[clustername]]$gini
        
        dims <- length( dim( post.mu[[varyname]] ) )
        
        if ( dims == 2 ) {
        
            signame <- paste( "Sigma_" , clustername , sep="" )
            ngroups <- nrow( post.mu[[ varyname ]] )
            cat( paste("\nGroup: ",clustername," (",ngroups," groups / imbalance: " , round(gini,2) , ")" , sep="" ) )
            
            S <- post.mu[[ signame ]]
            Scorr <- cov2cor(S)
            Sdisplay <- Scorr
            diag(Sdisplay) <- sqrt( diag(S) )
            Sdisplay[ upper.tri(Sdisplay) ] <- NA
            rnames <- attr(fit,"ranef")[[clustername]]$factors
            rnums <- paste( "(" , 1:ncol(Sdisplay) , ")" , sep="" )
            rnames <- paste( rnums , rnames )
            rownames( Sdisplay ) <- rnames
            colnames( Sdisplay ) <- rnums
            
            cat( "\n" )
            print( round( Sdisplay , digits ) )
            
        } else {
        
            signame <- paste( "sigma_" , clustername , sep="" )
            ngroups <- length( post.mu[[ varyname ]] )
            cat( paste("\nGroup: ",clustername," (",ngroups," groups / imbalance: " , round(gini,2) , ")" , sep="" ) )
            cat( paste("\n  (Intercept)  " , round(post.mu[[signame]],digits) , "  (SE " , round(post.se[[signame]],digits) , ")\n" , sep="" ) )
            
        }
            
    }
    
    # DIC
    if ( !is.null( attr(fit,"DIC") ) ) {
        w <- attr(fit,"DIC")
        cat( paste( "\nDIC: " , round(w$DIC,0) , "   pDIC: " , round(w$pD,1) , "   Deviance: " , round(w$Dhat,1) , "\n" , sep="" ) )
    }
    
    # WAIC
    if ( !is.null( attr(fit,"WAIC") ) ) {
        w <- attr(fit,"WAIC")
        cat( paste( "\nWAIC: " , round(w$WAIC,0) , "   pWAIC: " , round(w$pD,1) , "   -2*lppd: " , round(-2*w$lppd,1) , "\n" , sep="" ) )
    }
    
    # return summary stats
    invisible(list( means=post.mu , se=post.se ))
}

stantrace <- function( fit ) {
    rstan::traceplot( fit , ask=TRUE )
}

stanDIC <- function( fit ) {
    attr( fit , "DIC" )
}

stanWAIC <- function( fit ) {
    attr( fit , "WAIC" )
}

gitrstan <- function() {
    message( "Journey to mc-stan.org for a proper introduction to rstan." )
    message( "Attempting to install dependencies..." )
    install.packages('inline')
    install.packages('Rcpp')
    install.packages('RcppEigen')
    message( "Removing previous rstan installs (if any)..." )
    try( remove.packages('rstan') )
    message( "Now trying to download and compile rstan..." )
    options(repos = c(getOption("repos"), rstan = "http://wiki.stan.googlecode.com/git/R"))
    install.packages('rstan', type = 'source')
}

# make a ranef method that reports varying effect expectations in lme4 style
varef <- function( object , ... ) {
        
        if ( is.null( attr(object,"ranef") ) ) {
            stop( "Object does not appear to be fit by glmer2stan." )
        }
        
        unprefix <- function(astring) {
            astring <- gsub( vary_prefix , "" , astring , fixed=TRUE )
            astring
        }
        spot <- function(pattern,target) length( grep( pattern , target , fixed=TRUE ) )>0
        
        vary_prefix <- "vary_"
        
        post <- extract( object , permuted=TRUE )
        post.mu <- list()
        post.se <- list()
        for ( i in 1:length(post) ) {
            dims <- length( dim( post[[i]] ) )
            name <- names(post)[i]
            if ( name!="lp__" & spot(vary_prefix,name) ) {
                # dims check here redundant?
                if ( dims > 1 ) {
                    post.mu[[ unprefix(name) ]] <- as.matrix( apply( post[[i]] , 2:dims , mean ) )
                    post.se[[ unprefix(name) ]] <- as.matrix( apply( post[[i]] , 2:dims , sd ) )
                    # add column names
                    r <- attr(object,"ranef")[[unprefix(name)]]$factors
                    colnames(post.mu[[unprefix(name)]]) <- r
                    colnames(post.se[[unprefix(name)]]) <- r
                }
            }
        }
        
        # result
        list( expectation=post.mu , stddev=post.se )
        
    }#function

setMethod( "show" , signature( object = "list" ),
    function( object ) {
        cat( object$model )
    }
)
