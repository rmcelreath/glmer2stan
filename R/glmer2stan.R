# to do:
# (*) families: use built-in families and links (the closures) instead of strings
# (*) turbo option, to apply matt trick and use matrix operations for glm (design matrix)...reduces readability of model, but potentially large speed benefits
# (*) eventually make default chains=4, but need to monitor Rhat
# (*) make "conjugate" priors for traditional inv_wishart/inv_gamma...these priors suck, but people like them
# (*) check ordered response to ensure starts at 1
# (*) allow passing specific init values for specific named parameters - can just pass a named list and then merge it with the init list?
# (*) allow passing fit lme4 object
# (*) allow passing list produced when sample=FALSE 
# (1) add checks for variable types for zigamma...ensure zero outcome is integer
# (2) add checks for cluster variables to make into integer type
# (5) make DIC computation entirely post-sampling? would need to save parsed formulas.

parse_formula <- function( formula , data ) {
    ## take a formula and parse into fixed effect and varying effect lists
    nobars <- function(term) {
        if (!('|' %in% all.names(term))) return(term)
        if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
        if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) return(NULL)
        term[[2]] <- nb
        return(term)
        }
        nb2 <- nobars(term[[2]])
        nb3 <- nobars(term[[3]])
        if (is.null(nb2)) return(nb3)
        if (is.null(nb3)) return(nb2)
        term[[2]] <- nb2
        term[[3]] <- nb3
        term
    }
    findbars <- function(term) {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
        if (!is.call(term)) stop("term must be of class call")
        if (term[[1]] == as.name('|')) return(term)
        if (length(term) == 2) return(findbars(term[[2]]))
        c(findbars(term[[2]]), findbars(term[[3]]))
    }
    subbars <- function(term)
    ### Substitute the '+' function for the '|' function
    {
        if (is.name(term) || !is.language(term)) return(term)
        if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
        }
        stopifnot(length(term) >= 3)
        if (is.call(term) && term[[1]] == as.name('|'))
        term[[1]] <- as.name('+')
        for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
        term
    }
    hasintercept <- function(term) {
        attr( terms(term) , "intercept" )==1
    }
    
    # find fixed effects list by deleting random effects and expanding
    f_nobars <- nobars( formula )
    # catch implied intercept error -- happens when right side of formula is only () blocks
    if ( class(f_nobars)=="name" & length(f_nobars)==1 ) {
        f_nobars <- nobars( as.formula( paste( deparse(formula) , "+ 1" ) ) )
    }
    fixef <- colnames( model.matrix( f_nobars , data ) )
    
    # convert to all fixed effects and build needed model matrix
    mdat <- model.matrix( subbars( formula ) , data )
    outcome_name <- deparse( f_nobars[[2]] )
    # mdat <- cbind( data[[outcome_name]] , mdat )
    # colnames(mdat)[1] <- outcome_name
    outcome <- model.frame( f_nobars , data )[,1]
    if ( class(outcome)=="matrix" ) {
        # fix outcome name
        outcome_name <- colnames( outcome )[1]
    }
    
    # check for any varying effects
    if ( formula == nobars(formula) ) {
        # no varying effects
        ranef <- list()
    } else {
        # find varying effects list
        var <- findbars( formula )
        ranef <- list()
        for ( i in 1:length(var) ) {
            name <- all.vars( var[[i]][[3]] )
            
            # check that grouping vars are class integer
            if ( class( data[[name]] )!="integer" ) {
                stop( paste( "Grouping variables must be integer type. '" , name , "' is instead of type: " , class( data[[name]] ) , "." , sep="" ) )
            }
            
            # parse formula
            v <- var[[i]][[2]]
            if ( class(v)=="numeric" ) {
                # just intercept
                ranef[[ name ]] <- "(Intercept)"
            } else {
                # should be class "call"
                # need to convert formula to model matrix headers
                f <- as.formula( paste( "~" , deparse( v ) , sep="" ) )
                ranef[[ name ]] <- colnames( model.matrix( f , data ) )
            }
        }
    }
    
    # result sufficient information to build Stan model code
    list( y=outcome , yname=outcome_name , fixef=fixef , ranef=ranef , dat=as.data.frame(mdat) )
}

glmer2stan <- function( formula , data , family="gaussian" , varpriors="flat" , sample=TRUE , warmup=5000 , iter=10000 , chains=1 , initmethod="zero" , extract=FALSE , calcDIC=FALSE , verbose=TRUE , fixed_prefix="beta_" , vary_prefix="vary_" , ... ) {
    
    # hard-coded options
    indent <- "    " # 4 spaces
    REML <- TRUE # only used when calling lmer to get inits
    
    # table of vectorized densities
    vectorized <- c(1,1,1,1,0)
    names(vectorized) <- c( "normal" , "binomial" , "poisson" , "gamma" , "ordered" )
    
    # private functions
    class2type <- function( col ) {
        classes <- c("integer","numeric")
        types <- c("int","real")
        types[ which( classes==class(col) ) ]
    }
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring <- gsub( ":" , "_X_" , astring , fixed=TRUE )
        astring <- gsub( "(" , "" , astring , fixed=TRUE )
        astring <- gsub( ")" , "" , astring , fixed=TRUE )
        astring
    }
    make_vary_text <- function( fp , vary_prefix="vary_" , suffix="" , f_num ) {
        vterms <- c()
        for ( i in 1:length(fp$ranef) ) {
            name <- undot( names(fp$ranef)[i] )
            factors <- fp$ranef[[i]]
            jstart <- 0
            if ( f_num > 1 ) jstart <- sum( cluster_vars[[name]][ 1:(f_num-1) ] )
            total_factors <- sum( cluster_vars[[ name ]] )
            for ( j in 1:length(factors) ) {
                if ( length(factors)==1 & total_factors==1 ) {
                    # only varying intercept, and same cluster not used in another formula
                    aterm <- paste( vary_prefix , name , "[" , name , suffix , "[i]]" , sep="" )
                } else {
                    # need effect index
                    aterm <- paste( vary_prefix , name , "[" , name , suffix , "[i]," , j+jstart , "]" , sep="" )
                }
                if ( j > 1 ) {
                    aterm <- paste( aterm , " * " , undot(factors[j]) , suffix , "[i]" , sep="" )
                }
                vterms <- c( vterms , aterm )
            }
        }
        vterms
    }
    make_fixed_text <- function( fp , fixed_prefix="beta_" , suffix="" , drop_intercept=FALSE ) {
        fterms <- c()
        for ( i in 1:length(fp$fixef) ) {
            name <- undot( fp$fixef[i] )
            if ( name=="Intercept" ) {
                if ( !drop_intercept ) {
                    aterm <- paste( "Intercept" , suffix , sep="" )
                    fterms <- c( fterms , aterm )
                }
            } else {
                aterm <- paste( fixed_prefix , name , suffix , " * " , name , suffix , "[i]" , sep="" )
                fterms <- c( fterms , aterm )
            }
        }
        fterms
    }
    logistic <- function(x) {
        p <- 1/(1 + exp(-x))
        p <- ifelse(x == Inf, 1, p)
        p
    }
    logit <- function(x) {
        log( x / (1-x) )
    }
    dordlogit <- function(x, a, phi, log = FALSE) {
        a <- c(as.numeric(a), Inf)
        p <- logistic(a[x] - phi)
        na <- c(-Inf, a)
        np <- logistic(na[x] - phi)
        p <- p - np
        if (log == TRUE) p <- log(p)
        p
    }
    dgamma2 <- function(x, mu, scale = NULL, log = FALSE) {
        dgamma(x, shape = mu/scale, scale = scale, log = log)
    }
        
    # check params
    if ( warmup >= iter ) {
        stop( "Number of iterations (iter) must exceed number of warmup steps (warmup)." )
    }
    
    if ( class(family)!="list" ) {
        family <- list(family)
    }
    for ( f in 1:length(family) ) {
        legal_families <- c( "gaussian" , "binomial" , "poisson" , "ordered" , "gamma" , "zigamma" , "zipoisson" )
        if ( !(family[[f]] %in% legal_families ) ) {
            stop( paste( "Family '" , family , "' not recognized." , sep="" ) )
        }
    }
    
    if ( !(initmethod %in% c("lme4","random","zero")) ) {
        stop( paste( "Init method '" , initmethod , "' not recognized." , sep="" ) )
    }
    if ( !(varpriors %in% c("weak","flat")) ) {
        stop( paste( "Variance priors '" , varpriors , "' not recognized." , sep="" ) )
    }
    if ( missing(formula) || missing(data) ) {
        stop( "Must provide both formula and data objects." )
    }
    
    fixed_prefix <- undot( fixed_prefix )
    vary_prefix <- undot( vary_prefix )
    
    if ( sample==TRUE )
        require(rstan)
    
    ########################################################################
    # parse formulas
    if ( class(formula)!="list" ) {
        formula <- list(formula)
    }
    
    # automatic splitting of outcome for zero-inflated families
    if ( family[[1]]=="zigamma" | family[[1]]=="zipoisson" ) {
        if ( length(formula)==1 ) {
            # only one formula, so split outcome variable and duplicate formula to make mixture
            outvar <- deparse( formula[[1]][[2]] )
            y <- data[[ outvar ]]
            y1 <- ifelse( y==0 , 1 , 0 ) # zeros
            if ( family[[1]]=="zigamma" ) {
                y2 <- ifelse( y==0 , NA , y ) # non-zeros
            }
            if ( family[[1]]=="zipoisson" ) {
                y2 <- as.integer(y) # all counts; i.e. true "inflation"
            }
            zname <- paste( outvar , "_zeros" , sep="" )
            nzname <- paste( outvar , "_nonzero" , sep="" )
            if ( family[[1]] == "zipoisson" ) {
                nzname <- paste( outvar , "_count" , sep="" )
            }
            data[[ zname ]] <- as.integer( y1 )
            data[[ nzname ]] <- y2
            f1 <- as.formula( paste( zname , " ~ " , deparse( formula[[1]][[3]] ) , sep="" ) )
            f2 <- as.formula( paste( nzname , " ~ " , deparse( formula[[1]][[3]] ) , sep="" ) )
            formula <- list( f1 , f2 )
            if ( family[[1]] == "zigamma" ) {
                family <- list( "binomial" , "gamma" )
            }
            if ( family[[1]] == "zipoisson" ) {
                family <- list( "binomial" , "poisson" )
            }
        }
    }
    
    if ( family[[1]]=="zigamma" & length(formula)==2 ) {
        # two formulas, but need to split family into list
        family <- list( "binomial" , "gamma" )
    }
    
    if ( family[[1]]=="zipoisson" & length(formula)==2 ) {
        # two formulas, but need to split family into list
        family <- list( "binomial" , "poisson" )
    }
    
    # parse formulas
    fp <- list()
    for ( f in 1:length(formula) ) {
        fp[[f]] <- parse_formula( formula[[f]] , data )
    }
    num_formulas <- length(formula)
    
    # make suffix labels for each formula
    var_suffix <- ""
    if ( num_formulas > 1 ) {
        var_suffix <- paste( "_" , 1:num_formulas , sep="" )
    }
    if ( length(family)==2 ) {
        if ( family[[1]]=="binomial" & family[[2]]=="gamma" ) {
            var_suffix <- c( "_pi" , "_mu" )
        }
        if ( family[[1]]=="binomial" & family[[2]]=="poisson" ) {
            var_suffix <- c( "_z" , "_n" )
        }
    }
    
    ########################################################################
    # fit with lme4
    fit <- list()
    if ( initmethod=="lme4" ) {
    
        require(lme4)
        
        if ( verbose==TRUE ) {
            message( "Fitting lme4 model" )
            flush.console()
        }
        
        for ( f in 1:num_formulas ) {
        
            use.family <- family[[f]]
            
            if ( family[[f]]=="ordered" ) {
                use.family <- "gaussian"
            }
            if ( family[[f]]=="gamma" ) {
                use.family <- gaussian(link="log")
            }
            
            fit[[f]] <- glmer( formula[[f]] , data=data , family=use.family )
            
        }#f
    }# initmethod lme4
    
    ########################################################################
    # prep model specification
    m_data <- "data{\n"
    # add sample size variables
    for ( i in 1:num_formulas ) {
        m_data <- paste( m_data , indent , "int N" , var_suffix[i] , ";\n" , sep="" )
    }
    
    m_transdata1 <- "transformed data{\n"
    m_transdata2 <- ""
    flag_transdata <- FALSE
    
    m_pars <- "parameters{\n"
    
    m_transpars1 <- "transformed parameters{\n"
    m_transpars2 <- ""
    flag_transpars <- FALSE
    
    # m_model <- "model{\n    real vary;\n    real glm;\n"
    # m_model <- "model{\n    real vary[N];\n    real glm[N];\n"
    # need one vary[] and glm[] variable for each formula, because of different lengths and vectorizing the densities
    m_model <- "model{\n"
    for ( i in 1:num_formulas ) {
        m_model <- paste( m_model , indent , "real vary" , var_suffix[i] , "[N" , var_suffix[i] , "];\n" , sep="" )
        m_model <- paste( m_model , indent , "real glm" , var_suffix[i] , "[N" , var_suffix[i] , "];\n" , sep="" )
    }
    m_gq <- ""
    flag_intercept <- FALSE
    
    ########################################################################
    # data - add columns of data frame (with pre-multiplied interactions)
        
        # various housekeeping
        for ( f in 1:num_formulas ) {
        
            # prep binomial cbind outcome by splitting into success count and total
            if ( family[[f]]=="binomial" ) {
                if ( class(fp[[f]]$y)=="matrix" ) {
                    vname <- fp[[f]]$yname
                    bin_tot <- as.integer( fp[[f]]$y[,1] + fp[[f]]$y[,2] )
                    fp[[f]]$y <- as.integer( fp[[f]]$y[,1] ) # extract first column (success count)
                    fp[[f]]$dat[ , "bin_total" ] <- bin_tot
                } else {
                    # ensure int type
                    fp[[f]]$y <- as.integer( fp[[f]]$y )
                    # make bin_total column full of ones
                    fp[[f]]$dat[ , "bin_total" ] <- as.integer( rep( 1 , nrow(fp[[f]]$dat) ) )
                }
            }
            if ( family[[f]]=="poisson" || family[[f]]=="ordered" ) {
                # ensure outcome variable is integer
                fp[[f]]$y <- as.integer( fp[[f]]$y )
            }
            if ( family[[f]]=="ordered" ) {
                # must add data item for max integer in outcome var
                Kname <- paste( "K" , var_suffix[f] , sep="" )
                m_data <- paste( m_data , indent , "int<lower=2> " , Kname , ";\n" , sep="" )
            }
            # make sure varying effect cluster (index) variables are type integer
            if ( length( fp[[f]]$ranef ) > 0 ) {
                for ( i in 1:length( fp[[f]]$ranef ) ) {
                    fname <- names( fp[[f]]$ranef )[i]
                    fp[[f]]$dat[, fname ] <- as.integer( fp[[f]]$dat[ , fname ] )
                }
            }
            
        }#f
        
        # now add variable names to data block
        for ( i in 1:num_formulas ) {
            
            # add outcome variable to data block
            outvar <- undot( fp[[i]]$yname )
            vtype <- class2type( fp[[i]]$y )
            nname <- paste( "N" , var_suffix[i] , sep="" )
            m_data <- paste( m_data , indent , vtype , " " , outvar , var_suffix[i] , "[" , nname , "];\n" , sep="" )
            
            # other vars
            for ( j in 1:ncol(fp[[i]]$dat) ) {
                vtype <- class2type( fp[[i]]$dat[,j] )
                vname <- undot( colnames(fp[[i]]$dat)[j] )
                if ( vname!="Intercept" ) {
                    nname <- paste( "N" , var_suffix[i] , sep="" )
                    m_data <- paste( m_data , indent , vtype , " " , vname , var_suffix[i] , "[" , nname , "];\n" , sep="" )
                }
            }
        }#i
        
        # for each varying effect cluster, add sample size
        # these sizes span all formulas, because varying effects are sampled outside glm loops
        # first, catalog all cluster variables, across formulas
        # build a list of numbers of factors for each formula, so can later get indices correct
        cluster_vars <- list()
        cluster_size <- list()
        for ( i in 1:num_formulas ) {
            if ( length( fp[[i]]$ranef ) > 0 ) {
                for ( j in 1:length( fp[[i]]$ranef ) ) {
                    name <- names( fp[[i]]$ranef )[j]
                    cluster_name <- undot( name )
                    num_factors <- length( fp[[i]]$ranef[[j]] )
                    num_groups <- length( unique( fp[[i]]$dat[,name] ) )
                    cluster_vars[[ cluster_name ]] <- c( cluster_vars[[ cluster_name ]] , num_factors )
                    cluster_size[[ cluster_name ]] <- max( cluster_size[[ cluster_name ]] , num_groups )
                }
            }
        }
        # now add data variables for number of groups (individuals, households, etc.) in each cluster variable
        if ( length(cluster_vars) > 0 ) {
            for ( i in 1:length(cluster_vars) ) {
                var_name <- names( cluster_vars )[i]
                m_data <- paste( m_data , indent , "int N_" , var_name , ";\n" , sep="" )
            }
        }
        
    ########################################################################
    # parameters
        pars_list <- ""
        # add fixed effects
        for ( f in 1:num_formulas ) {
            for ( i in 1:length( fp[[f]]$fixef ) ) {
                var_name <- undot( fp[[f]]$fixef[i] )
                if ( var_name!="Intercept" ) {
                    var_name <- paste( fixed_prefix , var_name , var_suffix[f] , sep="" )
                    m_pars <- paste( m_pars , indent , "real " , var_name , ";\n" , sep="" )
                    pars_list <- c( pars_list , var_name )
                } else {
                    if ( family[[f]] != "ordered" ) { 
                    # ordered has no intercept, as cutpoints are the intercepts
                        flag_intercept <- TRUE
                        var_name <- paste( var_name , var_suffix[f] , sep="" )
                        m_pars <- paste( m_pars , indent , "real " , var_name , ";\n" , sep="" )
                        pars_list <- c( pars_list , var_name )
                    }
                }
            }
        
            # top level standard deviation for normal outcome variable
            if ( family[[f]]=="gaussian" ) {
                signame <- paste( "sigma" , var_suffix[f] , sep="" )
                m_pars <- paste( m_pars , indent , "real<lower=0> " , signame , ";\n" , sep="" )
                pars_list <- c( pars_list , signame )
            }
            # theta (rate) for gamma -- should rename this "lambda"?
            if ( family[[f]]=="gamma" ) {
                thetaname <- paste( "theta" , var_suffix[f] , sep="" )
                m_pars <- paste( m_pars , indent , "real<lower=0.001> " , thetaname , ";\n" , sep="" )
                pars_list <- c( pars_list , thetaname )
            }
            # add cut points vector
            if ( family[[f]]=="ordered" ) {
                cutsname <- paste( "cutpoints" , var_suffix[f] , sep="" )
                Kname <- paste( "K" , var_suffix[f] , sep="" )
                m_pars <- paste( m_pars , indent , "ordered[" , Kname , "-1] " , cutsname , ";\n" , sep="" )
                pars_list <- c( pars_list , cutsname )
            }
        }#f
        
        # add varying effects and their var-cov matrices
        # since varying effects span formulas, use the cluster_vars list from earlier
        if ( length(cluster_vars) > 0 ) {
            for ( i in 1:length(cluster_vars) ) {
                cluster_name <- names(cluster_vars)[i]
                num_effects <- sum( cluster_vars[[i]] )
                num_groups <- cluster_size[[i]]
                var_name <- paste( vary_prefix , cluster_name , sep="" )
                var_name <- paste( var_name , "[N_" , cluster_name , "]" , sep="" )
                ktype <- "real "
                if ( num_effects > 1 ) ktype <- paste( "vector[" , num_effects , "] " , sep="" )
                # add varying effect vector/matrix
                m_pars <- paste( m_pars , indent , ktype , var_name , ";\n" , sep="" )
                pars_list <- c( pars_list , paste( vary_prefix , cluster_name , sep="" ) )
                # add var-cov
                if ( num_effects == 1 ) {
                    # add standard deviation parameter
                    m_pars <- paste( m_pars , indent , "real<lower=0> sigma_" , cluster_name , ";\n" , sep="" )
                    pars_list <- c( pars_list , paste( "sigma_" , cluster_name , sep="" ) )
                } else {
                    # add var-cov matrix
                    if ( varpriors != "weak" ) {
                        m_pars <- paste( m_pars , indent , "cov_matrix[" , num_effects , "] Sigma_" , cluster_name , ";\n" , sep="" )
                        pars_list <- c( pars_list , paste( "Sigma_" , cluster_name , sep="" ) )
                    } else {
                        # weak priors, so Sigma_ goes in transformed parameters, and sigma_ and Rho_ go in parameters
                        m_pars <- paste( m_pars , indent , "vector<lower=0>[" , num_effects , "] sigma_" , cluster_name , ";\n" , sep="" )
                        m_pars <- paste( m_pars , indent , "corr_matrix[" , num_effects , "] Rho_" , cluster_name , ";\n" , sep="" )
                        m_transpars1 <- paste( m_transpars1 , indent , "cov_matrix[" , num_effects , "] Sigma_" , cluster_name , ";\n" , sep="" )
                        m_transpars2 <- paste( m_transpars2 , indent , "Sigma_" , cluster_name , " <- diag_matrix(sigma_" , cluster_name , ") * Rho_" , cluster_name , " * diag_matrix(sigma_" , cluster_name , ");\n" , sep="" )
                        flag_transpars <- TRUE
                        pars_list <- c( pars_list , paste( "Sigma_" , cluster_name , sep="" ) ) # monitor transformed var-cov matrix
                    }
                }
            }#i
        }
        
        # close parameters block
        m_pars <- paste( m_pars , "}\n" , sep="" )
        
    ############################################
    # model block
        
        m_model <- paste( m_model , indent , "// Priors\n" , sep="" )
        
        # add fixed effect priors
        for ( f in 1:num_formulas ) {
            for ( i in 1:length( fp[[f]]$fixef ) ) {
                name <- undot( fp[[f]]$fixef[i] )
                prefix <- fixed_prefix
                if ( name=="Intercept" ) {
                    if ( family[[f]] != "ordered" ) {
                        m_model <- paste( m_model , indent , name , var_suffix[f] , " ~ normal( 0 , 100 );\n" , sep="" )
                    }
                } else {
                    m_model <- paste( m_model , indent , prefix , name , var_suffix[f] , " ~ normal( 0 , 100 );\n" , sep="" )
                }
            }
        }#f
        
        # add variance priors
        if ( length(cluster_vars) > 0 ) {
            for ( i in 1:length(cluster_vars) ) {
                num_effects <- sum( cluster_vars[[i]] )
                name <- undot( names( cluster_vars )[i] )
                if ( num_effects == 1 ) {
                    if ( varpriors=="weak" )
                        m_model <- paste( m_model , indent , "sigma_" , name , " ~ gamma( 2 , 1e-4 );\n" , sep="" )
                    if ( varpriors=="flat" )
                        m_model <- paste( m_model , indent , "sigma_" , name , " ~ uniform( 0 , 100 );\n" , sep="" )
                } else {
                    if ( varpriors=="weak" ) {
                        # gamma priors for sd's
                        m_model <- paste( m_model , indent , "sigma_" , name , " ~ gamma( 2 , 1e-4 );\n" , sep="" )
                        # m_model <- paste( m_model , indent , "Sigma_" , name , " ~ inv_wishart( " , num_effects+1 , " , Omega_" , name , " );\n" , sep="" )
                        # lkj_corr eta parameter is 1 for flat prior. 1.5 is a slightly humped prior, excluding very large correlations
                        m_model <- paste( m_model , indent , "Rho_" , name , " ~ lkj_corr( 1.5 );\n" , sep="" )
                    }
                    if ( varpriors=="conjugate" ) {
                        m_model <- paste( m_model , indent , "Sigma_" , name , " ~ inv_wishart( " , num_effects+1 , " , Omega_" , name , " );\n" , sep="" )
                        # add Omega matrix to data block
                        m_data <- paste( m_data , indent , "cov_matrix[" , num_effects , "] Omega_" , name , ";\n" , sep="" )
                    }
                    #if ( varpriors=="flat" ) {
                        # do nothing
                    #}
                }
            }#i
        }
        
        # top-level priors
        for ( f in 1:num_formulas ) {
        
            # prior for top-level standard deviation for gaussian outcome
            if ( family[[f]]=="gaussian" ) {
                signame <- paste( "sigma" , var_suffix[f] , sep="" )
                if ( varpriors=="weak" )
                    m_model <- paste( m_model , indent , signame , " ~ gamma( 2 , 1e-4 );\n" , sep="" )
                if ( varpriors=="flat" )
                    m_model <- paste( m_model , indent , signame , " ~ uniform( 0 , 100 );\n" , sep="" )
            }
            # prior for theta in gamma
            if ( family[[f]]=="gamma" ) {
                thetaname <- paste( "theta" , var_suffix[f] , sep="" )
                m_model <- paste( m_model , indent , thetaname , " ~ uniform( 0.001 , 20 );\n" , sep="" )
                # m_model <- paste( m_model , "theta ~ gamma( 2 , 1e-4 );\n " , sep="" )
            }
        
        }#f
        
        # add varying effects sampling
        if ( length(cluster_vars) > 0 ) {
            m_model <- paste( m_model , indent , "// Varying effects\n" , sep="" )
            for ( i in 1:length(cluster_vars) ) {
                cluster_name <- undot( names( cluster_vars )[i] )
                num_effects <- sum( cluster_vars[[i]] )
                if ( num_effects == 1 ) {
                    m_model <- paste( m_model , indent , "for ( j in 1:N_" , cluster_name , " ) " , vary_prefix , cluster_name , "[j] ~ normal( 0 , sigma_" , cluster_name , " );\n" , sep="" )
                } else {
                    m_model <- paste( m_model , indent , "for ( j in 1:N_" , cluster_name , " ) " , vary_prefix , cluster_name , "[j] ~ multi_normal( zeros_" , cluster_name , " , Sigma_" , cluster_name , " );\n" , sep="" )
                    # also need to add zeros vector to transformed data block
                    m_transdata1 <- paste( m_transdata1 , indent , "vector[" , num_effects , "] zeros_" , cluster_name , ";\n" , sep="" )
                    m_transdata2 <- paste( m_transdata2 , indent , "for ( i in 1:" , num_effects , " ) zeros_" , cluster_name , "[i] <- 0;\n" , sep="" )
                    flag_transdata <- TRUE # mark modification of transformed data block
                }
            }
        }
        
        # now loop over each formula's N and sample fixed effects
        m_model <- paste( m_model , indent , "// Fixed effects\n" , sep="" )
        vary_text <- list()
        fixed_text <- list()
        for ( f in 1:num_formulas ) {
        
            # pretty formatting of linear models
            formula_collapse_text <- "\n                + "
            
            # compute varying component of GLM
            if ( length(fp[[f]]$ranef) > 0 ) {
                vary_terms <- make_vary_text( fp[[f]] , vary_prefix , var_suffix[f] , f )
                vary_text[[f]] <- paste( vary_terms , collapse=formula_collapse_text )
                vary_text[[f]] <- paste( indent , indent , "vary" , var_suffix[f] , "[i] <- " , vary_text[[f]] , ";\n" , sep="" )
            } else {
                # no varying effects for this formula
                vary_text[[f]] <- ""
            }
            
            # compute fixed component of GLM
            fixed_text[[f]] <- make_fixed_text( fp[[f]] , fixed_prefix , var_suffix[f] , drop_intercept=family[[f]]=="ordered" )
            fixed_text[[f]] <- paste( fixed_text[[f]] , collapse=formula_collapse_text )
            if ( length(fp[[f]]$ranef) > 0 ) {
                fixed_text[[f]] <- paste( indent , indent , "glm" , var_suffix[f] , "[i] <- vary" , var_suffix[f] , "[i] + " , fixed_text[[f]] , ";\n" , sep="" )
            } else {
                fixed_text[[f]] <- paste( indent , indent , "glm" , var_suffix[f] , "[i] <- " , fixed_text[[f]] , ";\n" , sep="" )
            }
            
            # put vary text and fixed text in the N for loop, but don't close loop
            # to do: recode glm to glm[N] to exploit vectorization of normal density
            m_model <- paste( m_model , indent , "for ( i in 1:N" , var_suffix[f] , " ) {\n" , vary_text[[f]] , fixed_text[[f]] , sep="" )
            
            # finally add top level distribution
            # there must be a better way to do this...
            # should make a table of families and links, so easy to edit
            out_var <- undot( fp[[f]]$yname )
            out_var <- paste( out_var , var_suffix[f] , sep="" )
            
            if ( family[[f]]=="binomial" ) {
                bintotname <- paste( "bin_total" , var_suffix[f] , sep="" )
                if ( FALSE ) {
                # having overflow issues in gen quant block with bernoulli
                # so FALSE'd out for now
                # fix later so can exploit vectorized version
                # although, vectorized code less readable for novices
                    # use bernoulli
                    m_model <- paste( m_model , "}\n" , sep="" )
                    m_model <- paste( m_model , indent , indent , out_var , " ~ bernoulli_logit( glm );" , sep="" )
                } else {
                    # use vectorized binomial
                    m_model <- paste( m_model , indent , indent , "glm" , var_suffix[f] , "[i] <- inv_logit( glm" , var_suffix[f] , "[i] );\n" , sep="" )
                    m_model <- paste( m_model , indent , "}\n" , sep="" )
                    m_model <- paste( m_model , indent , out_var , " ~ binomial( " , bintotname , " , glm" , var_suffix[f] , " );\n" , sep="" )
                }
            }
            
            if ( family[[f]]=="gaussian" ) {
                # vectorized normal
                signame <- paste( "sigma" , var_suffix[f] , sep="" )
                m_model <- paste( m_model , indent , "}\n" , sep="" )
                m_model <- paste( m_model , indent , out_var , " ~ normal( glm" , var_suffix[f] , " , " , signame , " );\n" , sep="" )
            }
            
            if ( family[[f]]=="poisson" ) {
                # vectorized poisson
                m_model <- paste( m_model , indent , indent , "glm[i] <- exp( glm" , var_suffix[f] , "[i] );\n" , sep="" )
                m_model <- paste( m_model , indent , "}\n" , sep="" )
                m_model <- paste( m_model , indent , out_var , " ~ poisson( glm" , var_suffix[f] , " );\n" , sep="" )
            }
            
            if ( family[[f]]=="ordered" ) {
                # can't be vectorized yet?
                cutsname <- paste( "cutpoints" , var_suffix[f] , sep="" )
                m_model <- paste( m_model , indent , indent , out_var , "[i] ~ ordered_logistic( glm" , var_suffix[f] , "[i] , " , cutsname , " );\n    }\n" , sep="" )
            }
            
            if ( family[[f]]=="gamma" ) {
                # vectorized gamma
                thetaname <- paste( "theta" , var_suffix[f] , sep="" )
                m_model <- paste( m_model , indent , indent , "glm" , var_suffix[f] , "[i] <- exp( glm" , var_suffix[f] , "[i] )*" , thetaname , ";\n" , sep="" )
                m_model <- paste( m_model , indent , "}\n" , sep="" )
                m_model <- paste( m_model , indent , out_var , " ~ gamma( glm" , var_suffix[f] , " , " , thetaname , " );\n" , sep="" )
            }
            
        }#f
        
        # close model block
        m_model <- paste( m_model , "}\n" , sep="" )
    
    ###########################################################################
    # build generated quantities, for deviance and DIC calculation
    # should consolidate a lot of this code into previous block to reduce code duplication
    
        if ( calcDIC == TRUE ) {
            m_gq <- "generated quantities{\n    real dev;\n"
            for ( i in 1:num_formulas ) {
                m_gq <- paste( m_gq , indent , "real vary" , var_suffix[i] , "[N" , var_suffix[i] , "];\n" , sep="" )
                m_gq <- paste( m_gq , indent , "real glm" , var_suffix[i] , "[N" , var_suffix[i] , "];\n" , sep="" )
            }
            m_gq <- paste( m_gq , indent , "dev <- 0;\n" , sep="" )
            for ( f in 1:num_formulas ) {
                m_gq <- paste( m_gq , indent , "for ( i in 1:N" , var_suffix[f] , " ) {\n" , sep="" )
                m_gq <- paste( m_gq , vary_text[[f]] , fixed_text[[f]] , sep="" )
                
                # add dev update
                out_var <- undot( fp[[f]]$yname )
                out_var <- paste( out_var , var_suffix[f] , sep="" )
                
                if ( family[[f]]=="binomial" ) {
                    bintotname <- paste( "bin_total" , var_suffix[f] , sep="" )
                    if ( FALSE ) {
                        # use bernoulli - having numerical issues with this function
                        m_gq <- paste( m_gq , indent , indent , "dev <- dev + (-2) * bernoulli_log( " , out_var , "[i] , inv_logit(glm[i]) );\n" , sep="" )
                    } else {
                        # use binomial
                        m_gq <- paste( m_gq , indent , indent , "dev <- dev + (-2) * binomial_log( " , out_var , "[i] , " , bintotname , "[i] , inv_logit(glm" , var_suffix[f] , "[i]) );\n" , sep="" )
                    }
                }
                
                if ( family[[f]]=="gaussian" ) {
                    signame <- paste( "sigma" , var_suffix[f] , sep="" )
                    m_gq <- paste( m_gq , indent , indent , "dev <- dev + (-2) * normal_log( " , out_var , "[i] , glm" , var_suffix[f] , "[i] , " , signame , " );\n" , sep="" )
                }
                
                if ( family[[f]]=="poisson" ) {
                    m_gq <- paste( m_gq , indent , indent , "dev <- dev + (-2) * poisson_log( " , out_var , "[i] , exp(glm" , var_suffix[f] , "[i]) );\n" , sep="" )
                }
                
                if ( family[[f]]=="ordered" ) {
                    cutsname <- paste( "cutpoints" , var_suffix[f] , sep="" )
                    m_gq <- paste( m_gq , indent , indent , "dev <- dev + (-2) * ordered_logistic_log( " , out_var , "[i] , glm" , var_suffix[f] , "[i] , " , cutsname , " );\n" , sep="" )
                }
                
                if ( family[[f]]=="gamma" ) {
                    thetaname <- paste( "theta" , var_suffix[f] , sep="" )
                    m_gq <- paste( m_gq , indent , indent , "dev <- dev + (-2) * gamma_log( " , out_var , "[i] , exp(glm" , var_suffix[f] , "[i])*" , thetaname , " , " , thetaname , " );\n" , sep="" )
                }
                
                # close loop for this formula
                m_gq <- paste( m_gq , indent , "}\n" , sep="" )
                
            }#f
            
            # close generated quantities block
            m_gq <- paste( m_gq , "}\n " , sep="" )
            
            # add dev to pars_list
            pars_list <- c( pars_list , "dev" )
            
        }#calcDIC
    
    #########################################################################
    # build data list to pass to stan
    
        data_list <- list()
        
        for ( f in 1:num_formulas ) {
        
            # add index length
            index_name <- paste( "N" , var_suffix[f] , sep="" )
            data_list[[ index_name ]] <- nrow(fp[[f]]$dat)
            
            # add outcome
            outvar <- undot( fp[[f]]$yname )
            outvar <- paste( outvar , var_suffix[f] , sep="" )
            data_list[[ outvar ]] <- fp[[f]]$y
            
            # add variables
            for ( i in 1:ncol(fp[[f]]$dat) ) {
                name <- undot(colnames(fp[[f]]$dat)[i])
                if ( name!="Intercept" ) {
                    col_name <- paste( name , var_suffix[f] , sep="" )
                    data_list[[ col_name ]] <- fp[[f]]$dat[,i]
                }
            }
            
            # add max value of ordered outcome
            if ( family[[f]]=="ordered" ) {
                Kname <- paste( "K" , var_suffix[f] , sep="" )
                data_list[[ Kname ]] <- as.integer( max( fp[[f]]$y ) )
            }
            
        }#f
        
        # add N's for size of each cluster variable
        if ( length(cluster_vars) > 0 ) {
            for ( i in 1:length(cluster_size) ) {
                vname <- paste( "N_" , undot(names(cluster_size)[i]) , sep="" )
                data_list[[ vname ]] <- cluster_size[[i]]
            }
        # add diagonal matrices for Sigma priors
        # not used anymore, but were used for wishart cov priors
            for ( i in 1:length(cluster_vars) ) {
                nterms <- sum( cluster_vars[[i]] )
                if ( nterms > 1 ) {
                    name <- undot(names(cluster_vars)[i])
                    Oname <- paste( "Omega_" , name , sep="" )
                    data_list[[ Oname ]] <- diag(nterms)
                }
            }
        }
    
    ###################################################################
    # build inits to pass to stan
    
        init_list <- list() # 1 chain
        
        # fixed effects in each formula
        for ( f in 1:num_formulas ) {
            for ( i in 1:length(fp[[f]]$fixef) ) {
            
                varname <- undot(fp[[f]]$fixef[i])
                prefix <- fixed_prefix
                
                init_value <- 0 # corresponds to initmethod="zero", the default
                
                if ( varname=="Intercept" ) {
                    prefix <- ""
                    
                    # need code here to initalize intercepts, depending upon family
                    if ( family[[f]]=="binomial" ) {
                        init_value <- 0 # 0.5 probability; could also use logit
                    }
                    if ( family[[f]]=="gaussian" ) {
                        init_value <- mean( fp[[f]]$y ) # mean of outcome
                    }
                    if ( family[[f]]=="poisson" ) {
                        init_value <- log( mean( fp[[f]]$y ) ) # log mean of outcome
                    }
                    if ( family[[f]]=="gamma" ) {
                        init_value <- log( mean( fp[[f]]$y ) ) # gamma mu
                    }
                }
                
                parname <- paste( prefix , varname , var_suffix[f] , sep="" )
                
                if ( initmethod=="lme4" ) {
                    # extract estimate provided by lme4
                    init_value <- as.numeric( fixef(fit[[f]])[i] )
                }
                
                init_list[[ parname ]] <- init_value
            }#i
        }#f
        
        # varying effects spanning formulas
        # var-cov matrices default to diagonal
        if ( length(cluster_vars) > 0 ) {
            for ( i in 1:length(cluster_vars) ) {
                name <- undot(names( cluster_vars )[i])
                vname <- paste( vary_prefix , name , sep="" )
                num_effects <- sum( cluster_vars[[i]] )
                num_groups <- cluster_size[[i]]
                init_values <- matrix( 0 , nrow=num_groups , ncol=num_effects )
                if ( num_effects==1 ) init_values <- as.numeric( init_values )
                
                name_varcov <- "sigma_"
                init_varcov <- 1
                if ( num_effects > 1 ) {
                    if ( varpriors == "weak" ) {
                        # vector of initial standard deviations
                        init_varcov <- rep( 1 , num_effects )
                        # corr_matrix
                        name_rho <- paste( "Rho_" , name , sep="" )
                        init_list[[ name_rho ]] <- diag(num_effects)
                    } else {
                        init_varcov <- diag(num_effects)
                        name_varcov <- "Sigma_"
                    }
                }
                name_varcov <- paste( name_varcov , name , sep="" )
                
                if ( initmethod=="lme4" ) {
                    if ( num_effects==1 ) {
                        # pull out standard deviation estimate of varying intercepts
                        # need to be wary of zero edge estimate
                        
                    } else {
                        # pull out estimated var-cov matrix
                        # need to be wary of edge estimates, so blank for now
                        
                    }
                }
                
                init_list[[ vname ]] <- init_values
                init_list[[ name_varcov ]] <- init_varcov
            }
        }
        
        # specific family inits
        for ( f in 1:num_formulas ) {
            
            if ( family[[f]]=="gaussian" ) {
                signame <- paste( "sigma" , var_suffix[f] , sep="" )
                sig <- sd( fp[[f]]$y ) # a decent guess
                if ( initmethod=="lme4" ) sig <- sigma(fit[[f]]) # sigma() in lme4
                if ( sig < 0.001 ) sig <- 1
                init_list[[ signame ]] <- sig
            }
            
            if ( family[[f]]=="ordered" ) {
                # cutpoint inits...just use a vague spread for now
                cutsname <- paste( "cutpoints" , var_suffix[f] , sep="" )
                K <- as.integer( max( fp[[f]]$y ) )
                init_list[[ cutsname ]] <- seq( from=-2 , to=1.5 , length.out=K-1 )
            }
            
            if ( family[[f]]=="gamma" ) {
                thetaname <- paste( "theta" , var_suffix[f] , sep="" )
                # need theta init
                # need a good way to guess this value!
                ftheta <- function( pars ) {
                    parmu <- mean( fp[[f]]$y )
                    -sum( dgamma2( x=fp[[f]]$y , mu=parmu , scale=pars[1] , log=TRUE ) )
                }
                o <- optim( c( 1 ) , ftheta , method="SANN" ) # SANN not fast, but robust
                init_list[[ thetaname ]] <- 1/o$par[1]
            }
            
        }#f
        
    #########################################################
    # compose model code
    
    # close data block
    m_data <- paste( m_data , "}\n" , sep="" ) 
    
    # check for transformed data block
    if ( flag_transdata == FALSE ) {
        m_transdata <- ""
    } else {
        m_transdata <- paste( m_transdata1 , m_transdata2 , "}\n\n" , sep="" )
    }
    
    # check for transformed parameters block
    if ( flag_transpars == FALSE ) {
        m_transpars <- ""
    } else {
        m_transpars <- paste( m_transpars1 , m_transpars2 , "}\n\n" , sep="" )
    }
    
    model <- paste( m_data , "\n" , m_transdata , m_pars , "\n" , m_transpars , m_model , "\n" , m_gq , sep="" )
    
    #################################
    # sample
    
    if ( sample==FALSE ) {
    
        result <- list( data=data_list , model=model , init=list(init_list) , pars=pars_list[ 2:length(pars_list) ] , clusters=cluster_vars , clusters_size=cluster_size )
        
    } else {
    
        passinit <- initmethod
        if ( initmethod=="zero" | initmethod=="lme4" ) {
            passinit <- list()
            for ( i in 1:chains ) {
                passinit[[i]] <- init_list
            }
        }
        
        pars_list <- pars_list[ 2:length(pars_list) ]
        
        if ( verbose==TRUE ) {
            message( "Starting Stan model" )
            flush.console()
            start.time <- Sys.time()
        }
        if ( num_formulas == 1 ) {
            modelname <- paste( deparse( formula[[1]] ) , collapse="" )
        } else {
            # need to work on this
            modelname <- paste( deparse( formula[[1]] ) , collapse="" )
            # modelname <- paste( "[" , 1:num_formulas , "]" , sapply(formula,deparse) , collapse=" // " )
        }
        modelname <- paste( modelname , " [" , family[[1]] , "]" , sep="" )
        
        # ignite stan
        result <- stan( model_name = modelname , model_code = model , data = data_list , init = passinit , iter = iter , warmup = warmup , chains = chains , pars = pars_list , save_dso=FALSE , ... )
        
        if ( verbose==TRUE ) {
            message( show( Sys.time()-start.time ) )
            flush.console()
        }
    }
    
    ################################
    # calculate DIC, after fitting
    if ( calcDIC==TRUE & sample==TRUE ) {
        if ( verbose==TRUE ) {
            message( "Computing DIC" )
            flush.console()
        }
        gc()
        
        post <- extract( result , permuted=TRUE ) # extract as list of arrays
        # now process into expectations
        postbar <- list()
        for ( i in 1:length(post) ) {
            dims <- dim( post[[i]] )
            if ( length(dims)==1 )
                postbar[[ names(post)[i] ]] <- mean( post[[i]] )
            else
                postbar[[ names(post)[i] ]] <- apply( post[[i]] , 2:length(dims) , mean )
        }
        gc()
        
        Dbar <- postbar$dev
        Dhat <- 0
        
        # for each formula, accumulate deviance at expected values of parameters
        for ( f in 1:num_formulas ) {
        
            # compute varying component of GLM
            N_all <- nrow(fp[[f]]$dat)
            vary <- rep( 0 , N_all )
            if ( length( fp[[f]]$ranef ) > 0 ) {
                for ( i in 1:length( fp[[f]]$ranef ) ) {
                    cluster_name <- undot( names(fp[[f]]$ranef)[i] )
                    parname <- paste( vary_prefix , cluster_name , sep="" )
                    nterms <- length( fp[[f]]$ranef[[i]] )
                    total_terms <- sum( cluster_vars[[ cluster_name ]] ) # across all formulas
                    fnames <- fp[[f]]$ranef[[i]]
                    for ( j in 1:nterms ) {
                        jstart <- 0
                        if ( f > 1 ) {
                            jstart <- sum( cluster_vars[[ cluster_name ]][ 1:(f-1) ] )
                        }
                        dox <- 1
                        if ( undot(fnames[j]) != "Intercept" ) {
                            xname <- paste( undot(fnames[j]) , var_suffix[f] , sep="" )
                            dox <- data_list[[ xname ]]
                        }
                        cname <- paste( cluster_name , var_suffix[f] , sep="" )
                        if ( total_terms > 1 ) {
                            vary <- vary + postbar[[ parname ]][ data_list[[cname]] , j+jstart ] * dox
                        } else {
                            vary <- vary + postbar[[ parname ]][ data_list[[cname]] ] * dox
                        }
                    } #j
                } #i
            }
        
            # compute fixed component of GLM
            nterms <- length( fp[[f]]$fixef )
            fnames <- fp[[f]]$fixef
            glm <- rep( 0 , N_all )
            for ( j in 1:nterms ) {
                prefix <- fixed_prefix
                if ( undot(fnames[j])=="Intercept" ) {
                    if ( family[[f]]!="ordered" ) {
                        parname <- paste( undot(fnames[j]) , var_suffix[f] , sep="" )
                        glm <- glm + postbar[[ parname ]]
                    }
                } else {
                    parname <- paste( prefix , undot(fnames[j]) , var_suffix[f] , sep="" )
                    xname <- paste( undot(fnames[j]) , var_suffix[f] , sep="" )
                    dox <- data_list[[ xname ]]
                    glm <- glm + postbar[[ parname ]] * dox
                }
            }
            
            # now compute deviance (Dhat) at top level
            
            outvar <- paste( undot( fp[[f]]$yname ) , var_suffix[f] , sep="" )
            
            if ( family[[f]]=="gaussian" ) {
                sigmavar <- paste( "sigma" , var_suffix[f] , sep="" )
                Dhat <- Dhat + (-2) * sum( dnorm(
                    x = data_list[[ outvar ]] ,
                    mean = glm + vary ,
                    sd = postbar[[ sigmavar ]] ,
                    log = TRUE
                ) )
            }
            if ( family[[f]]=="binomial" ) {
                bintotvar <- paste( "bin_total" , var_suffix[f] , sep="" )
                Dhat <- Dhat + (-2) * sum( dbinom(
                    x = data_list[[ outvar ]] ,
                    prob = logistic( glm + vary ) ,
                    size = data_list[[ bintotvar ]] ,
                    log = TRUE
                ) )
            }
            if ( family[[f]]=="poisson" ) {
                Dhat <- Dhat + (-2) * sum( dpois(
                    x = data_list[[ outvar ]] ,
                    lambda = exp( glm + vary ) ,
                    log = TRUE
                ) )
            }
            if ( family[[f]]=="ordered" ) {
                Dhat <- Dhat + (-2) * sum( dordlogit(
                    x = data_list[[ outvar ]] ,
                    phi = glm + vary ,
                    a = postbar[[ "cutpoints" ]] ,
                    log = TRUE
                ) )
            }
            if ( family[[f]]=="gamma" ) {
                thetavar <- paste( "theta" , var_suffix[f] , sep="" )
                Dhat <- Dhat + (-2) * sum( dgamma2(
                    x = data_list[[ outvar ]] ,
                    mu = exp( glm + vary ) ,
                    scale = 1/postbar[[ thetavar ]] ,
                    log = TRUE
                ) )
            }
            
        } #f
        
        pD <- Dbar - Dhat
        DIC <- Dbar + pD
        #attr(result,"DIC") <- DIC
        #attr(result,"pD") <- pD
        #attr(result,"Dhat") <- Dhat
        #attr(result,"Dbar") <- Dbar
        cat( paste( "Expected deviance =" , round(Dbar,2) ) )
        cat( paste( "\nDeviance of expectation =" , round(Dhat,2) ) )
        cat( paste( "\npD =" , round(pD,2) ) )
        cat( paste( "\nDIC =" , round(DIC,2) ) )
        cat( "\n" )
    }
    
    ##########################################
    # prep results
    if ( extract==TRUE & sample==TRUE ) {
        gc()
        result <- extract( result , permuted=TRUE )
    }
    
    if ( calcDIC==TRUE & sample==TRUE ) {
        attr(result,"DIC") <- DIC
        attr(result,"pD") <- pD
        attr(result,"Dhat") <- Dhat
        attr(result,"Dbar") <- Dbar
    }
    
    # add model structure attributes
    calcgini <- function (x, weights = rep(1, length = length(x))) {
        ox <- order(x)
        x <- x[ox]
        weights <- weights[ox]/sum(weights)
        p <- cumsum(weights)
        nu <- cumsum(weights * x)
        n <- length(nu)
        nu <- nu/nu[n]
        sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
    }
    ranefattr <- list()
    if ( length(cluster_vars) > 0 ) {
        for ( i in 1:length( cluster_vars ) ) {
            name <- names( cluster_vars )[i]
            n <- cluster_size[[ name ]]
            # find grouping var in original data
            data2 <- as.data.frame(data)
            udCols <- sapply( colnames(data2) , undot )
            iOrig <- which( udCols==name )
            n_per_j <- sapply( 1:n , function(j) sum( data2[,iOrig]==j ) )
            gini <- calcgini( n_per_j )
            factors <- c()
            for ( f in 1:num_formulas ) {
                if ( !is.null( fp[[f]]$ranef[[name]] ) )
                    factors <- c( factors , fp[[f]]$ranef[[name]] )
            }
            ranefattr[[i]] <- list( factors=factors , n=n , gini=gini )
            names( ranefattr )[i] <- name
        }
    }
    attr( result , "ranef" ) <- ranefattr
    
    # add parsed formulas, so can use later to compute predictions, etc.
    attr( result , "formulas" ) <- list(
        fp = fp ,
        cluster_vars = cluster_vars ,
        cluster_size = cluster_size ,
        var_suffix = var_suffix ,
        family = family ,
        formula = formula
    )
        
    # return result
    result
}

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
    
    if ( !is.null( attr(fit,"DIC") ) ) {
        cat( paste( "\nDeviance: " , round(attr(fit,"Dhat"),1) , "   DIC: " , round(attr(fit,"DIC"),0) , "\n" , sep="" ) )
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

gitrstan <- function() {
    message( "Journey to mc-stan.org for a proper introduction to rstan." )
    message( "Attempting to install dependencies..." )
    install.packages('inline')
    install.packages('Rcpp')
    install.packages('RcppEigen')
    message( "Removing previous rstan installs (if any)..." )
    remove.packages('rstan')
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
