# function to translate a "map" set of formulas into likelihood and prior for Stan
# map() is a function in rethinking package, which is not required to use this interface.

map2stan <- function( flist , data , start , fit=TRUE , chains=1 , ... ) {
    
    ########################################
    # check arguments
    if ( missing(flist) ) stop( "Formula required." )
    if ( class(flist) != "list" ) {
        if ( class(flist)=="formula" ) {
            flist <- list(flist)
        } else {
            stop( "Formula or list of formulas required." )
        }
    }
    if ( missing(start) ) stop( "'start' list required." )
    if ( missing(data) ) stop( "'data' required." )
    if ( !( class(data) %in% c("list","data.frame") ) ) {
        stop( "'data' must be of class list or data.frame." )
    }
    
    # substitution list
    # need to translate R density functions into Stan density names
    # also need to map parameters into Stan parameters
    # in each list:
    # (1) name of Stan function
    # (2) list with sub-lists: (1) maps position of each R argument to positions of Stan arguments, (2) C variable type of argument (real or int, usually)
    # (3) Stan data type for outcome
    density_sub_list <- list(
        dnorm = list('normal',list(mean=list(1,"real"),sd=list(2,"real")),"real"), # mean, sd
        dbinom = list('binomial',list(size=list(1,"int"),prob=list(2,"real")),"int"), # size, prob
        dpois = list('poisson',list(lamba=list(1,"real")),"int"), # lambda
        dbetabinom = list('betabinomial',list(1,2),"int"), 
        dgampois = list('gammapoisson',list(1,2),"int"),
        dcauchy = list('cauchy',list(location=list(1,"real"),scale=list(2,"real")),"real"), # location, scale
        dunif = list('uniform',list(min=list(1,"real"),max=list(2,"real")),"real"), # min, max
        dlaplace = list('laplace',list(1,2),"real")
    )
    
    # link function substitutes
    link_sub_list <- list(
        'logistic' = 'inv_logit'
    )
    
    # size of indent in Stan code
    indent <- "    "
    
    # private function to convert formula to text
    formula2text <- function( f ) {
        arg_list <- as.character( f[[3]] )[-1]
        arg_text <- paste( arg_list , collapse=" , " )
        paste( f[[2]] , " ~ " , f[[3]][[1]] , "( " , arg_text , " )" , sep="" )
    }
    
    ########################################
    # convert formulas to text in proper call order
    # check for vectorized priors and expand
    # also build list of parameters with priors as we go
    pars_with_priors <- list()
    if ( length(flist) > 0 ) {
        flag_flatten <- FALSE
        for ( i in 1:length(flist) ) {
            if ( !(class(flist[[i]])=="formula") )
                stop( "Input not a formula." )
            
            # check for combined priors
            LHS <- flist[[i]][[2]]
            if ( class(LHS)=="call" ) {
                if ( as.character(LHS[[1]])=="c" ) {
                    # found a vector prior
                    newflist <- list()
                    num_pars <- length(LHS) - 1
                    for ( j in 1:num_pars ) {
                        newflist[[j]] <- flist[[i]]
                        newflist[[j]][[2]] <- LHS[[j+1]]
                        pars_with_priors[[ as.character(LHS[[j+1]]) ]] <- 1
                    } #j
                    flist[[i]] <- newflist
                    flag_flatten <- TRUE
                } else {
                    # a call other than c() detected on left hand side
                    stop( paste("Invalid prior specification:",flist[[i]]) )
                }
            } else {
                # not a call, so just record name of parameter
                pars_with_priors[[as.character(LHS)]] <- 1
            }
        } #i
        
        # flatten embedded lists that arise from expanded priors
        if ( flag_flatten ) flist <- unlist(flist,recursive=FALSE)
        
        # now strip out argument names in flattened list
        # and translate function names
        outcome_type <- ""
        for ( i in 1:length(flist) ) {
            # convert function name
            dname <- as.character( flist[[i]][[3]][[1]] )
            # check for unknown distributions
            if ( is.null( density_sub_list[[ dname ]] ) ) {
                stop( paste( "Density" , dname , "not recognized as having an equivalent Stan function." ) )
            }
            flist[[i]][[3]][[1]] <- as.symbol( density_sub_list[[ dname ]][[1]] )
            # save likelihood distribution, so we can type the outcome variable later
            # also swap out likelihood parameters for temp variable names to be used in Stan code
            if ( i==1 ) {
                outcome_type <- density_sub_list[[ dname ]][[3]]
                lik_pars <- density_sub_list[[ dname ]][[2]]
                par_type <- lik_pars
                for ( j in 1:length(lik_pars) ) {
                    temp_name <- paste( "lik_" , names(lik_pars)[j] , sep="" )
                    lik_pars[[j]] <- flist[[i]][[3]][[j+1]]
                    names(lik_pars)[j] <- temp_name
                    flist[[i]][[3]][[j+1]] <- as.symbol( temp_name )
                    par_type[j] <- density_sub_list[[ dname ]][[2]][[j]][[2]]
                }
            }
            # convert arguments
            num_arguments <- length(flist[[i]][[3]])
            names(flist[[i]][[3]]) <- rep("",num_arguments)
        }
    }
    
    ######################################
    # build Stan model code
    model_block <- "model {\n"
    
    # add variable definitions
    for ( i in 1:length(lik_pars) ) {
        model_block <- paste(
            model_block ,
            indent ,
            par_type[i] , " " ,
            names(lik_pars)[i] ,
            "[N];\n" ,
            sep=""
        )
    }
    
    # add priors
    for ( i in length(flist):2 ) {
        model_block <- paste(
            model_block ,
            indent ,
            formula2text( flist[[i]] ) , 
            ";\n" ,
            sep=""
        )
    }
    
    # add parameter formulas
    model_block <- paste( model_block , indent , "for ( i in 1:N ) {\n" , sep="" )
    for ( i in 1:length(lik_pars) ) {
        # replace all variable names in formulas with var_name[i]
        for ( j in 1:length(data) ) {
            var_name <- names(data)[j]
            sub_name <- paste( names(data)[j] , "[i]" , sep="" )
            lik_pars[i] <- gsub( var_name , sub_name , lik_pars[i] )
        }
        # translate link function names
        for ( j in 1:length(link_sub_list) ) {
            lik_pars[i] <- gsub( names(link_sub_list)[j] , link_sub_list[[j]] , lik_pars[i] )
        }
        # append
        model_block <- paste(
            model_block ,
            indent , indent ,
            names(lik_pars)[i] ,
            "[i] <- " ,
            lik_pars[i] ,
            ";\n" ,
            sep=""
        )
    }
    model_block <- paste( model_block , indent , "}\n" , sep="" )
    
    # add likelihood
    model_block <- paste(
        model_block ,
        indent ,
        formula2text( flist[[1]] ) , 
        ";\n" ,
        sep=""
    )
    
    model_block <- paste( model_block , "}\n" , sep="" )
    
    ######################################
    # build Stan data code
    data_block <- "data {\n"
    data <- as.list(data)
    
    num_cases <- length(data[[1]])
    data_block <- paste( data_block , indent , 
        "int N;\n" , sep=""
    )
    
    num_variables <- length(data)
    var_names <- names(data)
    outcome_var <- as.character( flist[[1]][[2]] )
    for ( i in 1:num_variables ) {
        var_type <- "real"
        if ( var_names[i] == outcome_var ) {
            var_type <- outcome_type
        }
        data_block <- paste( data_block , indent , 
            var_type , " " , var_names[i] , "[N];\n" , sep=""
        )
    }
    
    data_block <- paste( data_block , "}\n" , sep="" )
    
    ######################################
    # build Stan parameter code
    
    par_block <- "parameters {\n";
    
    num_pars <- length(start)
    par_names <- names(start)
    for ( i in 1:num_pars ) {
        par_block <- paste( par_block , indent , 
            "real " , par_names[i] , ";\n" , sep="" 
        )
    }
    
    par_block <- paste( par_block , "}\n" , sep="" )
    
    ######################################
    # result
    stan_code <- paste(
        data_block ,
        par_block ,
        model_block ,
        sep=""
    )
    
    if ( fit==TRUE ) {
        require(rstan)
        inits <- vector( mode="list" , length=chains )
        for ( i in 1:chains ) inits[[i]] <- start
        dat <- unlist( list( N=length(data[[1]]) , as.list(data) ) , recursive=FALSE)
        stanfit <- try( stan( model_code=stan_code , data=dat , chains=chains , init=inits , ... ) )
        return( stanfit )
    } else {
        return( list( stan_code , flist ) )
    }
    
}

# TEST CODE
if ( FALSE ) {

data(cars)

flist1 <- list(
    dist ~ dnorm( mean=a , sd=sigma ) ,
    a ~ dnorm(0,10) , 
    sigma ~ dcauchy(0,1)
)

flist2 <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
    c(a,b) ~ dnorm(0,10) , 
    sigma ~ dcauchy(0,1)
)

# fit0 <- map( flist2 , start=list(a=40,b=0.1,sigma=20) , data=cars )

stan_code <- map2stan( flist1 , start=list(a=40,sigma=20) , data=cars )

stan_code <- map2stan( flist2 , start=list(a=40,b=0.1,sigma=20) , data=cars , fit=TRUE )

# binomial test
library(rethinking)
cars2 <- cars
p <- logistic( 0.2*(cars2$dist - mean(cars2$dist)) )
n <- length(p)
cars2$bin <- ifelse( runif(n)<p , 1 , 0 )
# plot( bin ~ speed, cars2 )

flist0 <- list(
    bin ~ dbinom( size=1 , prob=logistic(a+b*speed) ) ,
    c(a,b) ~ dnorm(0,10)
)

fit <- map2stan( flist0 , start=list(a=0,b=0) , data=cars2 , fit=TRUE )

} # TEST CODE


