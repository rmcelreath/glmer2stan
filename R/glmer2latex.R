# convert glmer formula to a latex formatted maths model
# use nested index notation by default

glmer2latex <- function( formula , data , family="gaussian" ) {
    
    ########################################################################
    # check params
    
    if ( missing(formula) || missing(data) ) {
        stop( "Must provide both formula and data objects." )
    }
    
    for ( f in 1:length(family) ) {
        legal_families <- c( "gaussian" , "binomial" , "poisson" , "ordered" , "gamma" , "zigamma" , "zipoisson" )
        if ( !(family[[f]] %in% legal_families ) ) {
            stop( paste( "Family '" , family , "' not recognized." , sep="" ) )
        }
    }
    
    ########################################################################
    # private functions
    
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
    
    make_fixed_text <- function( fp , fixed_prefix="\\beta_" , suffix="" , drop_intercept=FALSE ) {
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
    
    family_to_distribution <- function( family ) {
        dist <- "UNKNOWN"
        if ( family=="gaussian" ) dist <- "Normal"
        if ( family=="binomial" ) dist <- "Binomial"
        if ( family=="poisson" ) dist <- "Poisson"
        if ( family=="ordered" ) dist <- "Ordered"
        if ( family=="gamma" ) dist <- "Gamma"
        paste( "\\text{" , dist , "}" , sep="" )
    }
    
    family_to_parameters <- function( family ) {
        result <- ""
        if ( family=="gaussian" ) {
            result <- " \\mu_i , \\sigma "
        }
        if ( family=="binomial" ) {
            result <- " p_i , n "
        }
        result
    }
    
    family_to_link <- function( family ) {
        result <- ""
        if ( family=="gaussian" ) {
            result <- "\\mu_i = "
        }
        if ( family=="binomial" ) {
            result <- "\\log \\frac{p_i}{1-p_i} = "
            # result <- "\\text{logit}(p_i) = "
        }
        result
    }
    
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

    ################
    # start building latex code
    
    out <- ""
    
    # level 1 distributions (top-level outcomes)
    for ( i in 1:num_formulas ) {
        out <- paste( out , "\\text{" , fp[[i]]$yname , "}_{i} ~ " , sep="" )
        out <- paste( out , family_to_distribution(family[[i]]) , "(" , sep="" )
        out <- paste( out , family_to_parameters( family[[i]] ) , ")\n" , sep="" )
    }
    
    # links to glm
    for ( i in 1:num_formulas ) {
        out <- paste( out , family_to_link( family[[i]] ) , sep="" )
        # build glm
        
    }
    
    ##############
    # result
    
    # print
    # "\t" is a tab to R's cat, so replace
    #out_cat <- gsub( "\t" , "\\t" , out , fixed=TRUE )
    cat( out ) 
    
    # return invisibly
    invisible(out) 

}

# x <- parse_formula(pulled.left~(1|actor)+prosoc.left,d)
# xx <- glmer2latex( pulled.left~(1|actor)+prosoc.left , d , "binomial" )