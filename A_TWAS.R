A_TWAS <- function(genotype, expression, GWAS_summary, CHR, LD, option = c("BLasso", "Horseshoe", "Horseshoe+"), extra_eQTL_weights = NA) {
    library(bayesreg)
    library(ACAT)

    ### extract GWAS summary statistics on certain CHR
    GWAS_summary <- GWAS_summary[which(GWAS_summary$CHR == CHR),]

    ### extract the data with SNPs containing in both genotype data and GWAS summary data
    idx1 <- which(is.na(match(names(genotype), GWAS_summary$BP)) !=T )
    genotype <- genotype[,..idx1]
        
    idx2 <- which(is.na(match(GWAS_summary$BP, names(genotype))) !=T )
    GWAS_summary <- GWAS_summary[idx2,]

    idx3 <- which(is.na(match(names(LD), names(genotype))) !=T )
    LD <- LD[idx3,idx3]

    ### data standardization
    genotype <- scale(genotype, center = TRUE, scale = TRUE)
    expression <- scale(expression, center = TRUE, scale = TRUE)
    data <- data.frame(genotype, expression)  

    z.Bls <- NA
    z.hs <- NA 
    z.hsplus <- NA

    ### running TWAS
    # Bayesian Lasso
    if(("BLasso" %in% option) == TRUE){
        model.Bls <- bayesreg(expression~.,data,prior="lasso")
        omit <- capture.output(model.Bls.s <- summary(model.Bls))
        beta.Bls <- as.vector(model.Bls.s$mu.coef)[-length(model.Bls.s$mu.coef)]
        z.Bls <- GWAS_summary$zscore%*%beta.Bls/sqrt(beta.Bls%*%LD%*%beta.Bls)
        print("Bayesian Lasso Finish")
    }

    # Horseshoe
    if(("Horseshoe" %in% option) == TRUE){
        model.hs <- bayesreg(expression~.,data,prior="hs")
        omit <- capture.output(model.hs.s <- summary(model.hs))
        beta.hs <- as.vector(model.hs.s$mu.coef)[-length(model.hs.s$mu.coef)]
        z.hs <- GWAS_summary$zscore%*%beta.hs/sqrt(beta.hs%*%LD%*%beta.hs)
        print("Horseshoe Finish")
    }

    # Horseshoe+
    if(("Horseshoe+" %in% option) == TRUE){
        model.hsplus <- bayesreg(expression~.,data,prior="hs+") 
        omit <- capture.output(model.hsplus.s <- summary(model.hsplus))
        beta.hsplus <- as.vector(model.hsplus.s$mu.coef)[-length(model.hsplus.s$mu.coef)]
        z.hsplus <- GWAS_summary$zscore%*%beta.hsplus/sqrt(beta.hsplus%*%LD%*%beta.hsplus)
        print("Horseshoe+ Finish")
    }

    ### use ACAT to combine the result

    z_temp <- c(z.Bls, z.hs, z.hsplus)
    p_temp <- c(pnorm(-abs(z_temp))*2)
    p_temp <- na.omit(p_temp)
    p_value <- ACAT(p)
    A_TWAS <- data.frame(p_value)
    result <- list(A_TWAS)
    names(result)[length(result)] <- "A-TWAS"

    if(("BLasso" %in% option) == TRUE){
        p_value <- pnorm(-abs(z.Bls))*2
        zscore <- z.Bls
        Bayesian_Lasso_Model <- list(data.frame(p_value, zscore))
        result <- append(result, Bayesian_Lasso_Model)
        names(result)[length(result)] <- "Bayesian Lasso Model"
    }

    if(("Horseshoe" %in% option) == TRUE){
        p_value <- pnorm(-abs(z.hs))*2
        zscore <- z.hs
        Horseshoe_Model <- list(data.frame(p_value, zscore))
        result <- append(result, Horseshoe_Model)
        names(result)[length(result)] <- "Horseshoe Model"
    }

    if(("Horseshoe+" %in% option) == TRUE){
        p_value <- pnorm(-abs(z.hsplus))*2
        zscore <- z.hsplus
        Horseshoeplue_Model <- list(data.frame(p_value, zscore))
        result <- append(result, Horseshoeplue_Model)
        names(result)[length(result)] <- "Horseshoe+ Model"
    }
    
    if(class(extra_eQTL_weights) != "logical"){
        idx4 <- which(is.na(match(extra_eQTL_weights$BP, GWAS_summary$BP)) !=T )
        extra_eQTL_weights <- extra_eQTL_weights[idx4,]
        for(i in 2:dim(extra_eQTL_weights)[2]){
            weight <- as.vector(extra_eQTL_weights[,i])
            zscore <- GWAS_summary$zscore%*%weight/sqrt(weight%*%LD%*%weight)
            p_value <- pnorm(-abs(zscore))*2
            external_model <- list(data.frame(p_value, zscore))
            result <- append(result, external_model)
            names(result)[length(result)] <- names(extra_eQTL_weights[i])
        }
        print("External Model Finish")
    }
    
    return(result)
}