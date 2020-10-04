###################################################################
#### Pre-process data using processExpressionSetList functions ####
###################################################################
library(curatedBreastData)
library(ggplot2)
library(ggthemes)
library(ROCR)
library(foreach)
library(doParallel)
data("clinicalData")
data("curatedBreastDataExprSetList")

proc_curatedBreastDataExprSetList <- processExpressionSetList(exprSetList=curatedBreastDataExprSetList[1:34], outputFileDirectory = "./", numTopVarGenes=1000)
esets <- lapply(1:length(proc_curatedBreastDataExprSetList), function(x) exprs(proc_curatedBreastDataExprSetList[[x]]))

preliminary <- function(gene_num){
  studies <- c(12093, 16646, 17705, 20181, 20194, 2034, 25055, 25065)
  index  <- c()
  gene.sig4 <- list()
  expr.sig4 <- list()
  gene.full <- list()
  pvalue     <- list()
  
  count=1
  for(i in 1:34) {
    study_num <- curatedBreastDataExprSetList[[i]]@phenoData@data[["study_ID.x"]][1]
    if(study_num %in% studies){
      XX.temp2  <- t(esets[[i]])
      index.ESR1 <- which(colnames(XX.temp2) == "ESR1")
      if(length(index.ESR1)>0){
        XX.temp3  <- data.frame(cbind(XX.temp2[,-index.ESR1], XX.temp2[,index.ESR1]))
        colnames(XX.temp3)[ncol(XX.temp3)] <- "ESR1"
        
        if( nrow(XX.temp3)>0 ){
          XX.temp4  <- XX.temp3[!is.na(XX.temp3$ESR1),]
          
          too.much  <- c()
          for(cols in 1:ncol(XX.temp4)){
            if(mean(is.na(XX.temp4[,cols]))>=0.05){
              too.much <- append(too.much, cols)
            }
          }
          if(length(too.much)>0){
            XX.temp4  <- XX.temp4[,-too.much]
          }
          
          if((nrow(XX.temp4)>=30)&(ncol(XX.temp4)>0)){
            
            ngenes    <- ncol(XX.temp4)-1
            nlpvalue  <- rep(NA, ngenes)
            for(lm in 1:ngenes){
              df.temp <- XX.temp4[,c(lm,ncol(XX.temp4))]
              if(quantile(df.temp[,1],0.25,na.rm=TRUE)<quantile(df.temp[,1],0.75,na.rm = TRUE)){
                nlpvalue[lm] <- summary(lm(ESR1~.,df.temp))$coefficients[2,4]
              } else{
                nlpvalue[lm] <- 999
              }
            }  

            XX.temp5  <- XX.temp4[, order(nlpvalue)[1:gene_num]]
            
            gene.temp <- cbind(colnames(XX.temp5),rep(study_num, length(colnames(XX.temp5))), 
                               rep(i, length(colnames(XX.temp5))))
            
            XX.temp6  <- cbind(XX.temp5, ESR1=XX.temp4$ESR1)
            
            
            gene.sig4[[count]]  <- data.frame(gene.temp)
            expr.sig4[[count]]  <- data.frame(XX.temp6)
            gene.full[[count]]  <- data.frame(XX.temp4)
            pvalue[[count]]     <- sort(nlpvalue)
            count <- count + 1
            print(i)
          }
        }
      }
    }
  }
  return(list(gene.sig4  = gene.sig4, expr.sig4  = expr.sig4,
              gene.full  = gene.full, count_stud = count - 1,
              pvalue     = pvalue))
}

multi_study_comp <- function(gene_num){
  data_pre <- preliminary(gene_num)
  
  gene.sig4  = data_pre$gene.sig4
  expr.sig4  = data_pre$expr.sig4
  gene.full  = data_pre$gene.full 
  combination0 <- combn(data_pre$count_stud,3)
  combination <- c()
  for(k in 1:ncol(combination0)){
    temp <- matrix(c(combination0[,k][1], combination0[,k][2], 
                     combination0[,k][3], combination0[,k][2],
                     combination0[,k][3], combination0[,k][1],
                     combination0[,k][1], combination0[,k][3],
                     combination0[,k][2]), ncol = 3, nrow = 3)
    combination <- cbind(combination,temp )
  }

  MSE2.5.2    <- c()
  
  core.nonzero <- c()
  all.nonzero.inter  <- c()
  all.nonzero.uninter  <- c()
  all.nonzero.full <- c()
  
  core.sig.coef.final <- c()
  all.sig.coef.final <- c()
  
  count <- 1
  for(comb in 1:ncol(combination))  {
    core.nonzero.temp <- c()
    all.nonzero.inter.temp  <- c()
    all.nonzero.uninter.temp  <- c()
    all.nonzero.full.temp <- c()
    
    intersect1 <- intersect(gene.sig4[[combination[,comb][1]]][,1],
                            gene.sig4[[combination[,comb][2]]][,1])
    intersect2 <- intersect(intersect1,gene.sig4[[combination[,comb][3]]][,1])
    accuracy.temp <- rep(NA, 3)
    deviance.Train.KT.temp <- c()
    #print(intersect2)
    if(length(intersect2)>0){
      study.index <- combination[,comb]
      
      T1  <- expr.sig4[[study.index[1]]]
      T2  <- expr.sig4[[study.index[2]]]
      T3  <- expr.sig4[[study.index[3]]]
      
      intersect1 <- intersect(colnames(T1)[-ncol(T1)], colnames(T2)[-ncol(T2)])
      intersect2 <- intersect(intersect1, colnames(T3)[-ncol(T3)])
      
      #### Omitting
      if(length(intersect2)==1){
        X.omit.train <- c(T1[,intersect2], T2[, intersect2])
        Y.omit.train <- c(T1[,"ESR1"], T2[,"ESR1"])
        any.na          <- which(is.na(Y.omit.train))
        if(length(any.na)>0){
          Y.omit.train       <- Y.omit.train[-any.na]
          X.omit.train       <- X.omit.train[-any.na,]
        }
        
        X.omit.train[which(is.na(X.omit.train))] <- mean(X.omit.train,na.rm=TRUE)
        
        df.omit      <- data.frame(X=X.omit.train, Y= Y.omit.train)
        lm.omit     <- lm(Y~., df.omit)
        accuracy.temp[1] <- mean((predict(lm.omit, data.frame(X=T3[,intersect2]))-T3[,"ESR1"])^2)
      } else{
        X.omit.train <- as.matrix(rbind(T1[,intersect2], T2[, intersect2]))
        Y.omit.train <- as.matrix(c(T1[,"ESR1"], T2[,"ESR1"]))
        any.na          <- which(is.na(Y.omit.train))
        if(length(any.na)>0){
          Y.omit.train       <- Y.omit.train[-any.na]
          X.omit.train       <- X.omit.train[-any.na,]
        }
        for(countna in 1:ncol(X.omit.train)){
          X.omit.train[which(is.na(X.omit.train[,countna])),countna] <- mean(X.omit.train[,countna],na.rm=TRUE)
        }
        cv_fit     <- cv.glmnet(X.omit.train, Y.omit.train,  alpha = 1, 
                                lambda = 10^seq(3, -2, by = -.1), 
                                standardize = FALSE)
        opt_lambda      <- cv_fit$lambda.min
        glm.omit.lasso  <- glmnet(X.omit.train, Y.omit.train, alpha = 1, 
                                  lambda = opt_lambda,standardize = FALSE)
        accuracy.temp[1] <- mean((predict(glm.omit.lasso, as.matrix(T3[,intersect2]))-T3[,"ESR1"])^2,
                                 na.rm=TRUE)
        
      }
      
      #########################################
      ####### Polynomial Core Imputation ######
      #########################################

      T1.name <- colnames(T1)[-ncol(T1)]
      T2.name <- colnames(T2)[-ncol(T2)]
      T3.name <- colnames(T3)[-ncol(T3)]
      
      #### impute T1 and T3
      name.2      <- T2.name[!T2.name %in% intersect2]
      index       <- which(!(T2.name %in% intersect2))
      if(length(index) > 0){
        T1.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T1))
        T3.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T3))
        colnames(T1.imput) <- name.2
        colnames(T3.imput) <- name.2
        for(i in (1:length(index))){
          index.1         <- index[i]
          name.miss       <- T2.name[index.1]
          if( (name.miss %in% T1.name) & (name.miss %in% T3.name)){
            T1.imput[,i]  <- T1[,name.miss]
            T3.imput[,i]  <- T3[,name.miss]
          }
          if( !(name.miss %in% T1.name) & !(name.miss %in% T3.name) ){
            X.train         <- as.matrix((cbind(T2[, intersect2], T2[, intersect2]^2)))
            colnames(X.train) <- c(intersect2,paste0( intersect2,"sq"))
            Y.train         <- as.matrix((T2[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)
            
            core.nonzero.temp <- append(core.nonzero.temp, mean(as.vector(fit.poly$beta)!=0))
            
            X.for.pred1     <- as.matrix(cbind(T1[, intersect2], T1[, intersect2]^2))
            colnames(X.for.pred1) <- c(intersect2,paste0( intersect2,"sq"))
            T1.imput[,i]  <- predict(fit.poly, X.for.pred1)
            X.for.pred3     <- as.matrix(cbind(T3[, intersect2], T3[, intersect2]^2))
            colnames(X.for.pred3) <- c(intersect2,paste0( intersect2,"sq"))
            T3.imput[,i]  <- predict(fit.poly, X.for.pred3)
          }
          if( !(name.miss %in% T1.name) & (name.miss %in% T3.name) ){
            X.train         <- as.matrix((cbind(T2[, intersect2], T2[, intersect2]^2)))
            colnames(X.train) <- c(intersect2,paste0( intersect2,"sq"))
            Y.train         <- as.matrix((T2[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)
            core.nonzero.temp <- append(core.nonzero.temp, mean(as.vector(fit.poly$beta)!=0))
            
            X.for.pred1     <- as.matrix(cbind(T1[, intersect2], T1[, intersect2]^2))
            colnames(X.for.pred1) <- c(intersect2,paste0( intersect2,"sq"))
            T1.imput[,i]  <- predict(fit.poly, X.for.pred1)
            
            T3.imput[,i]  <- T3[,name.miss]
          }
          if( (name.miss %in% T1.name) & !(name.miss %in% T3.name) ){
            T1.imput[,i]      <- T1[,name.miss]
            
            X.train.temp1     <- as.matrix((cbind(T2[, intersect2], T2[, intersect2]^2)))
            X.train.temp2     <- as.matrix((cbind(T1[, intersect2], T1[, intersect2]^2)))
            X.train           <- rbind(X.train.temp1, X.train.temp2) 
            
            colnames(X.train) <- c(intersect2,paste0( intersect2,"sq"))
            Y.train         <- as.matrix(c((T2[, name.miss]), T1[,name.miss]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)
            core.nonzero.temp <- append(core.nonzero.temp, mean(as.vector(fit.poly$beta)!=0))
            
            X.for.pred3     <- as.matrix(cbind(T3[, intersect2], T3[, intersect2]^2))
            colnames(X.for.pred3) <- c(intersect2,paste0( intersect2,"sq"))
            T3.imput[,i]  <- predict(fit.poly, X.for.pred3)
          }
          
        }
      } 
      
      
      #### impute T2 and T3
      name.1      <- T1.name[!T1.name %in% intersect2]
      index       <- which(!(T1.name %in% intersect2))
      if(length(index) > 0){
        T2.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T2))
        T3.imput2  <- matrix(NA, ncol=length(index), nrow=nrow(T3))
        colnames(T2.imput) <- name.1
        colnames(T3.imput2) <- name.1
        for(i in (1:length(index))){
          index.1         <- index[i]
          name.miss       <- T1.name[index.1]
          if( (name.miss %in% T2.name) & (name.miss %in% T3.name)){
            T2.imput[,i]  <- T2[,name.miss]
            T3.imput2[,i]  <- T3[,name.miss]
          }
          if( !(name.miss %in% T2.name) & !(name.miss %in% T3.name) ){
            X.train         <- as.matrix((cbind(T1[, intersect2], T1[, intersect2]^2)))
            colnames(X.train) <- c(intersect2,paste0( intersect2,"sq"))
            Y.train         <- as.matrix((T1[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)
            
            core.nonzero.temp <- append(core.nonzero.temp, mean(as.vector(fit.poly$beta)!=0))
            
            X.for.pred1     <- as.matrix(cbind(T2[, intersect2], T2[, intersect2]^2))
            colnames(X.for.pred1) <- c(intersect2,paste0( intersect2,"sq"))
            T2.imput[,i]  <- predict(fit.poly, X.for.pred1)
            X.for.pred3     <- as.matrix(cbind(T3[, intersect2], T3[, intersect2]^2))
            colnames(X.for.pred3) <- c(intersect2,paste0( intersect2,"sq"))
            T3.imput2[,i]  <- predict(fit.poly, X.for.pred3)
          }
          if( !(name.miss %in% T2.name) & (name.miss %in% T3.name) ){
            X.train         <- as.matrix((cbind(T1[, intersect2], T1[, intersect2]^2)))
            colnames(X.train) <- c(intersect2,paste0( intersect2,"sq"))
            Y.train         <- as.matrix((T1[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)
            core.nonzero.temp <- append(core.nonzero.temp, mean(as.vector(fit.poly$beta)!=0))
            
            X.for.pred1     <- as.matrix(cbind(T2[, intersect2], T2[, intersect2]^2))
            colnames(X.for.pred1) <- c(intersect2,paste0( intersect2,"sq"))
            T2.imput[,i]  <- predict(fit.poly, X.for.pred1)
            
            T3.imput2[,i]  <- T3[,name.miss]
          }
          if( (name.miss %in% T2.name) & !(name.miss %in% T3.name) ){
            T2.imput[,i]      <- T2[,name.miss]
            
            X.train.temp1     <- as.matrix((cbind(T2[, intersect2], T2[, intersect2]^2)))
            X.train.temp2     <- as.matrix((cbind(T1[, intersect2], T1[, intersect2]^2)))
            X.train           <- rbind(X.train.temp1, X.train.temp2) 
            
            colnames(X.train) <- c(intersect2,paste0( intersect2,"sq"))
            Y.train         <- as.matrix(c((T2[, name.miss]), T1[,name.miss]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)
            core.nonzero.temp <- append(core.nonzero.temp, mean(as.vector(fit.poly$beta)!=0))
            
            X.for.pred3     <- as.matrix(cbind(T3[, intersect2], T3[, intersect2]^2))
            colnames(X.for.pred3) <- c(intersect2,paste0( intersect2,"sq"))
            T3.imput2[,i]  <- predict(fit.poly, X.for.pred3)
          }
          
        }
      } 
      core.nonzero <- append(core.nonzero, mean(core.nonzero.temp))
      
      check.T1.name <- name.1
      check.T1.imput.name <- colnames(T1.imput)
      duplicate.imput.T1.index  <- which(check.T1.imput.name %in% check.T1.name)
      lendup1 <- length(duplicate.imput.T1.index)
      
      check.T2.name <- name.2
      check.T2.imput.name <- colnames(T2.imput)
      duplicate.imput.T2.index  <- which(check.T2.imput.name %in% check.T2.name)
      lendup2 <- length(duplicate.imput.T2.index)
      
      
      check.T3.imput.name  <- colnames(T3.imput)
      check.T3.imput2.name <- colnames(T3.imput2)
      duplicate.imput.T3.index  <- which(check.T3.imput2.name %in% check.T3.imput.name)
      lendup3 <- length(duplicate.imput.T3.index)
      
      if(length(intersect2)==1){
        if(lendup1!=0){
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput[,-duplicate.imput.T1.index])
        } else{
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput)
        }
        if(lendup2!=0){
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput[,-duplicate.imput.T2.index])
        } else{
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput)
        }
        colnames(T1.KT)[1] <- intersect2 
        colnames(T2.KT)[1] <- intersect2 
        T1.KT   <- T1.KT[,sort(colnames(T1.KT))]
        T2.KT   <- T2.KT[,sort(colnames(T2.KT))]
        T1.KT   <- cbind(T1.KT, ESR1=T1[,"ESR1"])
        T2.KT   <- cbind(T2.KT, ESR1=T2[,"ESR1"])
        
        
        if(lendup3!=0){
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2[,-duplicate.imput.T3.index])
        } else{
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2)
        }          
        colnames(T3.KT)[1] <- intersect2 
        
        T3.KT   <- T3.KT[,sort(colnames(T3.KT))]
        T3.KT   <- cbind(T3.KT, ESR1=T3[,"ESR1"])
        
      } else{
        if(lendup1!=0){
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput[,-duplicate.imput.T1.index])
        } else{
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput)
        }
        if(lendup2!=0){
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput[,-duplicate.imput.T2.index])
        } else{
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput)
        }
        T1.KT   <- T1.KT[,sort(colnames(T1.KT))]
        T2.KT   <- T2.KT[,sort(colnames(T2.KT))]
        T1.KT   <- cbind(T1.KT, ESR1=T1[,"ESR1"])
        T2.KT   <- cbind(T2.KT, ESR1=T2[,"ESR1"])
        
        
        if(lendup3!=0){
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2[,-duplicate.imput.T3.index])
        } else{
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2)
        }        
        T3.KT   <- T3.KT[,sort(colnames(T3.KT))]
        T3.KT   <- cbind(T3.KT, ESR1=T3[,"ESR1"])
      }
      
      
      
      
      T12.poly <- rbind(T1.KT, T2.KT)
      
      X.poly.train <- as.matrix(T12.poly[,-ncol(T12.poly)])
      Y.poly.train <- as.matrix(T12.poly[,"ESR1"])
      for(row in 1:nrow(X.poly.train)){
        X.temp <- X.poly.train[row,]
        index.na <- which(is.na(X.temp))
        if(length(index.na)>0){
          X.poly.train[row,index.na] <- mean(X.temp, na.rm=TRUE)
        }
      }
      cv_fit     <- cv.glmnet(X.poly.train, Y.poly.train, alpha = 1, 
                              lambda = 10^seq(3, -2, by = -.1), 
                              standardize = FALSE)
      opt_lambda      <- cv_fit$lambda.min
      glm.in.lasso  <- glmnet(X.poly.train, Y.poly.train, alpha = 1, 
                              lambda = opt_lambda,standardize = FALSE)
      accuracy.temp[3] <- mean((predict(glm.in.lasso, as.matrix(T3.KT[,-ncol(T3.KT)]))-T3[,"ESR1"])^2,
                               na.rm=TRUE)
      core.sig.coef.final <- append(core.sig.coef.final, 
                                     mean(as.vector(glm.in.lasso$beta)!= 0))
      
      
      
      
      ##########################################
      ####### Polynomial All Imputation ########
      ##########################################
      T1.full <- gene.full[[study.index[1]]]
      ESR1.index1 <- which(colnames(T1.full)=="ESR1")
      if(length(ESR1.index1)>0){
        T1.full <- T1.full[,-ESR1.index1]
      }
      T1.full.gene.name <- colnames(T1.full)
      
      T2.full <- gene.full[[study.index[2]]]
      ESR1.index2 <- which(colnames(T2.full)=="ESR1")
      if(length(ESR1.index2)>0){
        T2.full <- T2.full[,-ESR1.index2]
      }
      T2.full.gene.name <- colnames(T2.full)
      
      
      T3.full <- gene.full[[study.index[3]]]
      ESR1.index3 <- which(colnames(T3.full)=="ESR1")
      if(length(ESR1.index3)>0){
        T3.full <- T3.full[,-ESR1.index3]
      }
      T3.full.gene.name <- colnames(T3.full)
      
      
      full.intersect1 <- intersect(colnames(T1.full), colnames(T2.full))
      full.intersect2 <- intersect(full.intersect1, colnames(T3.full))
      
      T1.name <- colnames(T1)[-ncol(T1)]
      T2.name <- colnames(T2)[-ncol(T2)]
      T3.name <- colnames(T3)[-ncol(T3)]
      
      #### impute T1 and T3
      KT.outside.coef.intersected   <- c()
      KT.outside.coef.nointersected <- c()
      name.2      <- T2.name[!T2.name %in% intersect2]
      index       <- which(!(T2.name %in% intersect2))
      if(length(index) > 0){
        T1.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T1))
        T3.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T3))
        colnames(T1.imput) <- name.2
        colnames(T3.imput) <- name.2
        for(i in (1:length(index))){
          index.1         <- index[i]
          name.miss       <- T2.name[index.1]
          if( (name.miss %in% T1.full.gene.name) & (name.miss %in% T3.full.gene.name)){
            T1.imput[,i]  <- T1.full[,name.miss]
            T3.imput[,i]  <- T3.full[,name.miss]
          }
          if( !(name.miss %in% T1.full.gene.name) & !(name.miss %in% T3.full.gene.name) ){
            X.train         <- as.matrix((cbind(T2.full[, full.intersect2], T2.full[, full.intersect2]^2)))
            colnames(X.train) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            Y.train         <- as.matrix((T2[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)

            index.inter     <- which(rownames(fit.poly$beta) %in% c(intersect2, paste0( intersect2,"sq")))
            all.nonzero.inter.temp   <- append(all.nonzero.inter.temp,
                                          mean(fit.poly$beta[index.inter] != 0))
            all.nonzero.uninter.temp   <- append(all.nonzero.uninter.temp,
                                        mean(fit.poly$beta[-index.inter] != 0))
            all.nonzero.full.temp <- append(all.nonzero.full.temp,
                                       mean( as.vector(fit.poly$beta) != 0))
            X.for.pred1     <- as.matrix(cbind(T1.full[, full.intersect2], T1.full[, full.intersect2]^2))
            colnames(X.for.pred1) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T1.imput[,i]  <- predict(fit.poly, X.for.pred1)
            X.for.pred3     <- as.matrix(cbind(T3.full[, full.intersect2], T3.full[, full.intersect2]^2))
            colnames(X.for.pred3) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T3.imput[,i]  <- predict(fit.poly, X.for.pred3)
          }
          if( !(name.miss %in% T1.full.gene.name) & (name.miss %in% T3.full.gene.name) ){
            X.train         <- as.matrix((cbind(T2.full[, full.intersect2], T2.full[, full.intersect2]^2)))
            colnames(X.train) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            Y.train         <- as.matrix((T2[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)

            index.inter     <- which(rownames(fit.poly$beta) %in% c(intersect2, paste0( intersect2,"sq")))
            all.nonzero.inter.temp   <- append(all.nonzero.inter.temp,
                                               mean(fit.poly$beta[index.inter] != 0))
            all.nonzero.uninter.temp   <- append(all.nonzero.uninter.temp,
                                                 mean(fit.poly$beta[-index.inter] != 0))
            all.nonzero.full.temp <- append(all.nonzero.full.temp,
                                            mean( as.vector(fit.poly$beta) != 0))
            
            X.for.pred1     <- as.matrix(cbind(T1.full[, full.intersect2], T1.full[, full.intersect2]^2))
            colnames(X.for.pred1) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T1.imput[,i]  <- predict(fit.poly, X.for.pred1)
            
            T3.imput[,i]  <- T3.full[,name.miss]
          }
          if( (name.miss %in% T1.full.gene.name) & !(name.miss %in% T3.full.gene.name) ){
            T1.imput[,i]      <- T1.full[,name.miss]
            
            X.train.temp1     <- as.matrix((cbind(T2.full[, full.intersect2], T2.full[, full.intersect2]^2)))
            X.train.temp2     <- as.matrix((cbind(T1.full[, full.intersect2], T1.full[, full.intersect2]^2)))
            X.train           <- rbind(X.train.temp1, X.train.temp2) 
            
            colnames(X.train) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            Y.train         <- as.matrix(c((T2[, name.miss]), T1.full[,name.miss]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)

            index.inter     <- which(rownames(fit.poly$beta) %in% c(intersect2, paste0( intersect2,"sq")))
            all.nonzero.inter.temp   <- append(all.nonzero.inter.temp,
                                               mean(fit.poly$beta[index.inter] != 0))
            all.nonzero.uninter.temp   <- append(all.nonzero.uninter.temp,
                                                 mean(fit.poly$beta[-index.inter] != 0))
            all.nonzero.full.temp <- append(all.nonzero.full.temp,
                                            mean( as.vector(fit.poly$beta) != 0))
            
            X.for.pred3     <- as.matrix(cbind(T3.full[, full.intersect2], T3.full[, full.intersect2]^2))
            colnames(X.for.pred3) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T3.imput[,i]  <- predict(fit.poly, X.for.pred3)
          }
          
        }
      } 
      
      
      #### impute T2 and T3
      name.1      <- T1.name[!T1.name %in% intersect2]
      index       <- which(!(T1.name %in% intersect2))
      if(length(index) > 0){
        T2.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T2))
        T3.imput2  <- matrix(NA, ncol=length(index), nrow=nrow(T3))
        colnames(T2.imput) <- name.1
        colnames(T3.imput2) <- name.1
        for(i in (1:length(index))){
          index.1         <- index[i]
          name.miss       <- T1.name[index.1]
          if( (name.miss %in% T2.full.gene.name) & (name.miss %in% T3.full.gene.name)){
            T2.imput[,i]  <- T2.full[,name.miss]
            T3.imput2[,i] <- T3.full[,name.miss]
          }
          if( !(name.miss %in% T2.full.gene.name) & !(name.miss %in% T3.full.gene.name) ){
            X.train         <- as.matrix((cbind(T1.full[, full.intersect2], T1.full[, full.intersect2]^2)))
            colnames(X.train) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            Y.train         <- as.matrix((T1[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)

            index.inter     <- which(rownames(fit.poly$beta) %in% c(intersect2, paste0( intersect2,"sq")))
            all.nonzero.inter.temp   <- append(all.nonzero.inter.temp,
                                               mean(fit.poly$beta[index.inter] != 0))
            all.nonzero.uninter.temp   <- append(all.nonzero.uninter.temp,
                                                 mean(fit.poly$beta[-index.inter] != 0))
            all.nonzero.full.temp <- append(all.nonzero.full.temp,
                                            mean( as.vector(fit.poly$beta) != 0))
            
            X.for.pred2     <- as.matrix(cbind(T2.full[, full.intersect2], T2.full[, full.intersect2]^2))
            colnames(X.for.pred2) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T2.imput[,i]    <- predict(fit.poly, X.for.pred2)
            X.for.pred3     <- as.matrix(cbind(T3.full[, full.intersect2], T3.full[, full.intersect2]^2))
            colnames(X.for.pred3) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T3.imput2[,i]   <- predict(fit.poly, X.for.pred3)
          }
          if( !(name.miss %in% T2.full.gene.name) & (name.miss %in% T3.full.gene.name) ){
            X.train         <- as.matrix((cbind(T1.full[, full.intersect2], T1.full[, full.intersect2]^2)))
            colnames(X.train) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            Y.train         <- as.matrix((T1[, index.1]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)

            index.inter     <- which(rownames(fit.poly$beta) %in% c(intersect2, paste0( intersect2,"sq")))
            all.nonzero.inter.temp   <- append(all.nonzero.inter.temp,
                                               mean(fit.poly$beta[index.inter] != 0))
            all.nonzero.uninter.temp   <- append(all.nonzero.uninter.temp,
                                                 mean(fit.poly$beta[-index.inter] != 0))
            all.nonzero.full.temp <- append(all.nonzero.full.temp,
                                            mean( as.vector(fit.poly$beta) != 0))
            
            X.for.pred2     <- as.matrix(cbind(T2.full[, full.intersect2], T2.full[, full.intersect2]^2))
            colnames(X.for.pred2) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T2.imput[,i]  <- predict(fit.poly, X.for.pred2)
            
            T3.imput2[,i]  <- T3.full[,name.miss]
          }
          if( (name.miss %in% T2.full.gene.name) & !(name.miss %in% T3.full.gene.name) ){
            T2.imput[,i]      <- T2.full[,name.miss]
            
            X.train.temp1     <- as.matrix((cbind(T2.full[, full.intersect2], T2.full[, full.intersect2]^2)))
            X.train.temp2     <- as.matrix((cbind(T1.full[, full.intersect2], T1.full[, full.intersect2]^2)))
            X.train           <- rbind(X.train.temp1, X.train.temp2) 
            
            colnames(X.train) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            Y.train         <- as.matrix(c((T1[, name.miss]), T2.full[,name.miss]))
            any.na          <- which(is.na(Y.train))
            if(length(any.na)>0){
              Y.train       <- Y.train[-any.na]
              X.train       <- X.train[-any.na,]
            }
            for(countna in 1:ncol(X.train)){
              X.train[which(is.na(X.train[,countna])),countna] <- mean(X.train[,countna],na.rm=TRUE)
            }
            
            cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                         lambda = 10^seq(3, -2, by = -.1), 
                                         standardize = FALSE)
            opt_lambda      <- cv_fit$lambda.min
            fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                      lambda = opt_lambda,standardize = FALSE)

            index.inter     <- which(rownames(fit.poly$beta) %in% c(intersect2, paste0( intersect2,"sq")))
            all.nonzero.inter.temp   <- append(all.nonzero.inter.temp,
                                               mean(fit.poly$beta[index.inter] != 0))
            all.nonzero.uninter.temp   <- append(all.nonzero.uninter.temp,
                                                 mean(fit.poly$beta[-index.inter] != 0))
            all.nonzero.full.temp <- append(all.nonzero.full.temp,
                                            mean( as.vector(fit.poly$beta) != 0))
            
            X.for.pred3     <- as.matrix(cbind(T3.full[, full.intersect2], T3.full[, full.intersect2]^2))
            colnames(X.for.pred3) <- c(full.intersect2,paste0( full.intersect2,"sq"))
            T3.imput2[,i]  <- predict(fit.poly, X.for.pred3)
          }
        }
      } 
      
      all.nonzero.full <- append(all.nonzero.full, mean(all.nonzero.full.temp))
      all.nonzero.inter <- append(all.nonzero.inter, mean(all.nonzero.inter.temp))
      all.nonzero.uninter <- append(all.nonzero.uninter, mean(all.nonzero.uninter.temp))
      check.T1.name <- name.1
      check.T1.imput.name <- colnames(T1.imput)
      duplicate.imput.T1.index  <- which(check.T1.imput.name %in% check.T1.name)
      lendup1 <- length(duplicate.imput.T1.index)
      
      check.T2.name <- name.2
      check.T2.imput.name <- colnames(T2.imput)
      duplicate.imput.T2.index  <- which(check.T2.imput.name %in% check.T2.name)
      lendup2 <- length(duplicate.imput.T2.index)
      
      
      check.T3.imput.name  <- colnames(T3.imput)
      check.T3.imput2.name <- colnames(T3.imput2)
      duplicate.imput.T3.index  <- which(check.T3.imput2.name %in% check.T3.imput.name)
      lendup3 <- length(duplicate.imput.T3.index)
      
      if(length(intersect2)==1){
        if(lendup1!=0){
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput[,-duplicate.imput.T1.index])
        } else{
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput)
        }
        if(lendup2!=0){
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput[,-duplicate.imput.T2.index])
        } else{
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput)
        }
        colnames(T1.KT)[1] <- intersect2 
        colnames(T2.KT)[1] <- intersect2 
        T1.KT   <- T1.KT[,sort(colnames(T1.KT))]
        T2.KT   <- T2.KT[,sort(colnames(T2.KT))]
        T1.KT   <- cbind(T1.KT, ESR1=T1[,"ESR1"])
        T2.KT   <- cbind(T2.KT, ESR1=T2[,"ESR1"])
        
        
        if(lendup3!=0){
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2[,-duplicate.imput.T3.index])
        } else{
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2)
        }          
        colnames(T3.KT)[1] <- intersect2 
        
        T3.KT   <- T3.KT[,sort(colnames(T3.KT))]
        T3.KT   <- cbind(T3.KT, ESR1=T3[,"ESR1"])
        
      } else{
        if(lendup1!=0){
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput[,-duplicate.imput.T1.index])
        } else{
          T1.KT   <- cbind(T1[,intersect2], T1[,name.1], T1.imput)
        }
        if(lendup2!=0){
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput[,-duplicate.imput.T2.index])
        } else{
          T2.KT   <- cbind(T2[,intersect2], T2[,name.2], T2.imput)
        }
        T1.KT   <- T1.KT[,sort(colnames(T1.KT))]
        T2.KT   <- T2.KT[,sort(colnames(T2.KT))]
        T1.KT   <- cbind(T1.KT, ESR1=T1[,"ESR1"])
        T2.KT   <- cbind(T2.KT, ESR1=T2[,"ESR1"])
        
        
        if(lendup3!=0){
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2[,-duplicate.imput.T3.index])
        } else{
          T3.KT   <- cbind(T3[,intersect2], T3.imput, T3.imput2)
        }        
        T3.KT   <- T3.KT[,sort(colnames(T3.KT))]
        T3.KT   <- cbind(T3.KT, ESR1=T3[,"ESR1"])
      }
      
      
      
      T12.poly <- rbind(T1.KT, T2.KT)
      
      X.poly.train <- as.matrix(T12.poly[,-ncol(T12.poly)])
      Y.poly.train <- as.matrix(T12.poly[,"ESR1"])
      for(row in 1:nrow(X.poly.train)){
        X.temp <- X.poly.train[row,]
        index.na <- which(is.na(X.temp))
        if(length(index.na)>0){
          X.poly.train[row,index.na] <- mean(X.temp, na.rm=TRUE)
        }
      }
      cv_fit     <- cv.glmnet(X.poly.train, Y.poly.train, alpha = 1, 
                              lambda = 10^seq(3, -2, by = -.1), 
                              standardize = FALSE)
      opt_lambda      <- cv_fit$lambda.min
      glm.poly.out  <- glmnet(X.poly.train, Y.poly.train, alpha = 1, 
                              lambda = opt_lambda,standardize = FALSE)
      accuracy.temp[2] <- mean((predict(glm.poly.out, as.matrix(T3.KT[,-ncol(T3.KT)]))-T3[,"ESR1"])^2,
                               na.rm=TRUE)
      all.sig.coef.final <- append(all.sig.coef.final, 
                                    mean(as.vector(glm.poly.out$beta)!= 0))


      MSE2.5.2   <- rbind(MSE2.5.2, accuracy.temp)


      count <- count +1
      
      print(comb)
    }
  }
  
  return(list(MSE2.5.2 = MSE2.5.2,
              core.nonzero = core.nonzero,
              all.nonzero.inter = all.nonzero.inter,
              all.nonzero.uninter = all.nonzero.uninter,
              all.nonzero.full = all.nonzero.full,
              core.sig.coef.final = core.sig.coef.final,
              all.sig.coef.final = all.sig.coef.final))
}

re.30 <- multi_study_comp(30)
write.csv(re.30$MSE2.5.2, file="revised_30.csv")
save(re.30, file="revised_30.RData")






