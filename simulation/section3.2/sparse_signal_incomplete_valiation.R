library(MASS)
library(glmnet)
####### generate signal variables
signal_data <- function(numx, numxs, size, XY, XSY, correlation){
  if(correlation == "1-2"){
    generate_XL_YL <- function(size, numx, numxs, XY, XSY){
      if ( XY =="Low" && XSY == "Low"){
        beta <- matrix( c( rep(0.5, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "Medium"){
        beta <- matrix( c( rep(0.5, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "High"){
        beta <- matrix( c( rep(0.5, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Low"){
        beta <- matrix( c( rep(1, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Medium"){
        beta <- matrix( c( rep(1, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "High"){
        beta <- matrix( c( rep(1, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Low"){
        beta <- matrix( c( rep(2, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Medium"){
        beta <- matrix( c( rep(2, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "High"){
        beta <- matrix( c( rep(2, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      
      
      mu <- rep(.1, numx) 
      sigma <- diag(1, numx)
      X <- mvrnorm(size, mu, sigma)
      Xstar <- matrix(nrow = size, ncol = numxs)
      for(S in 1:numxs){
        index <- sample(1:numx, 2, replace=FALSE)
        epsilon_Xstar <- mvrnorm(size, rep(0,1), diag(0.1, 1))
        Xstar[,S] <- sin(X[,index[1]]) +cos(X[,index[2]]) + epsilon_Xstar
      }
      
      XX.star <- cbind(X,Xstar)
      epsilon <- rnorm(size,0,.1) #generate noise
      
      Y <- XX.star %*% beta + epsilon #generate  return(list(X = X, Xstar = Xstar, Y = Y, XX.star = XX.star))
      df <- cbind(X,Xstar,Y)
      return(df)
    }
  }
  if(correlation == "1-4"){
    
    generate_XL_YL <- function(size, numx, numxs, XY, XSY){
      if ( XY =="Low" && XSY == "Low"){
        beta <- matrix( c( rep(0.5, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "Medium"){
        beta <- matrix( c( rep(0.5, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "High"){
        beta <- matrix( c( rep(0.5, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Low"){
        beta <- matrix( c( rep(1, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Medium"){
        beta <- matrix( c( rep(1, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "High"){
        beta <- matrix( c( rep(1, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Low"){
        beta <- matrix( c( rep(2, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Medium"){
        beta <- matrix( c( rep(2, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "High"){
        beta <- matrix( c( rep(2, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      
      mu <- rep(.1, numx) 
      sigma <- diag(1, numx)
      X <- mvrnorm(size, mu, sigma)
      Xstar <- matrix(nrow = size, ncol = numxs)
      for(S in 1:numxs){
        index <- sample(1:numx, 4, replace=FALSE)
        epsilon_Xstar <- mvrnorm(size, rep(0,1), diag(0.1, 1))
        Xstar[,S] <- sin(X[,index[1]]) + sin(X[,index[2]]) + cos(X[,index[3]]) + cos(X[,index[4]]) +epsilon_Xstar
      }
      
      XX.star <- cbind(X,Xstar)
      epsilon <- rnorm(size,0,.1) #generate noise
      
      Y <- XX.star %*% beta + epsilon #generate  return(list(X = X, Xstar = Xstar, Y = Y, XX.star = XX.star))
      df <- cbind(X,Xstar,Y)
      return(df)
    }
  }
  if(correlation == "1-1"){
    
    generate_XL_YL <- function(size, numx, numxs, XY, XSY){
      if ( XY =="Low" && XSY == "Low"){
        beta <- matrix( c( rep(0.5, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "Medium"){
        beta <- matrix( c( rep(0.5, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "High"){
        beta <- matrix( c( rep(0.5, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Low"){
        beta <- matrix( c( rep(1, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Medium"){
        beta <- matrix( c( rep(1, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "High"){
        beta <- matrix( c( rep(1, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Low"){
        beta <- matrix( c( rep(2, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Medium"){
        beta <- matrix( c( rep(2, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "High"){
        beta <- matrix( c( rep(2, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      
      mu <- rep(.1, numx) 
      sigma <- diag(1, numx)
      X <- mvrnorm(size, mu, sigma)
      Xstar <- matrix(nrow = size, ncol = numxs)
      for(S in 1:numxs){
        index <- sample(1:numx, 1, replace=FALSE)
        epsilon_Xstar <- mvrnorm(size, rep(0,1), diag(0.1, 1))
        Xstar[,S] <-  sin(X[,S]) + epsilon_Xstar
      }
      
      XX.star <- cbind(X,Xstar)
      epsilon <- rnorm(size,0,.1) #generate noise
      
      Y <- XX.star %*% beta + epsilon #generate  return(list(X = X, Xstar = Xstar, Y = Y, XX.star = XX.star))
      df <- cbind(X,Xstar,Y)
      return(df)
    }
  }
  if(correlation == "0-0"){
    generate_XL_YL <- function(size, numx, numxs, XY, XSY){
      if ( XY =="Low" && XSY == "Low"){
        beta <- matrix( c( rep(0.5, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "Medium"){
        beta <- matrix( c( rep(0.5, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Low" && XSY == "High"){
        beta <- matrix( c( rep(0.5, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Low"){
        beta <- matrix( c( rep(1, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "Medium"){
        beta <- matrix( c( rep(1, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="Medium" && XSY == "High"){
        beta <- matrix( c( rep(1, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Low"){
        beta <- matrix( c( rep(2, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "Medium"){
        beta <- matrix( c( rep(2, numx), rep(1, numxs)), ncol = 1 ) #coefficient
      } 
      if ( XY =="High" && XSY == "High"){
        beta <- matrix( c( rep(2, numx), rep(2, numxs)), ncol = 1 ) #coefficient
      } 
      
      
      mu <- rep(.1, (numx+numxs)) 
      sigma <- diag(1 , (numx + numxs)) #covariance matrix for X, X^\ast
      
      XX.star <- mvrnorm(size, mu, sigma) #Simulate X and X.star
      XX.star[,(numx+1) : (numx + numxs)] <- sin(XX.star[,(numx+1) : (numx + numxs)])
      X <- XX.star[, 1 : numx]
      Xstar <- XX.star[, (numx+1) : (numx + numxs)]
      epsilon <- rnorm(size,0,.1) #generate noise
      
      Y <- XX.star %*% beta + epsilon #generate  return(list(X = X, Xstar = Xstar, Y = Y, XX.star = XX.star))
      df <- cbind(X, Xstar, Y)
      return(df)
    }
  }
  
  data <- generate_XL_YL(size = size, numx = numx, 
                         numxs = numxs, XY = XY, XSY = XSY)
  colnames(data) <- c(paste0("X", 1:numx),paste0("Xstar",1:numxs),"Y")
  return(data)
}

noise_data <- function(number, study, study_size){
  noise <- matrix(nrow = study*study_size, ncol = number)

  mu <- rep(.1, number) 
  sigma <- matrix(0.2, ncol = number, nrow = number)
  sigma <- sigma + diag(0.8, number)
  
  noise <- mvrnorm(study*study_size, mu, sigma)
  
  return(data.frame(noise = noise))
}


#######################################
multi_study_comp <- function(gene_num){
  MSE  <- list()
  core.noise.name <- list()
  all.noise.name  <- list()
  core.noise.coef <- list()
  all.noise.coef  <- list()
  
  core.genes.name <- list()
  all.genes.name  <- list()
  core.genes.coef <- list()
  all.genes.coef  <- list()
  
  for(j in 1:300){
    set.seed(j)
    signal   <- signal_data(numx = 20, numxs = 20, size = 300, XY = "Low",
                            XSY = "High",correlation="1-2")
    
    noise    <- noise_data(number = 90 + 30, study  = 3, study_size = 100)
    
    full.data <- cbind(noise, signal)
    
    T1       <- full.data[1:100,]
    T2       <- full.data[101:200,]
    T3       <- full.data[201:300,]
    
    
    T1.available <- c(121:140, sample(141:160, 10), sample(c(1:120), 70) )
    T1.name      <- colnames(T1)[T1.available]
    T1.select    <- T1[, c(T1.available,161)]
    colnames(T1.select)[ncol(T1.select)] <- "ESR1"
    
    T2.available <- c(121:140, sample(141:160, 10), sample(c(1:120), 70) )
    T2.name      <- colnames(T2)[T2.available]
    T2.select    <- T2[, c(T2.available,161)]
    colnames(T2.select)[ncol(T2.select)] <- "ESR1"
    
    
    T3.available <- c(121:140, sample(141:160, 10), sample(c(1:120), 70) )
    T3.name      <- colnames(T3)[T3.available]
    T3.select    <- T3[, c(T3.available, 161)]
    colnames(T3.select)[ncol(T3.select)] <- "ESR1"
    
    
    T.name       <- cbind(T1.name, T2.name, T3.name)
    
    Study        <- cbind(rep("S1", gene_num), rep("S2", gene_num),
                          rep("S3", gene_num))
    gene.sig4 <- list()
    for(i in 1:3){
      gene.sig4[[i]] <- cbind(T.name[,i][1:gene_num], Study[,i])
    }
    
    expr.sig4 <- list()
    expr.sig4[[1]] <- T1.select[,c(1:gene_num, ncol(T1.select))]
    expr.sig4[[2]] <- T2.select[,c(1:gene_num, ncol(T2.select))]
    expr.sig4[[3]] <- T3.select[,c(1:gene_num, ncol(T3.select))]
    
    gene.full <- list()
    gene.full[[1]] <- T1.select
    gene.full[[2]] <- T2.select
    gene.full[[3]] <- T3.select
    
    
    combination <- combn(3,3)
    
    
    MSE2.5.2    <- c()
    
    count <- 1
    for(comb in 1:ncol(combination))  {
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
        ######## Polynomial Core Imputation #####
        #########################################
        
        T1.name <- colnames(T1)[-ncol(T1)]
        T2.name <- colnames(T2)[-ncol(T2)]
        T3.name <- colnames(T3)[-ncol(T3)]
        
        #### impute study-spcific genes of T2 in
        ####  T1 and Validation set
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
              
              X.for.pred3     <- as.matrix(cbind(T3[, intersect2], T3[, intersect2]^2))
              colnames(X.for.pred3) <- c(intersect2,paste0( intersect2,"sq"))
              T3.imput[,i]  <- predict(fit.poly, X.for.pred3)
            }
            
          }
        } 
        
        
        #### impute study-spcific genes of T1 in
        ####  T2 and Validation set
        name.1      <- T1.name[!T1.name %in% intersect2]
        index       <- which(!(T1.name %in% intersect2))
        if(length(index) > 0){
          T2.imput  <- matrix(NA, ncol=length(index), nrow=nrow(T2))
          T3.imput2 <- matrix(NA, ncol=length(index), nrow=nrow(T3))
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
              Y.train         <- as.matrix(c((T1[, name.miss]), T2[,name.miss]))
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
              
              X.for.pred3     <- as.matrix(cbind(T3[, intersect2], T3[, intersect2]^2))
              colnames(X.for.pred3) <- c(intersect2,paste0( intersect2,"sq"))
              T3.imput2[,i]  <- predict(fit.poly, X.for.pred3)
            }
            
          }
        } 
        
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
        
        
        #### Final predicting model
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
        
        non.zero          <- rownames(glm.in.lasso$beta)[which(glm.in.lasso$beta!=0)]
        first.cha         <- substr(non.zero,1, 1)
        noise.nonzero     <- non.zero[which(first.cha == "n")]
        genes.nonzero     <- non.zero[which(first.cha == "X")]
        
        core.noise.name[[j]] <- noise.nonzero
        core.noise.coef[[j]] <- glm.in.lasso$beta[noise.nonzero,]
        core.genes.name[[j]] <- genes.nonzero
        core.genes.coef[[j]] <- glm.in.lasso$beta[genes.nonzero,]
        
        
        
        ##########################################
        ########## Polynomial All Imputation #####
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
        
        #### impute study-spcific genes of T2 in
        ####  T1 and Validation set
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
              
              X.for.pred3     <- as.matrix(cbind(T3.full[, full.intersect2], T3.full[, full.intersect2]^2))
              colnames(X.for.pred3) <- c(full.intersect2,paste0( full.intersect2,"sq"))
              T3.imput[,i]  <- predict(fit.poly, X.for.pred3)
            }
            
          }
        } 
        
        
        #### impute study-spcific genes of T1 in
        ####  T2 and Validation set
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
              
              X.for.pred3     <- as.matrix(cbind(T3.full[, full.intersect2], T3.full[, full.intersect2]^2))
              colnames(X.for.pred3) <- c(full.intersect2,paste0( full.intersect2,"sq"))
              T3.imput2[,i]  <- predict(fit.poly, X.for.pred3)
            }
          }
        } 
        
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
        
        non.zero          <- rownames(glm.poly.out$beta)[which(glm.poly.out$beta!=0)]
        first.cha         <- substr(non.zero,1, 1)
        noise.nonzero     <- non.zero[which(first.cha == "n")]
        genes.nonzero     <- non.zero[which(first.cha == "X")]
        
        
        all.noise.name[[j]] <- noise.nonzero
        all.noise.coef[[j]] <- glm.poly.out$beta[noise.nonzero,]
        
        all.genes.name[[j]] <- genes.nonzero
        all.genes.coef[[j]] <- glm.poly.out$beta[genes.nonzero,]
        
        MSE2.5.2   <- rbind(MSE2.5.2, accuracy.temp)
        
        count <- count +1
        
        print(j)
      }
    }
    MSE[[j]]  <- MSE2.5.2
  }
  return(list(MSE = MSE,
              core.noise.name = core.noise.name,
              core.noise.coef = core.noise.coef,
              all.noise.name  = all.noise.name,
              all.noise.coef  = all.noise.coef,
              core.genes.name = core.genes.name,
              core.genes.coef = core.genes.coef,
              all.genes.name  = all.genes.name,
              all.genes.coef  = all.genes.coef))
}

top30 <- multi_study_comp(30)
top40 <- multi_study_comp(40)
top50 <- multi_study_comp(50)
top60 <- multi_study_comp(60)
top70 <- multi_study_comp(70)


save(top30, file="top30.RData")
save(top40, file="top40.RData")
save(top50, file="top50.RData")
save(top60, file="top60.RData")
save(top70, file="top70.RData")



