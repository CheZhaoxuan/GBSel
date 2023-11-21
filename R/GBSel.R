#' @title Perform GBSel to train and predict data sets
#'
#' @description
#' \code{GBSel}Perform GBSel to train and predict data sets
#'
#' @param y_train a phenotype vector of type numeric for training
#' @param X_train a genotype matrix of type numeric for training
#' @param X_test a genotype matrix of type numeric, the verification group(or no phenotype group) required to calculate GEBV
#' @param nIter a vector of number of iterations for three base models, Ridge, SVR, and KRR, in that order
#' @param fold the number of folds of cross-validation at the base layer
#' @param Parallel parallel or not
#' @param stacking whether to output the stacking results at the same time
#' 
#' @return A list of prediction values for train set and test set from GBSel and stacking(if stacking is TRUE)
#' @export
#' @importFrom foreach %dopar%
#' @importFrom stats "cor"
#' @importFrom stats "predict"
#' @importFrom stats "var"
#'
#' @examples
#' X <- genotype
#' X.FR <- Feature.reduce(X = X, method="GRM")
#' y.train <- phenotype[c(1:200),1]
#' X.train <- X.FR[c(1:200),]
#' X.test <- X.FR[-c(1:200),]
#' GBSel.fit <- GBSel(y_train=y.train,X_train=X.train,X_test=X.test,Parallel=FALSE)
#' 

GBSel <- function(
	y_train,
	X_train,
	X_test,
	nIter=c(10,10,10),
	fold=5,
	Parallel=FALSE,stacking=FALSE){
	# before run, need train.ML result of SVR, KRR, and Ridge
	#--------------------------------------------------------------------------------------------------------
	# Object: Perform GBSel to train and predict data sets
	#
	# Input:
	# y_train: a phenotype vector of type numeric for training
	# X_train: a genotype matrix of type numeric for training
	# X_test: a genotype matrix of type numeric, the verification group(or no phenotype group) required to calculate GEBV
	# nIter: a vector of number of iterations for three base models, Ridge, SVR, and KRR, in that order
	# fold: the number of folds of cross-validation at the base layer
	# Parallel: parallel or not
	# stacking: whether to output the stacking results at the same time
	#
	# Output:
	# A list of prediction values for train set and test set from GBSel and stacking(if stacking is TRUE)
	#--------------------------------------------------------------------------------------------------------
	#Create empty lists to store results   
	
	print("----------------------------------------< Welcome to GBSel >----------------------------------------")
	print("--------------------------------------- Environment loading... -------------------------------------")
	print("-------------------------------> Load sklearn and numpy from python.env --------------------------->")
	path <- reticulate::py_discover_config()
	py.lib <- python.env(env="python", env.path=path$python, module.check=TRUE)	
	sklearn <- py.lib$sklearn
	numpy <- py.lib$numpy
	print("---------------------------------------- Load successfully! ----------------------------------------")
	print("----------------------------------------------------------------------------------------------------")
	
	print("--------------------------------------- Inputdata checking... --------------------------------------")
	if(is.null(y_train)){stop("y_train must be required")}
	if(!is.numeric(y_train)){stop("y_train must be numeric")}
	if(!is.vector(y_train)){stop("y_train must be a vector")}
	if(is.null(X_train)){stop("X_train must be required")}
	if(!is.numeric(X_train)){stop("X_train must be numeric")}
	if(!is.matrix(X_train)){stop("X_train must be a matrix")}
	if(sum(is.na(X_train)) != 0){stop("NAs are not allowed in X_train")}
	if(length(y_train) != nrow(X_train)){stop("The individual dimensions of phenotype and genotype should be consistent")}
	if(is.null(X_test)){stop("X_test must be required")}
	if(!is.numeric(X_test)){stop("X_test must be numeric")}
	if(!is.matrix(X_test)){stop("X_test must be a matrix")}
	if(sum(is.na(X_test)) != 0){stop("NAs are not allowed in X_test")}

	print("----------------------------------------- Checking done! -------------------------------------------")
	print("----------------------------------------------------------------------------------------------------")
	
	print("--------------------------------------- GBSel starts running... ------------------------------------")
	print("----------------------------------------> Parameter tuning... ------------------------------------->")
	SVR.fit <- train.ML(y_train=y_train,X_train=X_train,X_test=X_test,method="SVR",py.lib=py.lib)
	KRR.fit <- train.ML(y_train=y_train,X_train=X_train,X_test=X_test,method="KRR",py.lib=py.lib)
	RR.fit <- train.ML(y_train=y_train,X_train=X_train,X_test=X_test,method="Ridge",py.lib=py.lib)
	Best_params <- list(best_params_SVR=SVR.fit$best_params, 
						best_params_KRR=KRR.fit$best_params, 
						best_params_Ridge=RR.fit$best_params)
						
	print("----------------------------------------> Training base model... ---------------------------------->")
	newtrain_Ridge <- vector("list", fold)
	newtrain_SVR <- vector("list", fold)
	newtrain_KRR <- vector("list", fold)
	newtest_Ridge <- vector("list", fold)
	newtest_SVR <- vector("list", fold)
	newtest_KRR <- vector("list", fold)
	nIter_m1=nIter[1]
	nIter_m2=nIter[2]
	nIter_m3=nIter[3]
	X1 <- cbind(y_train, X_train)
	fold1 <- cut(seq(1,nrow(X1)), breaks = fold, labels = FALSE)
	
	##stacking base learners
	if(Parallel){
		print("---------------------------------------> Parallel is choosen ----------------------------------->")	
		# Set up parallel processing
		num_cores <- parallel::detectCores()	
		print(paste0("------------------------------------> The number of CPU cores: ",num_cores," --------------------------------->"))	
		if(num_cores <= fold){cl <- parallel::makeCluster(num_cores)}else{cl <- parallel::makeCluster(fold)}
		print(paste0("----------------------------------------> Using cores: ",cl," ---------------------------------------->"))
		doParallel::registerDoParallel(cl)
		# Define the parallel foreach loop
		print("-------------------------------> The parallel foreach loop starts... ------------------------------>")
		Parallel_results <- foreach::foreach(j = 1:fold, .packages = c("foreach", "reticulate", "doParallel")) %dopar% {
			sklearn <- reticulate::import("sklearn") #加载sklearn
			numpy <- reticulate::import("numpy") #加载numpy
			testIndex <- which(fold1 == j, arr.ind = TRUE)
			basetrainID <- rownames(X_train[-testIndex,])
			basetestID <- rownames(X_train[testIndex,])
			y_train_base <- y_train[basetrainID]
			y_test_base <- y_train[basetestID]
			X_train_base <- X_train[basetrainID,]
			X_test_base <- X_train[basetestID,]
			T_val_Ridge <- matrix(NA, nrow = length(basetestID), ncol = (nIter_m1 + 1), dimnames = list(basetestID, paste0("Ridge_nIter", 0:nIter_m1)))
			T_val_SVR <- matrix(NA, nrow = length(basetestID), ncol = (nIter_m2 + 1), dimnames = list(basetestID, paste0("SVR_nIter", 0:nIter_m2)))
			T_val_KRR <- matrix(NA, nrow = length(basetestID), ncol = (nIter_m3 + 1), dimnames = list(basetestID, paste0("KRR_nIter", 0:nIter_m3)))
			T_test_Ridge <- matrix(NA, nrow = dim(X_test)[1], ncol = (nIter_m1 + 1), dimnames = list(rownames(X_test), paste0("Ridge_nIter", 0:nIter_m1)))
			T_test_SVR <- matrix(NA, nrow = dim(X_test)[1], ncol = (nIter_m2 + 1), dimnames = list(rownames(X_test), paste0("SVR_nIter", 0:nIter_m2)))
			T_test_KRR <- matrix(NA, nrow = dim(X_test)[1], ncol = (nIter_m3 + 1), dimnames = list(rownames(X_test), paste0("KRR_nIter", 0:nIter_m3)))
			y_train_boosted1 <- y_train_base
			y_train_boosted2 <- y_train_base
			y_train_boosted3 <- y_train_base
			#Ridge
			for(t1 in c(1:(nIter_m1+1))){
				Ridge_alpha <- Best_params$best_params_Ridge$alpha
				Ridge_solver <- Best_params$best_params_Ridge$solver
				model_Ridge <- sklearn$linear_model$Ridge(random_state=as.integer(0), alpha=Ridge_alpha, solver=Ridge_solver)
				model_Ridge$fit(X=X_train_base, y=y_train_boosted1)
				T_train_Ridge <- model_Ridge$predict(X_train_base)
				y_train_boosted1 <- y_train_boosted1 - T_train_Ridge
				T_val_Ridge[,t1] <- model_Ridge$predict(X_test_base)
				T_test_Ridge[,t1] <- model_Ridge$predict(X_test)
				if(t1 > 1){
					if(var(T_val_Ridge[,t1])==0){break}else{if(abs(cor(T_val_Ridge[,t1],T_val_Ridge[,t1-1]))>=0.99){break}}
				}
			}
			#SVR
			for(t2 in c(1:(nIter_m2+1))){
				SVR_kernel <- Best_params$best_params_SVR$kernel
				SVR_C <- Best_params$best_params_SVR$C
				SVR_gamma <- Best_params$best_params_SVR$gamma
				model_SVR <- sklearn$svm$SVR(kernel=SVR_kernel, C=SVR_C, gamma=SVR_gamma)
				model_SVR$fit(X=X_train_base, y=y_train_boosted2)
				T_train_SVR <- model_SVR$predict(X_train_base)
				y_train_boosted2 <- y_train_boosted2 - T_train_SVR
				T_val_SVR[,t2] <- model_SVR$predict(X_test_base)
				T_test_SVR[,t2] <- model_SVR$predict(X_test)
				if(t2 > 1){
					if(var(T_val_SVR[,t2])==0){break}else{if(abs(cor(T_val_SVR[,t2],T_val_SVR[,t2-1]))>=0.99){break}}
				}
			}
	
			#KRR
			for(t3 in c(1:(nIter_m3+1))){
				KRR_alpha <- Best_params$best_params_KRR$alpha
				KRR_kernel <- Best_params$best_params_KRR$kernel
				KRR_gamma <- Best_params$best_params_KRR$gamma
				model_KRR <- sklearn$kernel_ridge$KernelRidge(alpha=KRR_alpha, kernel=KRR_kernel, gamma=KRR_gamma, degree=as.integer(5))
				model_KRR$fit(X=X_train_base, y=y_train_boosted3)
				T_train_KRR <- model_KRR$predict(X_train_base)
				y_train_boosted3 <- y_train_boosted3 - T_train_KRR
				T_val_KRR[,t3] <- model_KRR$predict(X_test_base)
				T_test_KRR[,t3] <- model_KRR$predict(X_test)
				if(t3 > 1){
					if(var(T_val_KRR[,t3])==0){break}else{if(abs(cor(T_val_KRR[,t3],T_val_KRR[,t3-1]))>=0.99){break}}
				}
			}
			#Return the results
			list(T_val_Ridge = T_val_Ridge, T_val_SVR = T_val_SVR, T_val_KRR = T_val_KRR,
				T_test_Ridge = T_test_Ridge, T_test_SVR = T_test_SVR, T_test_KRR = T_test_KRR)
		}
		#Stop the parallel processing
		parallel::stopCluster(cl)
		print("-----------------------------------------> End parallel ------------------------------------------>")	
		for(j in 1:fold){
			newtrain_Ridge[[j]] <- Parallel_results[[j]]$T_val_Ridge
			newtrain_SVR[[j]] <-  Parallel_results[[j]]$T_val_SVR
			newtrain_KRR[[j]] <- Parallel_results[[j]]$T_val_KRR
			newtest_Ridge[[j]] <- Parallel_results[[j]]$T_test_Ridge
			newtest_SVR[[j]] <- Parallel_results[[j]]$T_test_SVR
			newtest_KRR[[j]] <- Parallel_results[[j]]$T_test_KRR
		}
	}else{
		for(j in 1:fold){
			testIndex <- which(fold1 == j, arr.ind = TRUE)
			basetrainID <- rownames(X_train[-testIndex,])
			basetestID <- rownames(X_train[testIndex,])
			y_train_base <- y_train[basetrainID]
			y_test_base <- y_train[basetestID]
			X_train_base <- X_train[basetrainID,]
			X_test_base <- X_train[basetestID,]
			T_val_Ridge <- matrix(NA,nrow=length(basetestID),ncol=(nIter_m1+1),dimnames=list(basetestID, paste0("Ridge_nIter",0:nIter_m1)))
			T_val_SVR <- matrix(NA,nrow=length(basetestID),ncol=(nIter_m2+1),dimnames=list(basetestID, paste0("SVR_nIter",0:nIter_m2)))
			T_val_KRR <- matrix(NA,nrow=length(basetestID),ncol=(nIter_m3+1),dimnames=list(basetestID, paste0("KRR_nIter",0:nIter_m3)))
			T_test_Ridge <- matrix(NA,nrow=dim(X_test)[1],ncol=(nIter_m1+1),dimnames=list(rownames(X_test), paste0("Ridge_nIter",0:nIter_m1)))
			T_test_SVR <- matrix(NA,nrow=dim(X_test)[1],ncol=(nIter_m2+1),dimnames=list(rownames(X_test), paste0("SVR_nIter",0:nIter_m2)))
			T_test_KRR <- matrix(NA,nrow=dim(X_test)[1],ncol=(nIter_m3+1),dimnames=list(rownames(X_test), paste0("KRR_nIter",0:nIter_m3)))
			y_train_boosted1 <- y_train_base
			y_train_boosted2 <- y_train_base
			y_train_boosted3 <- y_train_base
			#Ridge
			for(t1 in c(1:(nIter_m1+1))){
				Ridge_alpha <- Best_params$best_params_Ridge$alpha
				Ridge_solver <- Best_params$best_params_Ridge$solver
				model_Ridge <- sklearn$linear_model$Ridge(random_state=as.integer(0), alpha=Ridge_alpha, solver=Ridge_solver)
				model_Ridge$fit(X=X_train_base, y=y_train_boosted1)
				T_train_Ridge <- model_Ridge$predict(X_train_base)
				y_train_boosted1 <- y_train_boosted1 - T_train_Ridge
				T_val_Ridge[,t1] <- model_Ridge$predict(X_test_base)
				T_test_Ridge[,t1] <- model_Ridge$predict(X_test)
				if(t1 > 1){
					if(var(T_val_Ridge[,t1])==0){break}else{if(abs(cor(T_val_Ridge[,t1],T_val_Ridge[,t1-1]))>=0.99){break}}
				}
			}
			#SVR
			for(t2 in c(1:(nIter_m2+1))){
				SVR_kernel <- Best_params$best_params_SVR$kernel
				SVR_C <- Best_params$best_params_SVR$C
				SVR_gamma <- Best_params$best_params_SVR$gamma
				model_SVR <- sklearn$svm$SVR(kernel=SVR_kernel, C=SVR_C, gamma=SVR_gamma)
				model_SVR$fit(X=X_train_base, y=y_train_boosted2)
				T_train_SVR <- model_SVR$predict(X_train_base)
				y_train_boosted2 <- y_train_boosted2 - T_train_SVR
				T_val_SVR[,t2] <- model_SVR$predict(X_test_base)
				T_test_SVR[,t2] <- model_SVR$predict(X_test)
				if(t2 > 1){
					if(var(T_val_SVR[,t2])==0){break}else{if(abs(cor(T_val_SVR[,t2],T_val_SVR[,t2-1]))>=0.99){break}}
				}
			}
			#KRR
			for(t3 in c(1:(nIter_m3+1))){
				KRR_alpha <- Best_params$best_params_KRR$alpha
				KRR_kernel <- Best_params$best_params_KRR$kernel
				KRR_gamma <- Best_params$best_params_KRR$gamma
				model_KRR <- sklearn$kernel_ridge$KernelRidge(alpha=KRR_alpha, kernel=KRR_kernel, gamma=KRR_gamma, degree=as.integer(5))
				model_KRR$fit(X=X_train_base, y=y_train_boosted3)
				T_train_KRR <- model_KRR$predict(X_train_base)
				y_train_boosted3 <- y_train_boosted3 - T_train_KRR
				T_val_KRR[,t3] <- model_KRR$predict(X_test_base)
				T_test_KRR[,t3] <- model_KRR$predict(X_test)
				if(t3 > 1){
					if(var(T_val_KRR[,t3])==0){break}else{if(abs(cor(T_val_KRR[,t3],T_val_KRR[,t3-1]))>=0.99){break}}
				}
			}
			newtrain_Ridge[[j]] <- T_val_Ridge
			newtrain_SVR[[j]] <-  T_val_SVR
			newtrain_KRR[[j]] <- T_val_KRR
			newtest_Ridge[[j]] <- T_test_Ridge
			newtest_SVR[[j]] <- T_test_SVR
			newtest_KRR[[j]] <- T_test_KRR
			#if(t1==nIter_m1+1){print("Ridge not convergence")}else{print(paste0("Ridge convergence",": ",t1," iter"))}
			#if(t2==nIter_m2+1){print("SVR not convergence")}else{print(paste0("SVR convergence",": ",t2," iter"))}
			#if(t3==nIter_m3+1){print("KRR not convergence")}else{print(paste0("KRR convergence",": ",t3," iter"))}
		}
	}
	
	print("-------------------------------------> Gathering meta features... -------------------------------->")	
	##stacking meta learner
	#上面循环结束，产生新meta变量
	Metatrain_Ridge <- rbind(newtrain_Ridge[[1]],newtrain_Ridge[[2]],newtrain_Ridge[[3]],newtrain_Ridge[[4]],newtrain_Ridge[[5]])
	Metatrain_SVR <- rbind(newtrain_SVR[[1]],newtrain_SVR[[2]],newtrain_SVR[[3]],newtrain_SVR[[4]],newtrain_SVR[[5]])
	Metatrain_KRR <- rbind(newtrain_KRR[[1]],newtrain_KRR[[2]],newtrain_KRR[[3]],newtrain_KRR[[4]],newtrain_KRR[[5]])
	cbind_Ridge <- cbind(newtest_Ridge[[1]],newtest_Ridge[[2]],newtest_Ridge[[3]],newtest_Ridge[[4]],newtest_Ridge[[5]])
	cbind_SVR <- cbind(newtest_SVR[[1]],newtest_SVR[[2]],newtest_SVR[[3]],newtest_SVR[[4]],newtest_SVR[[5]])
	cbind_KRR <- cbind(newtest_KRR[[1]],newtest_KRR[[2]],newtest_KRR[[3]],newtest_KRR[[4]],newtest_KRR[[5]])
	Metatest_Ridge <- matrix(NA,nrow=dim(X_test)[1],ncol=nIter_m1+1)
	Metatest_SVR <- matrix(NA,nrow=dim(X_test)[1],ncol=nIter_m2+1)
	Metatest_KRR <- matrix(NA,nrow=dim(X_test)[1],ncol=nIter_m3+1)
	for(c in c(1:(nIter_m1+1))){
		cycleindex <- seq(c, by=nIter_m1+1,length.out=5)
		Metatest_Ridge[,c] <- Matrix::rowMeans(cbind_Ridge[,cycleindex],na.rm=TRUE)
	}
	for(c in c(1:(nIter_m2+1))){
		cycleindex <- seq(c, by=nIter_m2+1,length.out=5)
		Metatest_SVR[,c] <- Matrix::rowMeans(cbind_SVR[,cycleindex],na.rm=TRUE)
	}
	for(c in c(1:(nIter_m3+1))){
		cycleindex <- seq(c, by=nIter_m3+1,length.out=5)
		Metatest_KRR[,c] <- Matrix::rowMeans(cbind_KRR[,cycleindex],na.rm=TRUE)
	}
	Metatrain_Ridge <- Metatrain_Ridge[,apply(Metatrain_Ridge, 2, function(x){!any(is.na(x))})]
	Metatrain_SVR <- Metatrain_SVR[,apply(Metatrain_SVR, 2, function(x){!any(is.na(x))})]
	Metatrain_KRR <- Metatrain_KRR[,apply(Metatrain_KRR, 2, function(x){!any(is.na(x))})]
	Metatest_Ridge <- Metatest_Ridge[,1:dim(Metatrain_Ridge)[2]]
	Metatest_SVR <- Metatest_SVR[,1:dim(Metatrain_SVR)[2]]
	Metatest_KRR <- Metatest_KRR[,1:dim(Metatrain_KRR)[2]]
	MX_train <- cbind(Metatrain_Ridge, Metatrain_SVR, Metatrain_KRR)
	MX_test <- cbind(Metatest_Ridge, Metatest_SVR, Metatest_KRR)
	colnames(MX_test) <- colnames(MX_train)
	V=vector()
	for(col in 1:dim(MX_train)[2]){V[col] <- var(MX_train[,col])}
	index=(V!=0)
	MX_train <- MX_train[,index]
	MX_test <- MX_test[,index]
	
	print("---------------------------------------> Training meta model... ----------------------------------->")	
	#meta model
	model <- sklearn$linear_model$Lasso(random_state=as.integer(0))
	range_alpha <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0, 5, 10)
	lasso_parameters = list(alpha=range_alpha, max_iter=list(as.integer(10000000)))
	grid_search <- sklearn$model_selection$GridSearchCV(model, lasso_parameters, cv=as.integer(5))
	grid_search$fit(X=MX_train,y=y_train)
	MX_Te <- grid_search$predict(MX_test)
	MX_Tr <- grid_search$predict(MX_train)
	if(var(MX_Te)==0){
		lasso.model <- elasticnet::enet(as.matrix(MX_train),y_train,lambda=grid_search$best_params_$alpha)
		Tr <- predict(lasso.model,MX_train,type='fit')$fit
		Te <- predict(lasso.model,MX_test,type='fit')$fit
		MX_Tr <- Tr[,dim(Tr)[2]]
		MX_Te <- Te[,dim(Te)[2]]
	}
	
	if(stacking){
		MX_train_stac <- cbind(Metatrain_Ridge[,1], Metatrain_SVR[,1], Metatrain_KRR[,1])
		MX_test_stac <- cbind(Metatest_Ridge[,1], Metatest_SVR[,1], Metatest_KRR[,1])
		colnames(MX_test_stac) <- colnames(MX_train_stac)
		VV=vector()
		for(col in 1:dim(MX_train_stac)[2]){VV[col] <- var(MX_train_stac[,col])}
		index=(VV!=0)
		MX_train_stac <- MX_train_stac[,index]
		MX_test_stac <- MX_test_stac[,index]
		model <- sklearn$linear_model$Lasso(random_state=as.integer(0))
		range_alpha <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0, 5, 10)
		lasso_parameters = list(alpha=range_alpha, max_iter=list(as.integer(10000000)))
		grid_search <- sklearn$model_selection$GridSearchCV(model, lasso_parameters, cv=as.integer(5))
		grid_search$fit(X=MX_train_stac,y=y_train)
		MX_Te_stac <- grid_search$predict(MX_test_stac)
		MX_Tr_stac <- grid_search$predict(MX_train_stac)
		if(var(MX_Te_stac)==0){
			lasso.model.stac <- elasticnet::enet(as.matrix(MX_train_stac),y_train,lambda=grid_search$best_params_$alpha)
			Tr_stac <- predict(lasso.model.stac,MX_train_stac,type='fit')$fit
			Te_stac <- predict(lasso.model.stac,MX_test_stac,type='fit')$fit
			MX_Tr_stac <- Tr_stac[,dim(Tr_stac)[2]]
			MX_Te_stac <- Te_stac[,dim(Te_stac)[2]]
		}
		return(list(GBSel_PV_tr=MX_Tr,GBSel_PV_te=MX_Te,Stacking_PV_tr=MX_Tr_stac,Stacking_PV_te=MX_Te_stac))
	}else{
		return(list(GBSel_PV_tr=MX_Tr,GBSel_PV_te=MX_Te))
	}
	print("----------------------------------------- GBSel run over! ------------------------------------------")
	print("----------------------------------------------------------------------------------------------------")
	
}