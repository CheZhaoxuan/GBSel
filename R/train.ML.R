#' @title Train single ML model and tune its parameters
#'
#' @description
#' \code{train.ML}Train a single machine learning model and tune its parameters
#'
#' @param y_train a phenotype vector of type numeric for training
#' @param X_train a genotype matrix of type numeric for training
#' @param X_test a genotype matrix of type numeric, the verification group(or no phenotype group) required to calculate GEBV
#' @param method a character, options are Support Vector Regression (SVR), Kernel Ridge Regression (KRR), Ridge Regression (Ridge)
#' @param cv the number of folds of cross-validation, the default is 5
#' @param py.lib load python module sklearn and numpy
#' @param parameters a list of user-specified parameters, the default is NULL
#'
#' @return A list of prediction values from single ML model and its optimal parameters
#' @export
#'
#' @examples
#' data(example_dataset)
#' X <- genotype
#' X.FR <- Feature.reduce(X = X, method="GRM")
#' y.train <- phenotype[c(1:200),1]
#' X.train <- X.FR[c(1:200),]
#' X.test <- X.FR[-c(1:200),]
#' path <- reticulate::py_discover_config()
#' python.env(env="python", env.path=path$python, module.check=TRUE)
#' py.lib <- python.env(env="python", env.path=path$python, module.check=TRUE)
#' SVR.fit <- train.ML(y_train=y.train,X_train=X.train,X_test=X.test,method="SVR",py.lib=py.lib)
#'

train.ML <- function(y_train,X_train,X_test,method,cv=5,py.lib=py.lib,parameters=NULL){
	#--------------------------------------------------------------------------------------------------------
	# Object: Train single machine learning models and tune parameters
	#
	# Input:
	# y_train: a phenotype vector of type numeric for training
	# X_train: a genotype matrix of type numeric for training
	# X_test: a genotype matrix of type numeric, the verification group(or no phenotype group) required to calculate GEBV
	# method: a character, options are Support Vector Regression (SVR), Kernel Ridge Regression (KRR), Ridge Regression (Ridge)
	# cv: the number of folds of cross-validation, the default is 5
	# py.lib: load python module sklearn and numpy
	# parameters: a list of user-specified parameters, the default is NULL
	#
	# Output:
	# A list of prediction values from single ML model and its optimal 		
	#--------------------------------------------------------------------------------------------------------
	if(!method %in% c("SVR","KRR","Ridge")){stop("method shuold be one of 'SVR', 'KRR', 'Ridge'!")}
	sklearn <- py.lib$sklearn
	numpy <- py.lib$numpy
	# 转换数据为Python格式
	y_train_py <- numpy$array(y_train)
	X_train_py <- numpy$array(X_train)
	X_test_py <- numpy$array(X_test)
	if(is.null(parameters)){
	#SVR
	if(method=="SVR"){
	    #print("You choose SVR model based on python.")
		model <- sklearn$svm$SVR()
		range_C <- c(1000.0, 100.0, 10.0, 1.0, 0.1)
		SVR_parameters <- list(kernel=list('rbf'), C=range_C, gamma=c('scale', 'auto'))
		grid_search <- sklearn$model_selection$GridSearchCV(model, SVR_parameters, cv=as.integer(cv))
		grid_search$fit(X=X_train_py,y=y_train_py)
	}
	#KRR
	if(method=="KRR"){
		#print("You choose KRR model based on python.")
		model <- sklearn$kernel_ridge$KernelRidge()
		range_alpha <- c(10.0, 1.0, 1e-1, 1e-2, 1e-3)
		range_gamma <- c(1e-1, 1e-2, 1e-3)
		KRR_parameters <- list(alpha=range_alpha, kernel=list('rbf'), gamma=range_gamma)
		grid_search <- sklearn$model_selection$GridSearchCV(model, KRR_parameters, cv=as.integer(cv))
		grid_search$fit(X=X_train_py,y=y_train_py)
	}
	#Ridge
	if(method=="Ridge"){
		#print("You choose Ridge model based on python.")
		model <- sklearn$linear_model$Ridge(random_state=as.integer(0))
		range_alpha = c(5, 10)
		range_solver = c('svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag', 'saga')
		Ridge_parameters <- list(alpha=range_alpha, solver=range_solver)
		grid_search <- sklearn$model_selection$GridSearchCV(model, Ridge_parameters, cv=as.integer(cv))
		grid_search$fit(X=X_train_py,y=y_train_py)
	}		
	best_params <-  grid_search$best_params_
	PV_tr <- grid_search$predict(X_train_py)
	PV_te <- grid_search$predict(X_test_py)
	return(list(best_params=best_params,PV_tr=PV_tr,PV_te=PV_te))
	}else{
	#SVR
	if(method=="SVR"){
	    #print("You choose SVR model based on python.")
		SVR_kernel <- parameters$best_params_SVR$kernel
		SVR_C <- parameters$best_params_SVR$C
		SVR_gamma <- parameters$best_params_SVR$gamma
		model <- sklearn$svm$SVR()
		model_SVR <- sklearn$svm$SVR(kernel=SVR_kernel, C=SVR_C, gamma=SVR_gamma)
		model_SVR$fit(X=X_train_py,y=y_train_py)
		PV_tr <- model_SVR$predict(X_train_py)
		PV_te <- model_SVR$predict(X_test_py)
	}
	#KRR
	if(method=="KRR"){
		#print("You choose KRR model based on python.")
		KRR_alpha <- parameters$best_params_KRR$alpha
		KRR_kernel <- parameters$best_params_KRR$kernel
		KRR_gamma <- parameters$best_params_KRR$gamma
		model_KRR <- sklearn$kernel_ridge$KernelRidge(alpha=KRR_alpha, kernel=KRR_kernel, gamma=KRR_gamma, degree=as.integer(5))
		model_KRR$fit(X=X_train_py,y=y_train_py)
		PV_tr <- model_KRR$predict(X_train_py)
		PV_te <- model_KRR$predict(X_test_py)
	}
	#Ridge
	if(method=="Ridge"){
		#print("You choose Ridge model based on python.")
		Ridge_alpha <- parameters$best_params_Ridge$alpha
		Ridge_solver <- parameters$best_params_Ridge$solver
		model_Ridge <- sklearn$linear_model$Ridge(random_state=as.integer(0), alpha=Ridge_alpha, solver=Ridge_solver)
		model_Ridge$fit(X=X_train_py,y=y_train_py)
		PV_tr <- model_Ridge$predict(X_train_py)
		PV_te <- model_Ridge$predict(X_test_py)
	}			
	best_params = parameters
	return(list(best_params=best_params,PV_tr=PV_tr,PV_te=PV_te))	
	}
	
}
