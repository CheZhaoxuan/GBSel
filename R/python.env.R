#' @title Load the python environment, and import scikit-learn and numpy module
#'
#' @description
#' \code{python.env}Load the python environment, and import scikit-learn and numpy module
#'
#' @param env the environment of python, including python, conda, miniconda
#' @param env.path the path of python environment
#' @param module.check check that the required modules are installed
#'
#' @return environment and module
#' @export
#'
#' @examples
#' library(reticulate)
#' path <- reticulate::py_discover_config()
#' python.env(env="python", env.path=path$python, module.check=TRUE)
#' 


python.env <- function(env="python", env.path, module.check=TRUE){
	#--------------------------------------------------------------------------------------------------------
	# Object: Load the python environment, and import scikit-learn and numpy module
	#
	# Input:
	# env: the environment of python, including python, conda, miniconda
	# env.path: the path of python environment
	# module.check: check that the required modules are installed
	#
	# Output:
	# environment and module
	#--------------------------------------------------------------------------------------------------------
	if(!env %in% c("python","conda","miniconda")){stop("env shuold be one of 'python', 'conda', 'miniconda'")}
	if(is.null(reticulate::py_discover_config())){
		print("Check whether python is installed!")
	}else{
		if(env=="python"){reticulate::use_python(env.path)}
		if(module.check){
			if(reticulate::py_module_available("sklearn")){print("scikit-learn has been installed.")}else{stop("scikit-learn has not been installed!")}
			if(reticulate::py_module_available("numpy")){print("numpy has been installed.")}else{stop("numpy has not been installed!")}
		}
		sklearn <- reticulate::import("sklearn")
		numpy <- reticulate::import("numpy")
		#assign("sklearn", sklearn, envir = .GlobalEnv)
		#assign("numpy", numpy, envir = .GlobalEnv)
		return(list(sklearn=sklearn,numpy=numpy))		
	}
}
