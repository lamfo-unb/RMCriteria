
######################################  Gaussian Preference ########################################

kernelGaussian<-function(u){
  return((1/(sqrt(2*pi)))*exp(-0.5*(u^2)))
}

gaussianPreference<-function(delta,parms,COL){
  sigma2<-parms[COL,1]
  res<- ifelse(delta<0,0,(1-exp(-delta/(2*sigma2))))
  return(res)
}
integraFunctionGaussianPositive<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((gaussianPreference(val1-x,parms,COL)-gaussianPreference(val1-val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return positive flow
  return(sumKernel*gaussianPreference(x-val1,parms,COL))
}
integraFunctionGaussianNegative<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((gaussianPreference(-val1+x,parms,COL)-gaussianPreference(-val1+val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return negative flow
  return(sumKernel*gaussianPreference(-x+val1,parms,COL))
}



######################################  UsualPreference ############################################



usualPreference<-function(delta,parms){
  res<- ifelse(delta<0,0,1)
  return(res)
}
integraFunctionUsualPositive<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((usualPreference(val1-x,parms)-usualPreference(val1-val2,parms))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return positive flow
  return(sumKernel*usualPreference(x-val1,parms))
}
integraFunctionUsualNegative<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((usualPreference(-val1+x,parms)-usualPreference(-val1+val2,parms))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return negative flow
  return(sumKernel*usualPreference(-x+val1,parms))
}


######################################  UShapePreference ###########################################



ushapePreference<-function(delta,parms,COL){
  q<-parms[COL,1]
  res<- ifelse(delta<q,0,1)
  return(res)
}
integraFunctionUshapePositive<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((ushapePreference(val1-x,parms,COL)-ushapePreference(val1-val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return positive flow
  return(sumKernel*ushapePreference(x-val1,parms,COL))
}
integraFunctionUshapeNegative<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((ushapePreference(-val1+x,parms,COL)-ushapePreference(-val1+val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return negative flow
  return(sumKernel*ushapePreference(-x+val1,parms,COL))
}


######################################  VShapePreference ###########################################



vshapePreference<-function(delta,parms,COL){
  q<-parms[COL,1]
  if(delta<0){
    return(0)
  }
  else if(delta<=q){
    return(delta/q)
  }
  else{
    return(1)
  }
}
integraFunctionVshapePositive<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((vshapePreference(val1-x,parms,COL)-vshapePreference(val1-val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return positive flow
  return(sumKernel*vshapePreference(x-val1,parms,COL))
}
integraFunctionVshapeNegative<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((vshapePreference(-val1+x,parms,COL)-vshapePreference(-val1+val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return negative flow
  return(sumKernel*vshapePreference(-x+val1,parms,COL))
}


######################################  LevelPreference     ########################################



levelPreference<-function(delta,parms,COL){
  q<-parms[COL,1]
  p<-parms[COL,2]
  if(delta<q){
    return(0)
  }
  else if(delta<=p){
    return(0.5)
  }
  else{
    return(1)
  }
}
integraFunctionLevelPositive<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((levelPreference(val1-x,parms,COL)-levelPreference(val1-val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return positive flow
  return(sumKernel*levelPreference(x-val1,parms,COL))
}
integraFunctionLevelNegative<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((levelPreference(-val1+x,parms,COL)-levelPreference(-val1+val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return negative flow
  return(sumKernel*levelPreference(-x+val1,parms,COL))
}


###################################### Linear Preference ###########################################



linearPreference<-function(delta,parms,COL){
  q<-parms[COL,1]
  p<-parms[COL,2]
  if(delta<q){
    return(0)
  }
  else if(delta<=p){
    return((delta-q)/(p-q))
  }
  else{
    return(1)
  }
}

integraFunctionLinearPositive<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((linearPreference(val1-x,parms,COL)-linearPreference(val1-val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return positive flow
  return(sumKernel*linearPreference(x-val1,parms,COL))
}
integraFunctionLinearNegative<-function(x,ROW,COL,n,band,parms,datMat_temp){
  val1<-datMat_temp[ROW,COL]
  #Kernel Computation
  sumKernel<-0
  for(j in 1:n){
    val2<-datMat_temp[j,COL]
    sumKernel<-sumKernel+kernelGaussian((linearPreference(-val1+x,parms,COL)-linearPreference(-val1+val2,parms,COL))/band[COL,1])
  }
  sumKernel<-(1/(n*band[COL,1]))*sumKernel
  #Return negative flow
  return(sumKernel*linearPreference(-x+val1,parms,COL))
}


#' @title brutePrometheeIVKernel
#'
#' @description
#'   The PROMETHEE IV KERNEL method was developed by Albuquerque and Montenegro
#'   (2015),  as an alternative method to estimate PROMETHEE IV. It considers
#'   the empirical distribution of the criteria through kernel density
#'   estimation to evaluate alternatives.
#'
#'
#' @family RPromethee methods
#'
#' @aliases brutePrometheeIVKernel brutePrometheeIVKernel,RPrometheeArguments-method
#'
#' @param RPrometheeArguments An object with all RPromethee arguments. For
#' PROMETHEE IV KERNEL, the object must be supplied with a \code{band} argument,
#' for Kernel Density Estimation. See \code{\link{RPrometheeConstructor}} for
#' more information.
#'
#' @return
#'  \itemize{
#'   \item{PhiPlus} {The resulting PhiPlus from the alternatives for all
#'   criterias.}
#'   \item{PhiMinus} {The resulting PhiMinus from the alternatives for all
#'   criterias}
#'   \item{PhiNet} {The resulting PhiNet from the alternatives for all
#'   criterias}
#'  }
#'
#' @keywords decision-method mcda decision-analysis promethee
#'
#' @author Pedro Henrique Melo Albuquerque, \email{pedroa@@unb.br}
#' @author Gustavo Monteiro Pereira, \email{monteirogustavop@@gmail.com}
#' @export
brutePrometheeIVKernel<-function(datMat_temp, vecWeights, prefFunction, parms, band, normalize){
  #Step 1: Max ou Min orientation
  #inv<-(normalize==FALSE)
  #datMat_temp[,inv]<-(-1)*datMat_temp[,inv]
  n <- nrow(datMat_temp)
  datMat_temp <- datMat_temp
  matPlus<-matrix(NA,n,ncol(datMat_temp))
  matMinus<-matrix(NA,n,ncol(datMat_temp))
  #Step 2: Run the criterion
  for(COL in 1:ncol(datMat_temp)){
    #Gaussian Preference
    if(prefFunction[COL]==0){
      for(ROW in 1:n){
        matPlus[ROW,COL]<-integrate(integraFunctionGaussianPositive,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
        matMinus[ROW,COL]<-integrate(integraFunctionGaussianNegative,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
      }
    }
    else if(prefFunction[COL]==1){
      for(ROW in 1:n){
        matPlus[ROW,COL]<-integrate(integraFunctionUsualPositive,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
        matMinus[ROW,COL]<-integrate(integraFunctionUsualNegative,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
      }
    }
    else if(prefFunction[COL]==2){
      for(ROW in 1:n){
        q<-parms[COL,1]
        matPlus[ROW,COL]<-integrate(integraFunctionUshapePositive,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
        matMinus[ROW,COL]<-integrate(integraFunctionUshapeNegative,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
      }
    }
    else if(prefFunction[COL]==3){
      for(ROW in 1:n){
        p<-parms[COL,1]
        matPlus[ROW,COL]<-integrate(integraFunctionVshapePositive,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
        matMinus[ROW,COL]<-integrate(integraFunctionVshapeNegative,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
      }
    }
    else if(prefFunction[COL]==4){
      for(ROW in 1:n){
        q<-parms[COL,1]
        p<-parms[COL,2]
        matPlus[ROW,COL]<-integrate(integraFunctionLevelPositive,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
        matMinus[ROW,COL]<-integrate(integraFunctionLevelNegative,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
      }
    }
    else if(prefFunction[COL]==5){
      for(ROW in 1:n){
        q<-parms[COL,1]
        p<-parms[COL,2]
        matPlus[ROW,COL]<-integrate(integraFunctionLinearPositive,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
        matMinus[ROW,COL]<-integrate(integraFunctionLinearNegative,0,1,ROW=ROW,COL=COL,n=n,band=band,parms=parms,datMat_temp=datMat_temp)$value
      }
    }
  }

  #Step X: Create the dominance flow
  phiPlus<- rowSums(matPlus%*%diag(vecWeights))
  phiMinus<- rowSums(matMinus%*%diag(vecWeights))
  phiNet<-phiPlus-phiMinus

  ##Results
  res<-list()
  res[[1]]<-phiPlus
  res[[2]]<-phiMinus
  res[[3]]<-phiNet


  return(res)


}
