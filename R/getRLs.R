#' Convert gene expression profiles to a ranked list (data preprocessing)
#' 
#' Sorting the microarray probe-set identifiers according to the differential
#' expression values with respect to the untreated hybridization to obtaine a
#' ranked list. Gene-expression profiles in are represented in a nonparametric
#' fashion.
#' 
#' The genes on the array are rank-ordered according to their differential
#' expression relative to the control. First, control and treatment values less
#' than a primary threshold value (quartile) were set to that threshold value.
#' Finally, probe sets were ranked in descending order of d, where d is the
#' ratio of the corresponding treatment-to-control values. For probe sets where
#' d=1, a lower threshold was applied to the original difference values and a
#' new treatment to control ratio (d') calculated. These probe sets were then
#' sub-sorted in descending order of d.
#' 
#' @param control a matrix, including the vehicle control gene expression
#' profiles correspongding to the treatment gene expression profiles.
#' @param treatment a matrix, is composed of gene expression profiles.
#' @return A matrix is composed of ranked lists, a ranked list represents the
#' correspongding gene expression profiles.
#' @examples
#' 
#' 
#' if (require(GEOquery)){
#'   #treatment gene-expression profiles
#'   GSM118720 <- getGEO(filename=system.file("extdata/GSM118720.soft",package=
#'   "GeneExpressionSignature"))
#'   #control gene-expression profiles
#'   GSM118721 <- getGEO(filename=system.file("extdata/GSM118721.soft",package=
#'   "GeneExpressionSignature"))
#'   #data ranking according to the different expression values 
#'   control <- as.matrix(as.numeric(Table(GSM118721)[,2]))
#'   treatment <- as.matrix(as.numeric(Table(GSM118720)[,2]))
#'   ranked_list <-getRLs(control,treatment)
#' }
#' 
#' 
#' @export getRLs
getRLs <-
function(control,treatment){

   #check the feature of raw data (count, log, logratio)
   Check.Log <-function(control){
      logid=c()
      for (i in 1:ncol(control)){
          if (max(control[,i],na.rm=T)<50){
              if (abs(median(control[,i])-1)<0.01)
                  logid[i]='ratio'
              else if (abs(median(control[,i])) < 0.01)
                  logid[i] = 'logratio'
              else 
	          logid[i] = 'log';
           }
           else
               logid[i] = 'count'
       }
       return(logid)
   }
    
   #data preprocessing
   getrank <- function(control,treatment,logid){
      if (logid=='log' | logid=='logratio'){
           control=2^control
           treatment=2^treatment
       }
       control[control<0]=0
       treatment[treatment<0]=0 
       control_treatment=rbind(control,treatment)
       tmp=sort(control_treatment)
       th1=tmp[floor(length(control)/2)]
       th2=th1/10

       control_th=control
       treatment_th=treatment
       control_th[control<th1]=th1
       treatment_th[treatment_th<th1]=th1
       regulate=treatment_th/control_th
       tmp1=sort(regulate,decreasing=T,na.last=F)
       RL=order(regulate,decreasing=T,na.last=F)

       subRL=RL[which(tmp1==1)]
       subcontrol=control[subRL]
       subtreatment=treatment[subRL]
       subcontrol=pmax(subcontrol,th2)
       subtreatment=pmax(subtreatment,th2)
       subregulate=subtreatment/subcontrol
       sub_tmp1=sort(subregulate,decreasing=T)
       indx=order(subregulate,decreasing=T)
       subRL=subRL[indx]

       RL[which(tmp1==1)]=subRL
       RL=order(RL)
    }
    
    probe_num=nrow(control)
    treatment_num=ncol(control)
    RLs=matrix(0,probe_num,treatment_num)
    logid=Check.Log(control)
    for ( i in 1:treatment_num)
        RLs[,i]=getrank(control[,i],treatment[,i],logid[i])
    return(RLs)
}
