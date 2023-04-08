#https://github.com/Manigben/GAPIT3
# loading packages for GAPIT and GAPIT functions
source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

# loading data set
myY=read.table("mdp_traits.txt", head = TRUE)
myGD=read.table("mdp_numeric.txt",head=T)
myGM=read.table("mdp_SNP_information.txt",head=T)
#myG=read.table("mdp_genotype_test.hmp.txt", head = FALSE)


# performing simulation phenotype
set.seed(198521)
Para=list(h2=0.7,NQTN=20)
mysimulation<-GAPIT(Para=Para,GD=myGD,GM=myGM)
myY=mysimulation$Y


myGAPIT <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model=c("GLM","MLM","SUPER","MLMM","FarmCPU","Blink"),# choose model
  #model=c("FarmCPU"),
  PCA.total=3,                                          # set total PCAs
  NJtree.group=4,                                       # set the number of clusting group in Njtree plot
  QTN.position=mysimulation$QTN.position,
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T
)




#GLM
#The GAPIT uses Least Squares to solve the model. The GAPIT code for running a GLM is:

  myGAPIT_GLM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="GLM",
  PCA.total=5,
  file.output=T
  )

#MLM
#EMMA method is used in GAPIT, the code of MLM is:

  myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="MLM",
  PCA.total=5,
  file.output=T
  )

#CMLM
#Compress Mixed Linear Model is published by Zhang in 2010. The code of CMLM is:

  myGAPIT_CMLM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="CMLM",
  PCA.total=5,
  file.output=T
  )

#MLMM
#Multiple Loci Mixied linear Model is published by Segura in 2012. The code of MLMM in GAPIT is:

  myGAPIT_MLMM <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="MLMM",
  PCA.total=5,
  file.output=T
  )

#SUPER
#Settlement of MLM Under Progressively Exclusive Relation- ship is published by Qishan in 2014. The code of SUPER is:

  myGAPIT_SUPER <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="SUPER",
  PCA.total=5,
  file.output=T
  )

#Farm-CPU
#Fixed and random model Circulating Probability Unification (FarmCPU) is published by Xiaolei in 2016. The code of Farm-CPU in GAPIT is:

  myGAPIT_FarmCPU <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="FarmCPU",
  PCA.total=5,
  file.output=T
  )


#GS
#gBLUP
#gBLUP used marker kinship to replace the pedgree relationship matrix. The code is:

  myGAPIT_gBLUP <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="gBLUP",
  PCA.total=5,
  file.output=T
  )

#cBLUP
#cBLUP used group kinship to replace the individual matrix. The code is:

  myGAPIT_cBLUP <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="cBLUP",
  PCA.total=5,
  file.output=T
  )

#sBLUP
#sBLUP used SUPER method to build psedue QTN kinship matrix. The code is:

  myGAPIT_sBLUP <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model="sBLUP",
  PCA.total=5,
  file.output=T
  )

