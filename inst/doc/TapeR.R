### R code from vignette source 'TapeR.Rnw'

###################################################
### code chunk number 1: TapeR.Rnw:48-49 (eval = FALSE)
###################################################
## install.packages("TapeR")


###################################################
### code chunk number 2: TapeR.Rnw:55-56
###################################################
library(TapeR)


###################################################
### code chunk number 3: TapeR.Rnw:60-61 (eval = FALSE)
###################################################
## ?TapeR


###################################################
### code chunk number 4: TapeR.Rnw:80-108
###################################################

#load example data
data(DxHx.df)

#prepare the data (could be defined in the function directly)
Id = DxHx.df[,"Id"]
x = DxHx.df[,"Hx"]/DxHx.df[,"Ht"]#calculate relative heights
y = DxHx.df[,"Dx"]

#plot the example data
plot(x,y,pch=as.character(Id),xlab="Relative height (m)", ylab="Diameter (cm)")


#define the relative knot positions and order of splines

knt_x = c(0.0, 0.1, 0.75, 1.0) # B-Spline knots: fix effects
ord_x = 4 # ord = order (4 = cubic)
knt_z = c(0.0, 0.1, 1.0); ord_z = 4 # B-Spline knots: rnd effects

#fit the model
taper.model <- TapeR_FIT_LME.f(Id, x, y, knt_x, ord_x, knt_z, ord_z,
IdKOVb = "pdSymm")

#save model parameters for documentation or dissimination
#parameters can be load()-ed and used to predict the taper
#or volume using one or several measured dbh
spruce.taper.pars <- taper.model$par.lme
#save(spruce.taper.pars, file="spruce.taper.pars.rdata")##uncomment to save


###################################################
### code chunk number 5: Rhelp
###################################################
#example data
data(DxHx.df)

#taper curve parameters based on all measured trees
data(SK.par.lme)

#select data of first tree
Idi <- (DxHx.df[,"Id"] == unique(DxHx.df$Id)[1])
tree1 <- DxHx.df[Idi,]

## Predict the taper curve based on the diameter measurement in 2 m
## height and known height 
tc.tree1 <- E_DHx_HmDm_HT.f(Hx=1:tree1$Ht[1], Hm=tree1$Hx[3],
Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 0, par.lme = SK.par.lme) 

#plot the predicted taper curve
plot(tc.tree1$Hx, tc.tree1$DHx, type="l")
#lower CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,1], lty=2)
#upper CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,3], lty=2)
#lower prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,1], lty=3)
#upper prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,3], lty=3)
#add measured diameter
points(tree1$Hx[3], tree1$Dx[3], pch=3, col=2)
#add the observations
points(tree1$Hx, tree1$Dx)

## Add the population average taper curve (without calibration) to the
## plot (not of high practical interest but good to know how to get it).
Ht  = tree1$Ht[1]
Hx  = tree1$Hx
#get fixed-effects design matrix for the Splines
X = TapeR:::XZ_BSPLINE.f(x=Hx/Ht, knt = SK.par.lme$knt_x, ord = SK.par.lme$ord_x)
#Calculate population average taper curve
DHx_PA = X %*% SK.par.lme$b_fix
#add to plot
lines(tree1$Hx, DHx_PA, lwd=2, lty=4)


###################################################
### code chunk number 6: Rhelp-p2
###################################################
## Let's see how CI's change if there's some uncertainty in the height
## measurement 
tc.tree1 <- E_DHx_HmDm_HT.f(Hx=1:tree1$Ht[1], Hm=tree1$Hx[3],
Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 1, par.lme = SK.par.lme) 

#plot the predicted taper curve
plot(tc.tree1$Hx, tc.tree1$DHx, type="l", xlab="Height (m)",
ylab="Diameter (cm)")
#lower CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,1], lty=2)
#upper CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,3], lty=2)
#lower prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,1], lty=3)
#upper prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,3], lty=3)
#add measured diameter
points(tree1$Hx[3], tree1$Dx[3], pch=3, col=2)
#add the observations
points(tree1$Hx, tree1$Dx)


###################################################
### code chunk number 7: RhelpExactCI
###################################################
plot(tc.tree1$Hx, tc.tree1$DHx, type="l", xlab="Height (m)",
ylab="Diameter (cm)")
## Calculate "exact" CIs. Careful: This takes a while!
#library(pracma)# for numerical integration with gaussLegendre()
tc.tree1.exact <- E_DHx_HmDm_HT_CIdHt.f(Hx=1:tree1$Ht[1],
Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 1, par.lme =
SK.par.lme) 
#add exact confidence intervals to approximate intervals above - fits
#quite well
lines(tc.tree1.exact[,1], tc.tree1.exact[,2], lty=2,col=2)
lines(tc.tree1.exact[,1], tc.tree1.exact[,4], lty=2,col=2)


## Calculate the height given a certain diameter threshold, say 8.5 cm
ht.tree1.d8.5 <- E_HDx_HmDm_HT.f (Dx=8.5, Hm=tree1$Hx[3],
Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 1, par.lme = SK.par.lme) 
#add to above created plot
abline(v=ht.tree1.d8.5)

## Calculate the timber volume
#for the whole stem
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme)

#Calculate the timber volume for a selected section given a height
#threshold (0.3 - 5 m)
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, A=0.3, B=5, iDH = "H")

#Calculate the timber volume for a selected section given a diameter
#threshold (30 cm - 15 cm) (negative value if A<B)
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, A=30, B=15, iDH = "D")

#The variance estimate resulting from the tree height uncertainty using
#a Normal approximation takes much longer...
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, IA=FALSE)
proc.time() - ptm


#... than the calculation using a 2-point distribution...
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, IA=TRUE)
proc.time() - ptm

#...fastest if no height variance is assumed
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 0, par.lme = SK.par.lme, IA=FALSE)
proc.time() - ptm

#Also the number of supportive points for the numerical integration
#influences the calculation time
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 0, par.lme = SK.par.lme, IA=FALSE, nGL=10)
proc.time() - ptm



