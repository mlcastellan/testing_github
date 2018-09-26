rm(list=ls())
getwd()
setwd("home/project_07_image_registration")
path=getwd()
set=DSCXXXX
ida=1
idb=2
#------------------------------------------
library(jpeg)
library(imager)
library(OpenImageR)
library(raster)
library(rgdal)
library(dplyr)
library(sp)
library(plyr)
#------------------------------------------
nim=paste0(path,set,ida,".JPG")
ima=imresize(as.cimg(readJPEG(nim)[,1:3099,1]),scale=scl)
str(im)
nim=paste0(path,set,idb,".JPG")
imB=imresize(as.cimg(readJPEG(nim)[,1:3099,1]),scale=scl)
str(im)
#-----------------------------------------
#---------FUNCTIONS-------------
#HARRIS - Harris corner detector
Harris<-function(im,sigma=2){
  eps=1.e-10
  ix=imgradient(im,"x")
  iy=imgradient(im,"y")
  ix2=isoblur(ix*ix,sigma,gaussian = T)
  iy2=isoblur(iy*iy,sigma,gaussian = T)
  ixy=isoblur(ix*iy,sigma,gaussian = T)
  (ix2*iy2-ixy*ixy)/(ix2+iy2+eps)
}
#-------------------------------
#Detect Keypoints
get.centers <- function(im,thr="99%",sigma=3*scl,bord=30*scl){
  dt <- Harris(im,sigma) %>% imager::threshold(thr) %>% label
  as.data.frame(dt) %>% subset(value>0 ) %>% dplyr::group_by(value) %>% dplyr::summarise(mx=round(mean(x)),my=round(mean(y))) %>% subset(mx>bord & mx<width(im)-bord & my>bord & my<height(im)-bord)
}
#------------------------------
#Return 1 or 2 main global image orientations restricted to (-45,45)
get_orientations<-function(im){
  ix=imgradient(im,"x")
  iy=imgradient(im,"y")
  ita=atan(iy/ix)*180/pi
  iga=table(sample(round(ita*2)/2,200000))
  #plot(iga)
  ma1=max(iga)[1]
  m1=which(iga==ma1)
  theta_1=(as.numeric(names(m1)))
  iga[max((m1-20),0):min((m1+20),length(iga))]=0
  #plot(iga)
  ma2=max(iga)[1]
  m2=which(iga==ma2)
  theta_2=(as.numeric(names(m2)))
  if(theta_1>45) theta_1=theta_1-90
  if(theta_1<(-45))theta_1=theta_1+90
  if(theta_2>45) theta_2=theta_2-90
  if(theta_2<(-45))theta_2=theta_2+90
  if(abs(theta_1-theta_2)>5){
    return(c(theta_1,theta_2))
  }
  else{
    return(theta_1)
  }
}
#------------------------------
#Get oriented descriptors
get_descriptor_oriented<-function(im,theta,v){
  pm=get.stencil(im,stencil_ext,x=v[,1],y=v[,2])
  w=sqrt(length(pm))
  pm=as.cimg(pm,x=w,y=w)
  imr=imrotate(pm,-theta)
  ww=round(width(imr)/2)
  get.stencil(imr,stencil,x=ww,y=ww)
}
#------------------------------
# Estimate a homography h from points in P to points in p
est_homograph<-function(P,p){
  n=nrow(P)
  hh=NULL
  for(i in 1:n){
    a=t(c(p[i,],1))
    b=t(c(0,0,0))
    c=P[i,]
    d=-c%*%a
    hh=rbind(hh,cbind(rbind(c(a,b),c(b,a)),d))
  }
  h=t(matrix(svd(hh,nv=ncol(hh))$v[,9],nrow=3,ncol=3))
}
#------------------------------
#Apply homographyh to points in p
apply_homograph<-function(h,p){
  p1=t(cbind(p,1))
  q1=t(h%*%p1)
  q1=q1/q1[,3]
  q1[,1:2]
}
#------------------------------
#Robust homography estimation from p1 to p2. Return h and the list of inliers
ransac<-function(p1,p2,thresh=100,N=1000){
  n=nrow(p1)
  set.seed(12345)
  sn=c(1:n)
  flag=matrix(0,nrow=N,ncol=n)
  for(i in 1:N){
    smpl=sample(sn,4)
    pp1=p1[smpl,]
    pp2=p2[smpl,]
    h=est_homograph(pp2,pp1)
    p=apply_homograph(h,p1)
    d=rowSums((p-p2)^2)
    flag[i,]=as.numeric(d<thresh)
  }
  sinliers=rowSums(flag)
  sinliers=sinliers[!is.na(sinliers)]
  imax=which(sinliers==max(sinliers))[1]
  inliers=sn[flag[imax,]==1]
  h=est_homograph(p2[inliers,],p1[inliers,])
  list(h,inliers)
}
#------------------------------
map.affine <- function(x,y) {
    p=apply_homograph(hm1,cbind(x,y))
    list(x=p[,1],y=p[,2])
  }
#------------------------------
#------PROGRAMA----------------
sigma_b=6*scl
ima_bl=isoblur(ima,sigma_b,gaussian = T)
imb_bl=isoblur(imb,sigma_b,gaussian = T)
tha=get_orientations(ima_bl)
thb=get_orientations(imb_bl)
#------------------------------
par(mfrow=c(1,2))
plot(imrotate(ima,-tha[1]))
plot(imrotate(imb,-thb[1]))
#------------------------------
par(mfrow=c(1,2))
plot(ima)
# KEYPOINT DETECTION
stencil <- expand.grid(dx=round(seq(-20,20,5)*scl),dy=round(seq(-20,20,5)*scl))
stencil_ext <- expand.grid(dx=round(seq(-30*scl,30*scl,1)),dy=round(seq(-30*scl,30*scl,1)))
#-----------------
kpa=as.data.frame(ima %>% get.centers(sigma=3*scl,"98%"))[,2:3]
kpa %$% points(mx,my,col="red")
#-----------------
plot(imb)
kpb=as.data.frame(imb %>% get.centers(sigma=3*scl,"98%"))[,2:3]
kpb %$% points(mx,my,col="red")
#------------------------------
#--GUARDO EN FEATA(features de a) y FEATB(features de b)los point descriptos de ambas imagenes
feata=NULL
for(theta in tha){
  dfa<-alply(kpa,1,function(v){ ss=get_descriptor_oriented(ima_bl,theta,v)}) %>% do.call(rbind,.)
  dfa=as.data.frame(t(apply(dfa,1,scale)))
  feata <- rbind(feata,dfa)
}
featb=NULL
for(theta in thb){
  dfb<- alply(kpb,1,function(v){ ss=get_descriptor_oriented(imb_bl,theta,v)})  %>% do.call(rbind,.)
  dfb=as.data.frame(t(apply(dfb,1,scale)))
  featb <- rbind(featb,dfb)
}
#------------------------------
# agrupo los puntos anteriores usando knn
require(FNN)

kk<-get.knnx(data=feata, query=featb, k=2, algorithm ="kd_tree" )
if(length(thb)==1){
  lpb=c(1:nrow(kpb))
}else{
  lpb=c(c(1:nrow(kpb),c(1:nrow(kpb))))
}
if(length(tha)==2)kpa=rbind(kpa,kpa)

mask=(kk$nn.dist[,1]/kk$nn.dist[,2]<.8)
match=cbind(kk$nn.index[mask,1],lpb[mask])               

p1=as.matrix(kpa[match[,1],])
p2=as.matrix(kpb[match[,2],])
#------------------------------
par(mfrow=c(1,1))
plot(kk$nn.dist[,1],kk$nn.dis[,2],pch='.')
points(kk$nn.dist[mask,1],kk$nn.dis[mask,2],pch='o',col="red")
#------------------------------
# filtro los puntos unsando el algoritmo RANSAC (random sample consensus)
hh=ransac(p1[,1:2],p2[,1:2],100,5000)
h=hh[[1]]
inliers=hh[[2]]
#---------------
print(paste0("Number of inliers: ",length(inliers)))
print("h=")
print(h)
#---------------
par(mfrow=c(1,2))
plot(ima)
kpa %$% points(mx,my,col="red")
points(p1[inliers,],col="green")
plot(imb)
kpb %$% points(mx,my,col="red")
points(p2[inliers,],col="green")
#------------------------------
#-We apply the transformation to the first image and compare with the second one:-----
hm1=solve(h)
imat=imwarp(ima,map=map.affine,dir="backward")
#------------------------------
par(mfrow=c(1,2))
plot(imat)
plot(imb)
#-despues de ajustar las dos imagenes las comparamos entre si
d1=imat-imb
d2=(imat-imb)^2*(imat>0)
par(mfrow=c(1,2))
plot(d1)
plot(log(d2+.0001))
#------------------------------