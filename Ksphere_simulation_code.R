##functions Define

###Visualization을 위해 2차원에서 Sphere Point를 generate 
gen.sphere <- function(center, radius, n = 10000){
  intsec = seq(0,2*pi,length = n)
  Comp.1 = center[1] + radius*cos(intsec)
  Comp.2 = center[2] + radius*sin(intsec)
  return(data.frame(Comp.1,Comp.2))
}



###K-means를 기반으로 한 Ksphere
ksphere <- function(train, test, k, a = 0.05){
  Data.kmeans <- kmeans(train, centers = k, iter.max = 10000)
  
  c <- Data.kmeans$centers
  
  
  n = nrow(test)
  
  metric = matrix(0,nr = n,nc = k)
  
  
  for(i in 1:n){
    for(j in 1:k){
      metric[i,j] = sqrt(sum((c[j,]-test[i,])^2))
    }
  }
  
  
  metric.res = matrix(0,nr = n,nc = 2)
  
  for(i in 1:n){
    metric.res[i,1] = min(metric[i,])
    metric.res[i,2] = which.min(metric[i,])
  }
  
  #유의수준
  t = rep(0,k)
  
  
  for(i in 1:k){
    ind = which(metric.res[,2] == i)
    t[i] = quantile(metric.res[ind,1],1-a)
  }
  res = list()
  res$centers = c;res$radius = t
  
  return(res)
}


###n-dimension sphere explicit volume
Nvol = function(r,dimens){
  res = (sqrt(pi)*r)^dimens
  return(res/gamma(dimens/2 + 1))
}



###Volume Estimation by Importance Sampling
vol.est <- function(center, radius, sim.n = 1.0e+7, pi1 = rep(1,length(radius))){
  c = center;t = radius;k = dim(c)[1]
  di = dim(c)[2]
  C = Nvol(t,di)
  nan_ind = which(radius == 0)
  if(k == 1){
    sim.res = C
  } else {
    sim.res = 0
    
    
    C = pi1/C
    C[nan_ind] = 0
    
    for(i in 1:sim.n){
      a = as.vector(rmultinom(1, size = 1, prob = C))
      ran = rnorm(di)
      ran = ran/sqrt(sum(ran^2))
      ru = runif(1)*(t[as.logical(a)]^di)
      ru = ru^(1/di)
      sim = ran * ru
      aa = rep(1,k);aa = aa - a
      for(j in (1:k)[as.logical(aa)]){
        if(sqrt(sum((sim - c[j,])^2)) <= t[j])  a[j] = 1
      }
      sim.res = sim.res + 1/sum(a*C)
    }
    
    
    sim.res = sim.res/sim.n
  }
  
  return(sim.res)
}



###########simulation
set.seed(1)
library(ggplot2)
##generate Random Data


##setting 5 center
center = matrix(0,nr = 5,nc = 2)
center[2:5,]=cbind(c(4,-4,4,-4),c(4,4,-4,-4))
simd = matrix(nr = 5300,nc = 2)
label = vector(length = 5300)

##make save folder on wd
dir.create("./Ksphere")
dir.create("./Ksphere/result_per_k")

##generate normal number
simd[1:5000,] = cbind(rnorm(5000),rnorm(5000))
simd[5001:5300,] = cbind(runif(300,min = -10,max = 10),runif(300,min = -10,max = 10))

##labeling for simulation data(1:5000 -> normal from center, 5001:5003 : noise from uniform Dist)
for(i in 1:5){
  nm = paste0("Clust ", i)
  inn = seq(i*1000-999,i*1000)
  label[inn] = rep(nm,1000)
  simd[inn,] = simd[inn,] + matrix(rep(center[i,],each = 1000),nc = 2)
  
}
label[5001:5300] = "noise"

simdata = data.frame(simd,label)
colnames(simdata) = c("Comp.1","Comp.2","label")


###original data with label
viss <- ggplot(simdata,aes(x = Comp.1,y = Comp.2,color = label)) + geom_point()
ggsave("./Ksphere/Original_Data.png",plot = viss,width = 8,height = 7)


###split 2 data
trainind = rep(1000*0:4,each = 500) + rep(1:500,5)
trainind = c(trainind,sample((5001:5300),150))

train = simd[trainind,]
test = simd[-trainind,]


###record : volume of sphere on each cluster(k) and alpha(0.05, 0.1, 0.25)

record = matrix(0,nc = 10,nr = 3)
set.seed(1)
for(k in 1:10){
  A = ksphere.conform2(train,test,k,a = 0.05)
  B = ksphere.conform2(train,test,k,a = 0.1)
  C = ksphere.conform2(train,test,k,a = 0.25)
  pi1 = rep(1/k,k)
  ###set sim.n 1.0e+5 for speed so it might be unstable
  record[1,k] = vol.est(A$centers,A$radius,sim.n = 1.0e+5,pi1 = pi1)
  record[2,k] = vol.est(B$centers,B$radius,sim.n = 1.0e+5,pi1 = pi1)
  record[3,k] = vol.est(C$centers,C$radius,sim.n = 1.0e+5,pi1 = pi1)
  vis = ggplot(simdata,aes(x = Comp.1, y = Comp.2)) + geom_point()
  for(j in 1:k){
    vis = vis + geom_path(data = gen.sphere(A$centers[j,],A$radius[j]),color = "#db1aa8")
    vis = vis + geom_path(data = gen.sphere(B$centers[j,],B$radius[j]),color = "#b41adb")
    vis = vis + geom_path(data = gen.sphere(C$centers[j,],C$radius[j]),color = "#761ad9")
  }
  ###save each plot for cluster number on working directory
  filename = paste0("./Ksphere/result_per_k/clust",k,".png")
  ggsave(filename,plot = vis,width = 7,height = 7)
}


###visualize Volume vs K and Optimize K(K = 5)
Lab = rep(c("0.05","0.1","0.25"),10)
X = data.frame(rep(1:10,each = 3),as.vector(record),Lab)
colnames(X) = c("K","Volume","alpha")


###visualization of Volume vs K
vis2 = ggplot(X,aes(x = K,y = Volume,color = alpha)) + geom_line() + ggtitle("K vs Total Volume") + 
  theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))
vis2 <- vis2 + scale_x_continuous(breaks = 1:10)
vis2
ggsave("./Ksphere/volume_vs_K.png",plot = vis2)
record

###Kmeans on optimal K(for Check)
Kmeans.fit = kmeans(simd,5)
Kms.data = data.frame(simd,as.factor(Kmeans.fit$cluster))
colnames(Kms.data) = c("Comp.1","Comp.2","Cluster")


vis1 = ggplot(data = Kms.data,aes(x = Comp.1,y = Comp.2,color = Cluster)) + geom_point()
ggsave("./Ksphere/kmeans.png",plot = vis1,width = 7,height = 7)








#######################3
k = 5
A = ksphere.conform2(train,test,k,a = 0.1)
B = ksphere.conform2(train,test,k,a = 0.25)
C = ksphere.conform2(train,test,k,a = 0.5)
pi1 = rep(1/k,k)

vissa = viss
viss1 = viss
viss2 = viss
for(j in 1:k){
  vissa = vissa + geom_path(data = gen.sphere(A$centers[j,],A$radius[j]),color = "#b41adb")
  viss1 = viss1 + geom_path(data = gen.sphere(B$centers[j,],B$radius[j]),color = "#b41adb")
  viss2 = viss2 + geom_path(data = gen.sphere(C$centers[j,],C$radius[j]),color = "#761ad9")
}
filename = paste0("Ksphere/Final_clust",0.25,".png")
ggsave(filename,plot = vissa,width = 7,height = 7)

filename = paste0("Ksphere/Final_clust",0.25,".png")
ggsave(filename,plot = viss1,width = 7,height = 7)


filename = paste0("Ksphere/Final_clust",0.5,".png")
ggsave(filename,plot = viss2,width = 7,height = 7)

