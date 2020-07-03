
###K-means를 기반으로 한 Ksphere
ksphere.conform <- function(train, test, k, a = 0.05){
  Data.kmeans <- kmeans(train, centers = k, iter.max = 10000)
  
  tr_n = nrow(train)  
  py = Data.kmeans$size/tr_n
  sig = Data.kmeans$withinss/Data.kmeans$size
  c <- Data.kmeans$centers
  n = nrow(test)
  d = ncol(test)
  
  
  ####metric on testset
  metric = matrix(0,nr = n,nc = k)
  
  for(i in 1:n){
    for(j in 1:k){
      metric[i,j] = sum((c[j,]-test[i,])^2)/sig[j] + 
        d*log(sig[j]) - 2*log(py[j]) 
    }
  }
  
  ####metric to cluster
  
  metric.res = matrix(0,nr = n,nc = 2)
  
  for(i in 1:n){
    metric.res[i,1] = min(metric[i,])
    metric.res[i,2] = which.min(metric[i,])
  }
  
  ####same procedure on Train set
  ##trainset
  train.metric = matrix(0,nr = tr_n,nc = k)

  for(i in 1:tr_n){
    for(j in 1:k){
      train.metric[i,j] = sum((c[j,]-train[i,])^2)/sig[j] + 
        d*log(sig[j]) - 2*log(py[j]) 
    }
  }
  
  ####metric to cluster
  
  train.res = matrix(0,nr = tr_n,nc = 2)
  
  for(i in 1:tr_n){
    train.res[i,1] = min(train.metric[i,])
    train.res[i,2] = which.min(train.metric[i,])
  }
  
  

  
  #유의수준
  t = matrix(0,nr = k,nc = length(a))
  cluster.alpha = matrix(rep(c(train.res[,2],metric.res[,2]),length(a)),
                         nr = tr_n + n,nc = length(a))
  
  for(j in 1:length(a)){
    for(i in 1:k){
      
      ###testset
      ind = which(metric.res[,2] == i)
      t[i,j] = quantile(metric.res[ind,1],1-a[j])
      cluster.alpha[tr_n + ind,j] <- sapply(ind,function(x){ifelse(metric.res[x,1]>t[i,j],
                                                                   0,metric.res[x,2])})
      ind1 = which(train.res[,2] == i)
      cluster.alpha[ind1,j] <- sapply(ind1,function(x){ifelse(train.res[x,1]>t[i,j],
                                                             0,train.res[x,2])})
    }
  }
  ###restore radius
  radius = matrix(0,nr = k,nc = length(a))
  
  for(j in 1:length(a)){
    for(i in 1:k){
      radius[i,j] = sqrt(max(t[i,j] + 2*log(py[i]) - 
                               d*log(sig[i]),0))*sqrt(sig[i])
    }
  }
  
  merging <- matrix(nr = k,nc = length(a))
  
  ###marging clusters to result
  if (k > 1){
    for(al in 1:length(a)){
      clu = 1:k###make indexx
      for(i in 1:(k-1)){
        for(j in (i+1):k){
          if(sqrt(sum((c[i,]-c[j,])^2))< radius[i,al] + radius[j,al]){
            mi = min(clu[i],clu[j])
            ma = max(clu[i],clu[j])
            clu[clu == ma] <- mi
          }  
        }
      }
      temp = unique(clu);u = 1
      for(tt in temp){
        clu[clu==tt] <- u
        u <- u+1
      }
      merging[,al] = clu
      for(i in 1:k){
        ind = which(cluster.alpha[,al] == i)
        cluster.alpha[ind,al] = clu[i]
      }
    }
  } else {merging[,al] <- 1}
  
  
  nm = paste0(100*(1-a),"%")
  colnames(t) = nm
  colnames(cluster.alpha) = nm
  colnames(merging) = nm
  colnames(radius) = nm
  
  row.nm = 1:k
  rownames(merging) = row.nm
  rownames(radius) = row.nm
  rownames(c) = row.nm
  
  res = list()
  res$centers = c;res$radius = radius
  res$cluster.train = cluster.alpha[1:tr_n,]
  res$cluster.test = cluster.alpha[(tr_n+1):(tr_n+n),]
  res$merge = merging
  return(res)
}

