##### function: read single sample coverage files #####
readcov <- function(file.prefix=x)
{
  bifile <- NULL
  trifile <- NULL
  tetrafile <- NULL
  pentafile <- NULL
  if(file.exists(paste0(file.prefix,"_bicov.txt"))){
    bifile <- read.table(paste0(file.prefix,"_bicov.txt"),
                         col.names = c("CovA","CovB","isStrict","VarType","VarId",
                                       "VarNum","VarDis"),stringsAsFactors = FALSE,
                         colClasses = c("CovA"="numeric","CovB"="numeric"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_bicov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tricov.txt"))){
    trifile <- read.table(paste0(file.prefix,"_tricov.txt"),
                          col.names = c("CovA","CovB","CovC","isStrict","VarType","VarId",
                                        "VarNum","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tricov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tetracov.txt"))){
    tetrafile <- read.table(paste0(file.prefix,"_tetracov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","isStrict","VarType","VarId",
                                          "VarNum","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tetracov.txt")," ) does not exists !"))
  }
  if(file.exists(paste0(file.prefix,"_pentacov.txt"))){
    pentafile <- read.table(paste0(file.prefix,"_pentacov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","CovE","isStrict","VarType","VarId",
                                          "VarNum","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_pentacov.txt")," ) does not exists !"))
  }  
  
  
  cov <- NULL
  
  if(dim(bifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=bifile$CovA+bifile$CovB,VarNum=bifile$VarNum,VarSize=bifile$VarType))
  }
  if(dim(trifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=trifile$CovA+trifile$CovB+trifile$CovC,VarNum=trifile$VarNum,VarSize=trifile$VarType))
  }
  if(dim(tetrafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=tetrafile$CovA+tetrafile$CovB+tetrafile$CovC+tetrafile$CovD,VarNum=tetrafile$VarNum,VarSize=tetrafile$VarType))
  }
  if(dim(pentafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=pentafile$CovA+pentafile$CovB+pentafile$CovC+pentafile$CovD+pentafile$CovE,VarNum=pentafile$VarNum,VarSize=pentafile$VarType))
  }
  
  fre_all <- NULL
  bifre <- c()
  trifre <- c()
  tetrafre <- c()
  pentafre <- c()
  
  if(dim(bifile)[1]!=0)
  {
    bifre <- apply(bifile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(bifre[1,],bifre[2,]),VarNum=rep(bifile$VarNum,time=2),VarSize=bifile$VarType))
    
  }
  if(dim(trifile)[1]!=0)
  {
    trifre <- apply(trifile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(trifre[1,],trifre[2,],trifre[3,]),VarNum=rep(trifile$VarNum,time=3),VarSize=trifile$VarType))
    
  }
  if(dim(tetrafile)[1]!=0)
  {
    tetrafre <- apply(tetrafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(tetrafre[1,],tetrafre[2,],tetrafre[3,],tetrafre[4,]),VarNum=rep(tetrafile$VarNum,time=4),VarSize=tetrafile$VarType))
    
  }
  if(dim(pentafile)[1]!=0)
  {
    pentafile <- apply(pentafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])+as.numeric(x[5])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum,as.numeric(x[5])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(pentafile[1,],pentafile[2,],pentafile[3,],pentafile[4,],pentafile[5,]),VarNum=rep(pentafile$VarNum,time=5),VarSize=pentafile$VarType))
    
  }
  return (list(cov,fre_all))
}

##### function: read multiple samples coverage files #####

colour.readcov <- function(file.prefix=x)
{
  bifile <- NULL
  trifile <- NULL
  tetrafile <- NULL
  pentafile <- NULL
  if(file.exists(paste0(file.prefix,"_bicov.txt"))){
    bifile <- read.table(paste0(file.prefix,"_bicov.txt"),
                         col.names = c("CovA","CovB","color","isStrict","VarType","VarId",
                                       "VarNum","coe","VarDis"),stringsAsFactors = FALSE,
                         colClasses = c("CovA"="numeric","CovB"="numeric"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_bicov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tricov.txt"))){
    trifile <- read.table(paste0(file.prefix,"_tricov.txt"),
                          col.names = c("CovA","CovB","CovC","color","isStrict","VarType","VarId",
                                        "VarNum","coe","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tricov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tetracov.txt"))){
    tetrafile <- read.table(paste0(file.prefix,"_tetracov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","color","isStrict","VarType","VarId",
                                          "VarNum","coe","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tetracov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_pentacov.txt"))){
    pentafile <- read.table(paste0(file.prefix,"_pentacov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","CovE","color","isStrict","VarType","VarId",
                                          "VarNum","coe","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_pentacov.txt")," ) does not exists !"))
  }
  
  
  
  cov <- NULL
  
  if(dim(bifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=bifile$CovA+bifile$CovB,VarNum=bifile$VarNum,coe=bifile$coe,color=bifile$color,VarSize=bifile$VarType))
  }
  if(dim(trifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=trifile$CovA+trifile$CovB+trifile$CovC,VarNum=trifile$VarNum,coe=trifile$coe,color=trifile$color,VarSize=trifile$VarType))
  }
  if(dim(tetrafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=tetrafile$CovA+tetrafile$CovB+tetrafile$CovC+tetrafile$CovD,VarNum=tetrafile$VarNum,coe=tetrafile$coe,color=tetrafile$color,VarSize=tetrafile$VarType))
  }
  if(dim(pentafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=pentafile$CovA+pentafile$CovB+pentafile$CovC+pentafile$CovD+pentafile$covE,VarNum=pentafile$VarNum,coe=pentafile$coe,color=pentafile$color,VarSize=pentafile$VarType))
  }  
  
  fre_all <- NULL
  bifre <- c()
  trifre <- c()
  tetrafre <- c()
  pentafre <- c()
  
  if(dim(bifile)[1]!=0)
  {
    bifre <- apply(bifile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(bifre[1,],bifre[2,]),VarNum=rep(bifile$VarNum,time=2),coe=rep(bifile$coe,time=2),color=rep(bifile$color,time=2),VarSize=bifile$VarType))
    
  }
  if(dim(trifile)[1]!=0)
  {
    trifre <- apply(trifile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(trifre[1,],trifre[2,],trifre[3,]),VarNum=rep(trifile$VarNum,time=3),coe=rep(trifile$coe,time=3),color=rep(trifile$color,time=3),VarSize=trifile$VarType))
    
  }
  if(dim(tetrafile)[1]!=0)
  {
    tetrafre <- apply(tetrafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(tetrafre[1,],tetrafre[2,],tetrafre[3,],tetrafre[4,]),VarNum=rep(tetrafile$VarNum,time=4),coe=rep(tetrafile$coe,time=4),color=rep(tetrafile$color,time=4),VarSize=tetrafile$VarType))
    
  }
  if(dim(pentafile)[1]!=0)
  {
    pentafre <- apply(pentafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])+as.numeric(x[5])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum,as.numeric(x[5])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(pentafre[1,],pentafre[2,],pentafre[3,],pentafre[4,],pentafre[5,]),VarNum=rep(pentafile$VarNum,time=5),coe=rep(pentafile$coe,time=5),color=rep(pentafile$color,time=5),VarSize=pentafile$VarType))
    
  }
  return (list(cov,fre_all))
  
}

#####SNJ17#####
SNJ17dir <- ""
list_read <- readcov(SNJ17dir) #39.3
cov <- 39.3
p <- 2
v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}


frequency <- list_read[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_num1 <- frequency[frequency$VarNum==1&frequency$VarSize<=10,]
frequency_num1$type <- "VarNum=1&VarSize<=10"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_num1)

frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


coverage <- list_read[[1]]
coverage$type <- "all"
coverage_num5 <- coverage[coverage$VarNum<=5&coverage$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
coverage_num1 <- coverage[coverage$VarNum==1&coverage$VarSize<=10,]
coverage_num1$type <- "VarNum=1&VarSize<=10"
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_num1)

coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=40),  
                            axis.text.y = element_text(size=40),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 30),
                            strip.text.y = element_text(size = 30),
                            legend.position='none')

p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov*(p-1),cov*(p+1))),linetype=4)+
  theme_bw()+facet_grid(~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=40),  
                          axis.text.y = element_text(size=40),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 30),
                          strip.text.y = element_text(size = 30),
                          legend.position='none')


all <- c(1.00005,0.356266,0.300164,0.215398,0.177853,0.150212,0.130034,0.114602,0.102419)
s11n6 <- c(0.990932,0.356502,0.298067,0.21587,0.177891,0.150271,0.130087,0.114653,0.102468)
s11n2 <- c(1.09362,0.362962,0.300708,0.2179,0.17947,0.151635,0.131299,0.115754,0.103485)
filter.method <- rep(c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"),each=9)
p <- 2
ploidy <- rep(2:10,times=3)
dt <- data.frame(avg.loglikelihood=c(all,s11n6,s11n2),filter=filter.method,ploidy=ploidy)
p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(color=color,yintercept=y),data.frame(y=c(all[p-1],s11n6[p-1],s11n2[p-1]),color=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10")),linetype=3)+
  theme_classic()+
  theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
        axis.title = element_text(size=50),
        axis.text.x = element_text(size=40),  
        axis.text.y = element_text(size=40),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40),
        legend.position = c(0.72,0.8),
        legend.title =element_text(size=40),
        legend.text =element_text(size=35),
        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))


#####LSX118
list_read <- read_cov("/data1/sunmz_data/trim_data/cyc/LS/LSX118/PloidyFrost_output/LSX118") #18
cov <- 18
p <- 4
v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}


frequency <- list_read[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_num1 <- frequency[frequency$VarNum==1&frequency$VarSize<=10,]
frequency_num1$type <- "VarNum=1&VarSize<=10"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_num1)

frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


coverage <- list_read[[1]]
coverage$type <- "all"
coverage_num5 <- coverage[coverage$VarNum<=5&coverage$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
coverage_num1 <- coverage[coverage$VarNum==1&coverage$VarSize<=10,]
coverage_num1$type <- "VarNum=1&VarSize<=10"
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_num1)

coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=40),  
                            axis.text.y = element_text(size=40),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 30),
                            strip.text.y = element_text(size = 30),
                            legend.position='none')

p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov*(p-1),cov*(p+1))),linetype=4)+
  theme_bw()+facet_grid(~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=40),  
                          axis.text.y = element_text(size=40),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 30),
                          strip.text.y = element_text(size = 30),
                          legend.position='none')

all <- c(0.227424,0.293669,0.363846,0.288452,0.283436,0.291667,0.315141,0.240532,0.237677)
s11n6 <- c(0.210516,0.295988,0.365156,0.279061,0.281015,0.315482,0.316616,0.237769,0.234488)
s11n2 <- c(0.219527,0.308364,0.380471,0.284347,0.290132,0.307597,0.330356,0.246293,0.243295)
filter.method <- rep(c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"),each=9)
p <- 4
ploidy <- rep(2:10,times=3)
dt <- data.frame(avg.loglikelihood=c(all,s11n6,s11n2),filter=filter.method,ploidy=ploidy)
p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(color=color,yintercept=y),data.frame(y=c(all[p-1],s11n6[p-1],s11n2[p-1]),color=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10")),linetype=3)+
  theme_classic()+
  theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
        axis.title = element_text(size=50),
        axis.text.x = element_text(size=40),  
        axis.text.y = element_text(size=40),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40),
        legend.position = c(0.72,0.8),
        legend.title =element_text(size=40),
        legend.text =element_text(size=35),
        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))

#####8ana

Fana <- readcov("/data1/sunmz_data/trim_data/ploidyfrost/8ana/PloidyFrost_output/Fana") #cov52
cov <- 52
p <- 8
v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}


frequency <- Fana[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_num1 <- frequency[frequency$VarNum==1&frequency$VarSize<=10,]
frequency_num1$type <- "VarNum=1&VarSize<=10"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_num1)

frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


coverage <- Fana[[1]]
coverage$type <- "all"
coverage_num5 <- coverage[coverage$VarNum<=5&coverage$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
coverage_num1 <- coverage[coverage$VarNum==1&coverage$VarSize<=10,]
coverage_num1$type <- "VarNum=1&VarSize<=10"
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_num1)

coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=30),  
                            axis.text.y = element_text(size=30),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 30),
                            strip.text.y = element_text(size = 30),
                            legend.position='none')

p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov*(p-1),cov*(p+1))),linetype=4)+
  theme_bw()+facet_grid(~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=30),  
                          axis.text.y = element_text(size=30),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 30),
                          strip.text.y = element_text(size = 30),
                          legend.position='none')

all <- c(0.16978,0.181245,0.248488,0.218651,0.245583,0.225221,0.274842,0.212538,0.205768)
s11n6 <- c(0.143511,0.165359,0.236592,0.205225,0.22885,0.211037,0.267179,0.197518,0.189292)
s11n2 <- c(0.181215,0.173444,0.272752,0.226054,0.268693,0.235357,0.305577,0.238818,0.229267)
avg.loglikelihood <- c(
  0.16978,0.181245,0.248488,0.218651,0.245583,0.225221,0.274842,0.212538,0.205768,
  0.143511,0.165359,0.236592,0.205225,0.22885,0.211037,0.267179,0.197518,0.189292,
  0.181215,0.173444,0.272752,0.226054,0.268693,0.235357,0.305577,0.238818,0.229267,
)
filter.method <- rep(c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"),each=9)
p <- 8
ploidy <- rep(2:10,times=3)
dt <- data.frame(avg.loglikelihood=avg.loglikelihood,filter=filter.method,ploidy=ploidy)
p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(color=color,yintercept=y),data.frame(y=c(all[p-1],s11n6[p-1],s11n2[p-1]),color=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10")),linetype=3)+
  theme_classic()+
  theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
        axis.title = element_text(size=50),
        axis.text.x = element_text(size=40),  
        axis.text.y = element_text(size=40),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40),
        legend.position = c(0.7,0.3),
        legend.title =element_text(size=40),
        legend.text =element_text(size=35),
        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))


#####SNJ#####
name <- c("SNJ17","SNJ18","SNJ110")
list_read <- colour.readcov("/data1/sunmz_data/trim_data/cyc/SNJ/SNJ/PloidyFrost_output/SNJ")
cov1 <- 39.3
cov2 <- 40
cov3 <- 45.8

p <- 2
kcov <- 200
v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}
frequency <- list_read[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_dayu <- frequency[frequency$coe>=0.25,]
frequency_dayu$type <- "Cramer's V >= 0.25"
frequency_xiaoyu <- frequency[frequency$coe<0.25,]
frequency_xiaoyu$type <- "Cramer's V < 0.25"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_dayu)
frequency <- rbind(frequency,frequency_xiaoyu)
frequency$color[frequency$color==0] <- "SNJ17"
frequency$color[frequency$color==1] <- "SNJ18"
frequency$color[frequency$color==2] <- "SNJ110"
frequency$color <- factor(frequency$color,levels=name)
frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25","Cramer's V < 0.25"))


coverage <- list_read[[1]]
coverage$type <- "all"

coverage_num5 <- coverage[coverage$VarNum<=5&coverage$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
coverage_dayu <- coverage[coverage$coe>=0.25,]
coverage_dayu$type <- "Cramer's V >= 0.25"
coverage_xiaoyu <- coverage[coverage$coe<0.25,]
coverage_xiaoyu$type <- "Cramer's V < 0.25"
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_dayu)
coverage <- rbind(coverage,coverage_xiaoyu)

dim(list_read[[1]][list_read[[1]]$color==0,])
dim(list_read[[1]][list_read[[1]]$color==1,])
dim(list_read[[1]][list_read[[1]]$color==2,])

dim(coverage_num5[coverage_num5$color==0,])
dim(coverage_num5[coverage_num5$color==1,])
dim(coverage_num5[coverage_num5$color==2,])

dim(coverage_dayu510[coverage_dayu510$color==0,])
dim(coverage_dayu510[coverage_dayu510$color==1,])
dim(coverage_dayu510[coverage_dayu510$color==2,])

dim(list_read[[1]][list_read[[1]]$color==0&(list_read[[1]]$coverage<cov1|list_read[[1]]$coverage>cov1*3),])
dim(list_read[[1]][list_read[[1]]$color==1&(list_read[[1]]$coverage<cov2|list_read[[1]]$coverage>cov2*3),])
dim(list_read[[1]][list_read[[1]]$color==2&(list_read[[1]]$coverage<cov3|list_read[[1]]$coverage>cov3*3),])

dim(coverage_xiaoyu[coverage_xiaoyu$color==0&(coverage_xiaoyu$coverage<cov1|coverage_xiaoyu$coverage>cov1*3),])
dim(coverage_xiaoyu[coverage_xiaoyu$color==1&(coverage_xiaoyu$coverage<cov2|coverage_xiaoyu$coverage>cov2*3),])
dim(coverage_xiaoyu[coverage_xiaoyu$color==2&(coverage_xiaoyu$coverage<cov3|coverage_xiaoyu$coverage>cov3*3),])

coverage_dayu510 <- list_read[[1]][list_read[[1]]$VarNum>5|list_read[[1]]$VarSize>10,]
dim(coverage_dayu510[coverage_dayu510$color==0&(coverage_dayu510$coverage<cov1|coverage_dayu510$coverage>cov1*3),])
dim(coverage_dayu510[coverage_dayu510$color==1&(coverage_dayu510$coverage<cov2|coverage_dayu510$coverage>cov2*3),])
dim(coverage_dayu510[coverage_dayu510$color==2&(coverage_dayu510$coverage<cov3|coverage_dayu510$coverage>cov3*3),])

coverage$color[coverage$color==0] <- "SNJ17"
coverage$color[coverage$color==1] <- "SNJ18"
coverage$color[coverage$color==2] <- "SNJ110"
coverage$color <- factor(coverage$color,levels=name)
coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25","Cramer's V < 0.25"))


p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(color~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=34),  
                            axis.text.y = element_text(size=28),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 28),
                            strip.text.y = element_text(size = 28),
                            legend.position='none')


p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov1*(p-1),cov1*(p+1),cov2*(p-1),cov2*(p+1),cov3*(p-1),cov3*(p+1)),color=rep(name,each=2)),linetype=4)+
  theme_bw()+facet_grid(color~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=34),  
                          axis.text.y = element_text(size=28),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 28),
                          strip.text.y = element_text(size = 28),
                          legend.position='none')

sample <- c("SNJ17","SNJ18","SNJ110")
s17 <- c(1.05808,0.360986,0.302128,0.217343,0.179301,0.151464,0.131149,0.11562,0.103364,
         1.04461,0.360759,0.300284,0.21741,0.179128,0.151333,0.131032,0.115513,0.103265,
         1.04415,0.362348,0.301035,0.218622,0.179922,0.152025,0.13165,0.11608,0.103794)

s18 <- c(1.07401,0.361564,0.304226,0.217727,0.179728,0.151821,0.131468,0.115912,0.103636,
         1.05997,0.361304,0.302159,0.21775,0.179504,0.151648,0.131314,0.115771,0.103506,
         1.0586,0.362517,0.302596,0.218636,0.180067,0.152143,0.131755,0.116176,0.103883)

s110 <- c(1.1237,0.361371,0.308822,0.21689,0.179666,0.151726,0.131384,0.115833,0.10356,
          1.1067,0.360687,0.307132,0.21685,0.179402,0.151521,0.131201,0.115666,0.103406,
          1.11396,0.362532,0.308665,0.218251,0.180375,0.152366,0.131955,0.116358,0.104051)

filter.method <- rep(rep(c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25"),each=9),times=3)
p <- 2
ploidy <- rep(2:10,times=9)
dt <- data.frame(avg.loglikelihood=c(s17,s18,s110),filter=filter.method,ploidy=ploidy,sample=rep(sample,each=27))
dt$sampe <- factor(sample,levels=name,labels=name)
p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(yintercept=y,color=filter),data.frame(y=c(s17[p-1],s17[p-1+9],s17[p-1+18],
                                                           s18[p-1],s18[p-1+9],s18[p-1+18],
                                                           s110[p-1],s110[p-1+9],s110[p-1+18]),
                                                       filter=rep(c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25"),times=3),
                                                       sample=factor(rep(sample,each=3),levels=name,labels=name)),
             linetype=3)+
  facet_grid(factor(sample,levels=name,labels=name)~.)+
  theme_classic()+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                        axis.title = element_text(size=50),
                        axis.text.x = element_text(size=34),  
                        axis.text.y = element_text(size=28),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        strip.text.x = element_text(size = 40),
                        strip.text.y = element_text(size = 40),
                        legend.position = c(0.7,0.9),
                        legend.title =element_text(size=40),
                        legend.text =element_text(size=35),
                        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))

