# rm(list = ls())
# gc(reset = T)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('../Evaluation.R')
source('../Forward_Learning_with_outlier.R')
source('../RGSEM_algorithm.R')

if(!require(tidyverse))install.packages('tidyverse');require(tidyverse)
if(!require(ggmap))install.packages('ggmap');require(ggmap)
if(!require(reshape2))install.packages('reshape2');require(reshape2)
if(!require(maptools))install.packages('maptools');require(maptools)
if(!require(data.table))install.packages('data.table');require(data.table)
if(!require(sp))install.packages('sp');require(sp)
if(!require(maptools))install.packages('maptools');require(maptools)

library("accelerometry")

korea_map_shp = readShapePoly('SHP/KR.shp')
korea_map = fortify(korea_map_shp)

Region_all = fread("Korea.csv",encoding = 'UTF-8')



covid = read.csv("covid_KR.csv")
covid = t(covid)
date_idx =row.names(covid)[-1]
date_idx = as.Date(paste0(substr(date_idx,2,5),'-' ,substr(date_idx,8,9),'-' ,substr(date_idx,12,13)))
covid[covid=='-']= 0
colnames(covid) = covid[1,] ; covid = covid[-1,] ;  covid = apply(covid,2, as.numeric)

covid_view = cbind(as.Date(date_idx),as.data.frame(covid))

# View(covid_view)


covid[covid<0]= 0
row.names(covid)= date_idx
covid_raw = covid
library("accelerometry")

roll = 3

for(j in 2:ncol(covid)){
  covid[roll:nrow(covid) ,j] = movingaves(covid[,j], window = roll)
}
covid = covid[-(1: (roll-1) ),]
date_idx = date_idx[c(-1,-length(date_idx))]

covid_korea = function(date_idx_for_analysis= date_idx_for_analysis ,alpha , C,restricted_infection, method,tikzedge = F){
  
  covid_for_analysis = covid[date_idx_for_analysis,]
  valid_region_idx = (colSums(covid_for_analysis)>restricted_infection)
  for(j in 1:ncol(covid_for_analysis)){
    covid_for_analysis[,j] = covid_for_analysis[,j]/sum(covid_for_analysis[,j])
  }
  
  Region = Region_all[valid_region_idx,]
  Region$address = as.character(Region$address)
  covid_for_analysis = covid_for_analysis[,valid_region_idx]
  covid_for_analysis[covid_for_analysis<0] = 0
  
  ###### Proposed#############
  
  if(sum(colnames(covid_for_analysis) == Region$address) != ncol(covid_for_analysis)){
    message("correct name of regions required.");
  }
  
  X =  covid_for_analysis
  # covid_new_temp  = covid_new_temp*100
  DAG_covid = GSEM_Algorithm(data = X, method = method, alpha = alpha, C = C)
  
  print(colnames(X)[DAG_covid$Ordering])
  
  DAG = DAG_covid$DAG#DAG=Ordering_covid$DAG
  DAG = as.data.frame(DAG)
  colnames(DAG) = rownames(DAG) = Region$address
  DAG = cbind(rownames(DAG), DAG) ; colnames(DAG)[1] = 'Address'
  
  direction = melt(DAG,id.vars = 'Address') %>% filter(value ==1) %>% dplyr::select(c(1,2))
  
  colnames(direction) = c('To', 'From')
  direction$From = as.character(direction$From);direction$To = as.character(direction$To)
  
  direction= direction %>% left_join(Region, by = c('To' = 'address' ))
  colnames(direction)[3:4] = c('To_Y', 'To_X')
  
  direction = direction %>% left_join(Region, by = c('From' = 'address' ))
  colnames(direction)[5:6] = c('From_Y', 'From_X')
  
  
  if(tikzedge == T){
    temp = direction[order(direction$From),c(2,1) ]
    
    temp2 = paste0("\\draw [edge] (",temp[,1],') -- (', temp[,2],"); \n") 
    for(i in 1:length(temp2)){
      cat(temp2[i])
    }
  }
  q= ggplot() + 
    geom_polygon(data = korea_map, aes(x = long, y = lat, group = group), fill = "#FFFFFF", colour = "#000000")+
    geom_segment(data = direction, aes(x = From_X, y= From_Y, xend = To_X, yend = To_Y, color= From),
                 arrow = arrow(length = unit(0.2,'cm'),type = 'closed', angle = 20) , size = 0.8)+
    geom_text(data= Region_all[valid_region_idx,] ,aes(x = Longitude, y= Latitude+0.1, label = address))+
    coord_fixed(ratio=1)
  
  return(q)
}



idx = (date_idx>='2020-03-04' &date_idx<='2020-06-04')

res_inf = 23
alpha = 0.05

result1 = covid_korea(idx , alpha =alpha , C = 6,restricted_infection = res_inf, method = 'Cols2', tikzedge = F) 
result2 = covid_korea(idx , alpha = alpha, C = 6,restricted_infection = res_inf, method = 'Cols1', tikzedge = F) 

