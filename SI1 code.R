#Supplementary Information 1.
#Multi-taxa ecological responses to habitat loss and fragmentation in western Amazonia 

#1)	R script to calculate ClassStat metrics of SDMtools package.

##Install packages##
install.packages(c("foreach", "tidyverse", "doSNOW", "raster",
                   "rgdal", "SDMTools", "sp", "rgeos" ))

library(foreach)
library(tidyverse)
library(doSNOW)
library(raster)
library(rgdal)
library(rgeos)
library(sp)


####Inputs#### 

##Upload raster file## 
#(Cobertura Natural = clase 1 --- Intervenido = Clase 0 --- Otros = NoData)
Raster <- raster("Inputs/MapaBinario_rst/CobSueloFnF.tif")

##Upload rasterized grid  
ResolutionMask<-raster("Inputs/ResolutionMask_rst/MaskRSL.tif")

##Upload grid shapefile##
shp_files <- list.files("Inputs/Malla_shp/", ".shp$", full.names = TRUE)[1] 

## Upload -FragIndex- ##
FragIndex_shp <- shapefile("Inputs/FragIndex_shp/FragIndex.shp")

##Upload plot location##
par <- shapefile("Inputs/Parcelas_shp/Parcelas.shp")

##Run classStatt function##
source("Scripts/FunsClassMetrics.R")

####Paralel processing t####

##Configurar R para ejecutar calculo en paralelo##
#el numero de nucleos debera ser modificado dependiendo de la  capacidad de procesamiento del computador#
cluster <- makeCluster(3)
registerDoSNOW(cluster)

##Proceso en paralelo##
foreach(shp = shp_files, 
        .packages = c("raster", "tidyverse", "foreach")) %do% {
          #vectores
          Sectors <- shapefile(shp)
          
          # Analysis  ----                   
          out <- foreach(i = seq_along(Sectors)) %dopar%
            frag_by_sector(Raster, Sectors[i,], ID = "id", 
                           Background = 0,
                           output = paste0("Outputs/csv/",
                                           gsub(".shp", ".csv",basename(shp))),
                           latlon = FALSE) %>%
            bind_rows()
          
          # Join by atributes
          Sectors@data <- Sectors@data %>%
            left_join(out, by = c("id"="OBJECTID"))
          
          
          
          #guardar datos a shapefile y/o csv
          shapefile(Sectors, filename=paste0("Outputs/shp/", 
                                             basename(shp)), overwrite=TRUE)
        }

##Parar Cluster al terminar el proceso paralelo## 
#Si el proceso se interrumpe se debe parar el cluster para luego volverlo a registrar#
stopCluster(cluster)  

####Ajuste y seleccion de metricas####

teselas_shp  <- readOGR("Outputs/shp", 
                        strsplit(basename(shp_files), split = "[.]")[[1]][1])

teselas_shp@data  <- teselas_shp@data %>% 
  mutate(prp_l = ifelse(is.nan(prp_l), 0, prp_l),
         ptch_ = ifelse(prp_l == 1, 0, ptch_),
         ttl_d =ifelse(prp_l == 1, 0, ttl_d),
         edg_d =ifelse(prp_l == 1, NA, edg_d),
         lnd__ =ifelse(prp_l == 1, 1, lnd__),
         lrg__ =ifelse(prp_l == 1, 1, lrg__),
         men_ptch_r =ifelse(near(prp_l, 1), ttl_r, men_ptch_r),
         sd_p_ =ifelse(prp_l == 1, NA, sd_p_),
         min_ptch_r =ifelse(prp_l == 1, ttl_r, min_ptch_r),
         mx_p_ =ifelse(prp_l == 1, NA, mx_p_),
         pr___ =ifelse(prp_l == 1, NA, pr___),
         men_prm_r_ =ifelse(prp_l == 1, NA, men_prm_r_),
         sd_pr__ =ifelse(prp_l == 1, NA, sd_pr__),
         min_prm_r_ =ifelse(prp_l == 1, NA, min_prm_r_),
         mx_pr__ =ifelse(prp_l == 1, NA, mx_pr__),
         men_shp_nd =ifelse(prp_l == 1, 1, men_shp_nd),
         sd_s_ =ifelse(prp_l == 1, NA, sd_s_),
         min_shp_nd =ifelse(prp_l == 1, NA, min_shp_nd),
         mx_s_ =ifelse(prp_l == 1, NA, mx_s_),
         men_frc_dm=ifelse(prp_l == 1, NA, men_frc_dm),
         sd_f__=ifelse(prp_l == 1, NA, sd_f__),
         min_frc_dm=ifelse(prp_l == 1, NA, min_frc_dm),
         mx_f__=ifelse(prp_l == 1, NA, mx_f__),
         ttl__=ifelse(prp_l == 1, NA, ttl__),
         prp_ln_=ifelse(prp_l == 1, 1, prp_ln_),
         men_ptch_c=ifelse(prp_l == 1, NA, men_ptch_c),
         sd_pt__=ifelse(prp_l == 1, NA, sd_pt__),
         min_ptch_c=ifelse(prp_l == 1, NA, min_ptch_c),
         mx_pt__=ifelse(prp_l == 1, NA, mx_pt__),
         prp_lk_=ifelse(prp_l == 1, 1, prp_lk_),
         aggr_=ifelse(prp_l == 1, NA, aggr_),
         lns__=ifelse(prp_l == 1, 0, lns__),
         splt_=ifelse(prp_l == 1, NA, splt_),
         eff__=ifelse(prp_l == 1, NA, eff__),
         ptc__=ifelse(prp_l == 1, NA, ptc__)) %>%  
  dplyr::select(prp_l, lnd__,lrg__, min_ptch_r, men_ptch_r,men_shp_nd, lns__,prp_lk_,  prp_ln_)  
# Variables seleccionadas para el proyecto piloto#

####Incorporacion de datos de parcela####
#spatial join#
teselas_tmp <- teselas_shp
teselas_tmp@data <- teselas_tmp@data %>% mutate(ID = 1:n())
overlay <- over(teselas_tmp, par) %>% mutate(ID = as.character(1:n())) %>% 
  filter(!is.na(Id))


####Rasterizacion de varibles seleccionadas####


##Ajuste de resoluci?n del mapa de cobertura del suelo##  
Resample <- resample(Raster, ResolutionMask, method='ngb')
values(Resample)[!is.na(values(Resample))] <- 1

##Rasterizacion##

colNames <- names(teselas_shp)
Rst_teselas <- lapply(colNames, function(x)  rasterize(teselas_shp, Resample, field = x ) )
rst_mult <- lapply(Rst_teselas, function(x, rstIN){ x*rstIN}, Resample )
names(rst_mult) <- colNames
oldgw <- getwd()
setwd("Outputs/Raster/")
#exporta lista de rasters
for(i in seq(length(colNames))) { 
  writeRaster(rst_mult[[i]], filename = paste0(colNames[i], ".tif"), 
              format= "GTiff", bylayer=T, 
              suffix='names', overwrite=T) 
}
setwd(oldgw) 

setwd("Outputs/Raster/")
FragIndex_rst <- rasterize(FragIndex_shp, ResolutionMask, field=FragIndex_shp$FragIndex)

writeRaster(FragIndex_rst, filename="FragIndex.tif", format="GTiff", overwrite=TRUE)

setwd(oldgw) 

rst_mult[[ length(colNames) + 1  ]] <-  FragIndex_rst
colNames[ length(colNames) +1 ] <- "FragIndex"

names(rst_mult) <- colNames


####Extraccion y transferencia de datos a centroides####
#creacion de centroides#
centroids <- gCentroid(teselas_shp, byid = T)
#Estraccion de datos en raster y calculo de coordenadas#
centroidsData <- lapply(rst_mult,function(x, y){SDMTools::extract.data(y, x) }, 
                        coordinates(centroids) ) %>% bind_cols %>% 
  mutate(Xcoord = coordinates(centroids)[, 1],
         Ycoord = coordinates(centroids)[, 2], 
         ID = as.character(1:n()) ) %>% left_join(overlay, by = "ID") %>% 
  select(-ID, -Id)

####Output .csv####
##sobrescribir csv##
write.csv(centroidsData, 
          paste0("Outputs/csv/",strsplit(basename(shp_files), 
                                         split = "[.]")[[1]][1], ".csv") ,na = "" )


#############################################################################################

###################### Class metrics function#############
frag_by_sector <- function(rst = Raster, Sector = test, 
                           ID = "id",
                           cs = res(rst)[1], 
                           Background = 0, 
                           output = NA, ...){
  library(SDMTools)
  library(tidyverse)
  library(raster)
  selected = Sector
  rst_msk = rst %>% crop(selected) %>% mask(selected)
  cl = ClassStat(rst_msk, cellsize = cs, bkgd = Background, ...)
  #cl = ClassStat(rst_msk, cellsize = cs, bkgd = Background)
  
  if(!is.null(cl)) {
    cl$OBJECTID = selected@data[1,ID]
    
  } else {
    cl = data.frame("class" = 1, 
                    "n.patches" = NA, 
                    "total.area" = round(sum(rst_msk[rst_msk != Background])*cs^2, 4),
                    "prop.landscape" = round((cellStats(rst_msk, sum)*cs^2)/
                                               (sum(rst_msk[rst_msk != Background])*cs^2),4),
                    "patch.density" = NA, 
                    "total.edge" = NA,
                    "edge.density" = NA, 
                    "landscape.shape.index" = 0,
                    "largest.patch.index" = 0, 
                    "mean.patch.area" = 0,
                    "sd.patch.area" = NA, 
                    "min.patch.area" = 0, 
                    "max.patch.area" = NA,
                    "perimeter.area.frac.dim" = NA, 
                    "mean.perim.area.ratio" = NA,
                    "sd.perim.area.ratio" = NA, 
                    "min.perim.area.ratio" = NA,
                    "max.perim.area.ratio" = NA, 
                    "mean.shape.index" = 0,
                    "sd.shape.index" = NA, 
                    "min.shape.index" = NA, 
                    "max.shape.index" = NA,
                    "mean.frac.dim.index" = NA, 
                    "sd.frac.dim.index" = NA, 
                    "min.frac.dim.index" = NA,
                    "max.frac.dim.index" = NA, 
                    "total.core.area" = NA, 
                    "prop.landscape.core" = 0,
                    "mean.patch.core.area" = NA, 
                    "sd.patch.core.area" = NA, 
                    "min.patch.core.area" = NA,
                    "max.patch.core.area" = NA, 
                    "prop.like.adjacencies" = 0,
                    "aggregation.index" = NA, 
                    "lanscape.division.index" = 1,
                    "splitting.index" = NA, 
                    "effective.mesh.size" = NA, 
                    "patch.cohesion.index" = NA,
                    "OBJECTID" = selected@data[1,ID])
    
  }
  if(!is.null(output)) {
    if(!file.exists(output)) {
      write_csv(cl, output)
    } else {
      write_csv(cl, output, append = TRUE, col_names = FALSE)}
  }
  cl
}

