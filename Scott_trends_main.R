library(raster)
library(rgdal)
library(plyr)
library(stringr)
library(tidyverse)
library(readxl)
library(foreign)
library(raster)
library(here)
library(osfr)

# set to 1 if first time this code is run (HF processing, )
first <- 0
write <- 1

owd <- getwd()


if(first){
  
  trends <- read.csv(here("data/BBS_trend_eBird_all_spp.csv"),
                     stringsAsFactors=F)
  
  
  cand <- read.csv(here("data/ERD2016_PROD_species_list.csv"),
                   stringsAsFactors=F)
  names(cand)[2] <- "Species.Name"
  
  comb <- join(cand, trends, by = "Species.Name")
  
  cand2 <- read_excel(here("data/Species listing and status_Feb 2.xlsx"))
  names(cand2)[1] <- "Species.Name"
  
  comb <- join(comb, cand2, by = "Species.Name")
  
  comb <- comb[order(comb$SPECIES_CODE),]
  
  #keep only birds that match with comb list
  comb.red <- comb[!is.na(comb$Include) & comb$Include == 1,]
  #comb.red <- comb.red[!is.na(comb.red$Ebird.Analysis) & comb.red$Ebird.Analysis == 1,]
  
  nms <- comb.red$Species.Name
  
  #This is where the STEM models downloaded from OSF are
  
  STEM_dir <- "eBird_relative_abundance_data"
  if (file.exists(STEM_dir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(here(), STEM_dir))
    setwd(file.path(here(), STEM_dir))
    osf_retrieve_file("https://osf.io/3qhuv/") %>% 
      osf_download()
    fl <- list.files()
    unzip(fl)
    file.remove(fl)
    
    osf_retrieve_file("https://osf.io/qk73j/") %>% 
      osf_download()
  }
  
  basedir <- getwd()
  
  trim_r <- raster("srd2016_wh_raster_template_3xagg.tif")
  
  
  birds <- list.files()
  fl.splt <- str_split_fixed(birds, "-", 2)
  
  birds <- birds[fl.splt[,1] %in% comb.red$SPECIES_CODE]
  #View(cbind(comb.red,birds))
  
  bird.df <- data.frame(Species.Name=nms,eBird=birds, stringsAsFactors=F)
  bird.df <- join(bird.df,comb,by="Species.Name")
  
  fl.splt <- str_split_fixed(birds, "-", 2)
  
  foldstr <- "abundance_umean"
  nb.list <- list()
  
  #Dates
  dts <- c('01-04', '01-11', '01-18', '01-25', 
           '02-01', '02-08', '02-15', '02-22', '02-29', 
           '03-07', '03-14', '03-21', '03-28', 
           '04-04', '04-11', '04-18', '04-25', 
           '05-02', '05-09', '05-16', '05-23', '05-30', 
           '06-06', '06-13', '06-20', '06-27', 
           '07-04', '07-11', '07-18', '07-25', 
           '08-01', '08-08', '08-15', '08-22', '08-29', 
           '09-05', '09-12', '09-19', '09-26', 
           '10-03', '10-10', '10-17', '10-24', '10-31', 
           '11-07', '11-14', '11-21', '11-28', 
           '12-05', '12-12', '12-19', '12-26')
  
  #November 14‐March 14
  dts_red <- c('01-04', '01-11', '01-18', '01-25', 
               '02-01', '02-08', '02-15', '02-22', '02-29', 
               '03-07', '03-14', '11-14', '11-21', '11-28', 
               '12-05', '12-12', '12-19', '12-26')
  
  for (ii in 1:length(bird.df$eBird)){
    setwd(sprintf("%s/%s/%s",basedir,bird.df$eBird[ii],foldstr))
    
    print(bird.df$Species.Name[ii])
    flush.console()
    
    rast <- list.files()
    tt <- str_split(str_split(rast, "_week_", simplify = TRUE)[,2],
                    "_presentation", simplify = TRUE)[,1]
    
    jj <- 1
    for(rcand in rast[tt %in% dts_red]){
      tmp.r <- try(extend(raster(rcand), extent(trim_r)))
      
      if (class(tmp.r) == "RasterLayer"){
        if(jj == 1){
          spp.nbst <- tmp.r
          jj <- 0
        } else {
          spp.nbst <- stack(spp.nbst,tmp.r)
        }
        
      }
      
      rm(tmp.r)
      
    }
    
    nb.list[[ii]]<-mean(spp.nbst,na.rm=TRUE)
    names(nb.list[[ii]])<- sprintf("%s",bird.df$Species.Name[ii])
    # writeRaster(spp.b,file="spp.nonbreed.tif") # write all to working or new folder?
    
    setwd(basedir)
  }
  
  nb.list.red <- nb.list[sapply(nb.list, function(x) class(x) == "RasterLayer")]
  
  nb.stack <- stack(nb.list.red)
  
  save.image(here("/raster_setup.RData"))
  
} else{
  
  load("raster_setup.RData")
  first <- 0
  
}


nms.nb <- names(nb.stack)

if(first){
  loc <- here("data","rast_all_spp")
  for(ii in 1:length(nb.list.red)){
    writeRaster(nb.list.red[[ii]], filename= sprintf("%s/%s.tif", loc, nms.nb[ii]), format="GTiff")
  }
}

##########################
### Diversity index
##########################
loc2 <- here("out_all_spp")

if(first){
  library(vegan)
  nb.val <- values(nb.stack)

  nb.val[is.na(nb.val)] <- 0

  nb.div <- diversity(nb.val, index="shannon")

  nb.div.rast <- nb.list[[1]]

  nb.div.rast[] <- nb.div

  writeRaster(nb.div.rast, filename=paste(loc2,"/NB_Shannon.tif",sep=""), format="GTiff", overwrite=TRUE)
}

##########################
### HF 
##########################
# out_path <- here(paste0(owd, "/raster_other/")
HF2009 <- raster(here("data/raster_other","HF2009.tif"))
HF1993 <- raster(here("data/raster_other","HF1993.tif"))
HFchange <- HF2009 - HF1993

##########################
### Analysis part I: BBS trend maps
##########################
nb.stack.tmp <- nb.stack

comb.red2 <- comb.red


for(ii in 1:nlayers(nb.stack.tmp)){
  
  print(names(nb.stack.tmp[[ii]]))
  flush.console()
  
  #vv <- b.stack.tmp[[ii]][]
  xx <- nb.stack.tmp[[ii]][]
  
  #vv[!is.na(vv)] <- comb.red2$Annual.Trend..1966.2015.[ii]
  xx[!is.na(xx)] <- comb.red2$Annual.Trend..1966.2015.[ii]
  
  #b.stack.tmp[[ii]][] <- vv
  nb.stack.tmp[[ii]][] <- xx
  
  #  rm(vv,xx)
  rm(xx)
}

nb.mean <- mean(nb.stack.tmp, na.rm = TRUE)

writeRaster(nb.mean, filename=paste(loc2,"NB_Mean_trend.tif",sep="/"), format="GTiff", overwrite=TRUE)

nb.df <- data.frame(nb.stack.tmp[[1]][])
for(ii in 2:nlayers(nb.stack.tmp)){
  nb.df[,ii] <- nb.stack.tmp[[ii]][]
}
names(nb.df) <- names(nb.stack.tmp)

nb.median <- apply(nb.df,1, function(x) median(x, na.rm = TRUE))
qu <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)

nb.quant <- apply(nb.df,1, function(x) quantile(x, probs = qu, na.rm = TRUE))
nb.quant <- t(nb.quant)

nb.tmp <- nb.mean

nb.tmp[] <- nb.median
writeRaster(nb.tmp, filename=paste(loc2,"NB_Median_trend.tif",sep="/"), format="GTiff", overwrite=TRUE)

for(ii in 1:length(qu)){
  nb.tmp[] <- nb.quant[,ii]
  
  writeRaster(nb.tmp, filename=paste(loc2,"/NB_Quant_",qu[ii],"_trend.tif",sep=""), format="GTiff", overwrite=TRUE)
}


#what trend cutoff should we use
trend.cut <- -1

nb.df2 <- nb.df
nb.df2[is.na(nb.df2)] <- 0
nb.df2[nb.df2 > trend.cut] <- 0
nb.df2[nb.df2 <= trend.cut] <- 1

nb.cut.val <- rowSums(nb.df2, na.rm = TRUE)

nb.tmp[] <- nb.cut.val

writeRaster(nb.tmp, filename=paste(loc2,"NB_Sum_cutoff_minus_1.tif",sep="/"), format="GTiff", overwrite=TRUE)


##########################
### Analysis part I: PIF status
##########################
nb.stack.tmp.PIF <- nb.stack

comb.red2 <- comb.red

for(ii in 1:nlayers(nb.stack.tmp.PIF)){
  
  print(names(nb.stack.tmp.PIF[[ii]]))
  flush.console()
  
  xx <- nb.stack.tmp.PIF[[ii]][]
  
  xx[!is.na(xx)] <- comb.red2$`PIF status`[ii]
  
  nb.stack.tmp.PIF[[ii]][] <- xx
  
  rm(xx)
}

nb.mean.PIF <- mean(nb.stack.tmp.PIF, na.rm = TRUE)

writeRaster(nb.mean.PIF, filename=paste(loc2,"NB_Mean_PIF_score.tif",sep="/"), format="GTiff", overwrite=TRUE)


nb.df.PIF <- data.frame(nb.stack.tmp.PIF[[1]][])
for(ii in 2:nlayers(nb.stack.tmp.PIF)){
  nb.df.PIF[,ii] <- nb.stack.tmp.PIF[[ii]][]
}
names(nb.df.PIF) <- names(nb.stack.tmp.PIF)

nb.median.PIF <- apply(nb.df.PIF,1, function(x) median(x, na.rm = TRUE))

qu <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)

nb.quant.PIF <- apply(nb.df.PIF,1, function(x) quantile(x, probs = qu, na.rm = TRUE))
nb.quant.PIF <- t(nb.quant.PIF)

nb.tmp.PIF <- nb.mean.PIF

nb.tmp.PIF[] <- nb.median.PIF
writeRaster(nb.tmp.PIF, filename=paste(loc2,"NB_Median_PIF_score.tif",sep="/"), format="GTiff", overwrite=TRUE)

for(ii in 1:length(qu)){
  nb.tmp.PIF[] <- nb.quant.PIF[,ii]
  
  writeRaster(nb.tmp.PIF, filename=paste(loc2,"/NB_Quant_",qu[ii],"_PIF_score.tif",sep=""), format="GTiff", overwrite=TRUE)
  
}

##########################
### Analysis part II: HF change
### For part II we also need to re-run part I with just the 1993-2009 trends
##########################
nb.stack.tmp2 <- nb.stack

nb.HF.mean <- vector()

for(ii in 1:nlayers(nb.stack.tmp2)){
  
  print(names(nb.stack.tmp2[[ii]]))
  flush.console()
  
  xx <- nb.stack.tmp2[[ii]][]
  
  xx[!is.na(xx)] <- HFchange[!is.na(xx)]
  
  
  nb.HF.mean[ii] <- mean(xx, na.rm = TRUE)
  
  names(nb.HF.mean)[ii] <- names(nb.stack.tmp2[[ii]])
  
  rm(xx)
}

HF.fr <- data.frame(bird = gsub("B_","",names(nb.HF.mean)), 
                    HF_change_NB = nb.HF.mean,
                    BBS_trend_1996_2015 = comb.red2$Annual.Trend..1966.2015.)

write.csv(HF.fr, paste0(loc2,"/HF_change.csv"), row.names = FALSE)



##########################
### ESA landcover analysis
##########################
# memory.limit(200000)
# 
# in_path2 <- "D:/Work/LiberEro/_Project/ESA/"
# in_path3 <- "D:/Work/LiberEro/_Project/ESA/LC_detail_paper2/"
# out_path <- "D:/Work/LiberEro/_Project/_MS2_Scott_trends/Scott_trends/out_all_spp/"
# 
# ESA.leg <- read.csv(paste0(in_path2,"ESACCI-LC-Legend.csv"),sep=";")
# 
# ESA.raw <- raster(paste(in_path2,"ESACCI-LC-L4-LCCS-Map-300m_2015_clip.tif",sep=""))
# 
# if(first){
#   system(paste("gdalwarp -r near"
#                ,paste0(out_path,"NB_Shannon.tif")
#                ,paste0(in_path3,"diversity.tif"),sep=" "))
#   
#   system(paste("gdalwarp -r near"
#                ,paste0(out_path,"NB_Median_trend.tif")
#                ,paste0(in_path3,"BBS.tif"),sep=" "))
#   
#   system(paste("gdalwarp -r near"
#                ,paste0(out_path,"NB_Median_PIF_score.tif")
#                ,paste0(in_path3,"PIF.tif"),sep=" "))
#   
# }
# 
# Shannon <- raster(paste0(in_path3,"diversity.tif"))
# BBS <- raster(paste0(in_path3,"BBS.tif"))
# PIF <- raster(paste0(in_path3,"PIF.tif"))
# 
# LC <- stack(c(ESA.raw,Shannon,BBS,PIF))
# 
# out_l <- list()
# for(dd in 1:nlayers(LC)){
#   rr <- LC[[dd]]
#   nr <- nrow(rr)
#   nc <- ncol(rr)
#   #mm <- matrix(nrow = nr, ncol = nc)
#   cnt <- 1
#   for(ii in 1:nrow(rr)){
#     
#     if(!(ii %% 10000)){
#       vec_tmp <- getValues(rr, (ii-9999), 10000)
#       print(ii)
#       flush.console()
#       if(cnt == 1){
#         vec <- vec_tmp
#         cnt <- 2
#       } else{
#         vec <- c(vec, vec_tmp)
#       }
#     }
#     
#   }
#   vec <- c(vec,getValues(rr, (ii-ii %% 10000), (ii %% 10000)))
#   
#   out_l[[dd]] <- vec
#   names(out_l)[dd] <- names(rr)
# }
# 
# rm(list=ls()[! ls() %in% c("out_l","ESA.leg")])
# gc()
# 
# out_l_red <- list()
# for(ll in 1:2){
#   out_l_red[[ll]] <- out_l[[ll]][!is.na(out_l[[1]]) & out_l[[1]] != 210]
# }
# 
# names(out_l_red) <- c("NB_LAB", "idx")
# out_df <- as_data_frame(out_l_red)
# 
# outdf_red <- out_df[out_df$PIF > 0,]
# 
# Div_q <- quantile(outdf_red$Shannon, probs = c(0.2, 0.4, 0.6, 0.8))
# BBS_q <- quantile(outdf_red$BBS, probs = c(0.2, 0.4, 0.6, 0.8))
# PIF_q <- quantile(outdf_red$PIF, probs = c(0.2, 0.4, 0.6, 0.8))
# 
# lc_df <- outdf_red %>% group_by(NB_LAB) %>% 
#   dplyr::summarise(
#     cells = n(),
#     Div_low = sum(Shannon <= Div_q[1]),
#     Div_med_low = sum(Shannon > Div_q[1] & Shannon <= Div_q[2]),
#     Div_med = sum(Shannon > Div_q[2] & Shannon <= Div_q[3]),
#     Div_med_high = sum(Shannon > Div_q[3] & Shannon <= Div_q[4]),
#     Div_high = sum(Shannon > Div_q[4]),
#     
#     BBS_low = sum(BBS <= BBS_q[1]),
#     BBS_med_low = sum(BBS > BBS_q[1] & BBS <= BBS_q[2]),
#     BBS_med = sum(BBS > BBS_q[2] & BBS <= BBS_q[3]),
#     BBS_med_high = sum(BBS > BBS_q[3] & BBS <= BBS_q[4]),
#     BBS_high = sum(BBS > BBS_q[4]),
#     
#     PIF_low = sum(PIF <= PIF_q[1]),
#     PIF_med_low = sum(PIF > PIF_q[1] & PIF <= PIF_q[2]),
#     PIF_med = sum(PIF > PIF_q[2] & PIF <= PIF_q[3]),
#     PIF_med_high = sum(PIF > PIF_q[3] & PIF <= PIF_q[4]),
#     PIF_high = sum(PIF > PIF_q[4])
#     
#   )
# lc_df <- join(lc_df, ESA.leg, by="NB_LAB")
# 
# write.csv(lc_df, "./out_all_spp/LC_div_BBS_PIF_df.csv", row.names = FALSE)
# 
# tt <- lc_df[,c(2:7,18)]
# tt2 <- gather(tt, Class, value, -LCCOwnLabel, -cells)
# tt2$perc <- tt2$value/tt2$cells*100
# 
# ggplot(data=tt2, aes(x=LCCOwnLabel, y=perc, fill=Class)) +
#   geom_bar(stat="identity") + coord_flip()
# 
# 
# tt <- lc_df[,c(2,8:12,18)]
# tt2 <- gather(tt, Class, value, -LCCOwnLabel, -cells)
# tt2$perc <- tt2$value/tt2$cells*100
# 
# ggplot(data=tt2, aes(x=LCCOwnLabel, y=perc, fill=Class)) +
#   geom_bar(stat="identity") + coord_flip()
# 
# 
# tt <- lc_df[,c(2,13:18)]
# tt2 <- gather(tt, Class, value, -LCCOwnLabel, -cells)
# tt2$perc <- tt2$value/tt2$cells*100
# 
# ggplot(data=tt2, aes(x=LCCOwnLabel, y=perc, fill=Class)) +
#   geom_bar(stat="identity") + coord_flip()
# 
# ##########################
# ### Ecoregion analysis
# ##########################
# ecor_path <- "D:/Work/LiberEro/_Project/Ecoregions/rast/"
# out_path <- "D:/Work/LiberEro/_Project/_MS2_Scott_trends/Scott_trends/out_all_spp/"
# 
# ECO_ID <- read.dbf("D:/Work/LiberEro/_Project/Ecoregions/ecoregions_clipped_to_study_extent.dbf")
# ECO_ID_NAME <- ECO_ID %>% group_by(ECO_ID) %>% dplyr::summarise(ECO_NAME = first(ECO_NAME))
# 
# 
# shann <- raster(paste0(out_path,"NB_Shannon.tif"))
# BBS <- raster(paste0(out_path,"NB_Median_trend.tif"))             
# PIF <- raster(paste0(out_path,"NB_Median_PIF_score.tif"))
# 
# 
# 
# if(first){
#   ecor <- raster(paste0(ecor_path,"LE_ecoregions.tif"))
#   ecor[] <- NA
#   writeRaster(ecor, filename=paste0(ecor_path,"LE_ecoregions.tif"), format="GTiff", overwrite=TRUE)             
#   
#   system(paste("gdalwarp -r near"
#                ,paste0(ecor_path,"ecoregions_clipped_to_study_extent_diss_by_ECO_ID_proj.tif")
#                ,paste0(ecor_path,"LE_ecoregions.tif"),sep=" "))
# }
# 
# ecor <- raster(paste0(ecor_path,"LE_ecoregions.tif"))
# 
# stack_df <- as.data.frame(stack(c(ecor,shann,BBS,PIF)))
# stdf_red <-drop_na(stack_df)
# 
# names(stdf_red) <- c("ECO_ID", "Shannon", "BBS", "PIF")
# 
# 
# Div_q <- quantile(stdf_red$Shannon, probs = c(0.2, 0.4, 0.6, 0.8))
# BBS_q <- quantile(stdf_red$BBS, probs = c(0.2, 0.4, 0.6, 0.8))
# PIF_q <- quantile(stdf_red$PIF, probs = c(0.2, 0.4, 0.6, 0.8))
# 
# ecor_df <- stdf_red %>% group_by(ECO_ID) %>% 
#   dplyr::summarise(
#     cells = n(),
#     Div_low = sum(Shannon <= Div_q[1]),
#     Div_med_low = sum(Shannon > Div_q[1] & Shannon <= Div_q[2]),
#     Div_med = sum(Shannon > Div_q[2] & Shannon <= Div_q[3]),
#     Div_med_high = sum(Shannon > Div_q[3] & Shannon <= Div_q[4]),
#     Div_high = sum(Shannon > Div_q[4]),
#     
#     BBS_low = sum(BBS <= BBS_q[1]),
#     BBS_med_low = sum(BBS > BBS_q[1] & BBS <= BBS_q[2]),
#     BBS_med = sum(BBS > BBS_q[2] & BBS <= BBS_q[3]),
#     BBS_med_high = sum(BBS > BBS_q[3] & BBS <= BBS_q[4]),
#     BBS_high = sum(BBS > BBS_q[4]),
#     
#     PIF_low = sum(PIF <= PIF_q[1]),
#     PIF_med_low = sum(PIF > PIF_q[1] & PIF <= PIF_q[2]),
#     PIF_med = sum(PIF > PIF_q[2] & PIF <= PIF_q[3]),
#     PIF_med_high = sum(PIF > PIF_q[3] & PIF <= PIF_q[4]),
#     PIF_high = sum(PIF > PIF_q[4])
#     
#   )
# ecor_df <- join(ecor_df, ECO_ID_NAME, by="ECO_ID")
# 
# write.csv(ecor_df, "./out_all_spp/Ecoregion_df.csv", row.names = FALSE)
# 
# tt <- ecor_df[1:20,c(2,13:18)]
# tt2 <- gather(tt, Class, value, -ECO_NAME, -cells)
# tt2$perc <- tt2$value/tt2$cells*100
# 
# ggplot(data=tt2, aes(x=ECO_NAME, y=perc, fill=Class)) +
#   geom_bar(stat="identity") + coord_flip()


#####
# Comb metric
#####

shann <- raster(here("out_all_spp","NB_Shannon.tif"))
BBS <- raster(here("out_all_spp","NB_Median_trend.tif"))             
PIF <- raster(here("out_all_spp","NB_Median_PIF_score.tif"))
stack_df <- as.data.frame(stack(c(shann,BBS,PIF)))

names(stack_df) <- c("Shannon", "BBS", "PIF")

stack_df$Shannon[is.na(stack_df$BBS)] <- NA

stack_df$Shannon <- scale(stack_df$Shannon)
stack_df$BBS <- scale(stack_df$BBS)
stack_df$PIF <- scale(stack_df$PIF)
stack_df$PIF <- stack_df$PIF * -1


stack_df$Shannon <- stack_df$Shannon + abs(min(stack_df$Shannon, na.rm = TRUE))
stack_df$BBS <- stack_df$BBS + abs(min(stack_df$BBS, na.rm = TRUE))
stack_df$PIF <- stack_df$PIF + abs(min(stack_df$PIF, na.rm = TRUE))

stack_df$Shannon <- stack_df$Shannon/max(stack_df$Shannon, na.rm = TRUE)
stack_df$BBS <- stack_df$BBS/max(stack_df$BBS, na.rm = TRUE)
stack_df$PIF <- stack_df$PIF/max(stack_df$PIF, na.rm = TRUE)

hist(stack_df$Shannon)
hist(stack_df$BBS)
hist(stack_df$PIF)

stack_df$comb_score <- 3 * (stack_df$Shannon * stack_df$BBS * stack_df$PIF) / 
  (stack_df$Shannon + stack_df$BBS + stack_df$PIF)

hist(stack_df$comb_score)
comb <- shann
comb[] <- stack_df$comb_score
writeRaster(comb, filename=here("out_all_spp","Sh_BBC_PIF_comb.tif"), format="GTiff", overwrite=TRUE)

#####
# Comb metric Shannon BBS
#####


stack_df$BBS <- (stack_df$BBS - 1) * -1

hist(stack_df$Shannon)
hist(stack_df$BBS)
hist(stack_df$PIF)

stack_df$comb_score_Div_BBS <- 2 * (stack_df$Shannon * stack_df$BBS) / 
  (stack_df$Shannon + stack_df$BBS)

hist(stack_df$comb_score_Div_BBS)

stack_df$sqrt_Div_BBS <- sqrt(stack_df$Shannon * stack_df$BBS) 
hist(stack_df$sqrt_Div_BBS)

comb <- shann
comb[] <- stack_df$sqrt_Div_BBS
writeRaster(comb, filename=here("out_all_spp","Sh_BBS_neg_sqrt.tif"), format="GTiff", overwrite=TRUE)

##########################
### setup csv for Scott
##########################


########
#diversity
library(vegan)
nb.val <- values(nb.stack)
nb.val[is.na(nb.val)] <- 0
nb.div <- diversity(nb.val, index="shannon")

########
#median trend
nb.stack.tmp <- nb.stack
comb.red2 <- comb.red

for(ii in 1:nlayers(nb.stack.tmp)){
  xx <- nb.stack.tmp[[ii]][]
  xx[!is.na(xx)] <- comb.red2$Annual.Trend..1966.2015.[ii]
  nb.stack.tmp[[ii]][] <- xx
  rm(xx)
}

nb.mean <- mean(nb.stack.tmp, na.rm = TRUE)

nb.df <- data.frame(nb.stack.tmp[[1]][])
for(ii in 2:nlayers(nb.stack.tmp)){
  nb.df[,ii] <- nb.stack.tmp[[ii]][]
}
names(nb.df) <- names(nb.stack.tmp)

nb.median <- apply(nb.df,1, function(x) median(x, na.rm = TRUE))


ss <- as.data.frame(stack(nb.df))
ss2 <- nb.df
ss2[is.na(ss2)] <- 0

neg.spp <- apply(ss2, 1, function(row)  sum(row < 0))
n.spp <- apply(ss2, 1, function(row)  sum(row != 0))

########
### HF
out_path <- paste0(owd, "/raster_other/")

HF93 <- raster("D:/Work/LiberEro/_Project/HF/Maps/HFP1993.tif")
HF09 <- raster("D:/Work/LiberEro/_Project/HF/Maps/HFP2009.tif")

if(first){
  
  HF.delta <- HF09 - HF93
  
  writeRaster(HF.delta, filename=paste0(out_path,"HF_delta_1km.tif"), format="GTiff", overwrite=TRUE)
  
  
  system(paste("gdalwarp -r average"
               ,paste0(out_path,"HF_delta_1km.tif")
               ,paste0(out_path,"HF_change_mean.tif"),sep=" "))
  
  system(paste("gdalwarp -r med"
               ,paste0(out_path,"HF_delta_1km.tif")
               ,paste0(out_path,"HF_change_median.tif"),sep=" "))
  
}

HF.delta.mean <- raster(paste0(out_path,"HF_change_mean.tif"))[]
HF.delta.mean[HF.delta.mean < -100] <- 0
HF.delta.median <- raster(paste0(out_path,"HF_change_median.tif"))[]
HF.delta.median[HF.delta.median < -100] <- 0

gc()
########
#LC
# memory.limit(200000)
# in_path2 <- "D:/Work/LiberEro/_Project/ESA/"
# ESA.raw <- raster(paste(in_path2,"ESACCI-LC-L4-LCCS-Map-300m_2015_clip.tif",sep=""))
# 
# idx.r <- nb.stack[[1]]
# idx.r[] <- 1:length(nb.stack[[1]][])
# 
# if(first){
#   writeRaster(idx.r, filename=paste0(out_path,"Index_raster.tif"), format="GTiff", overwrite=TRUE)
#   
#   
#   system(paste("gdalwarp -r near"
#                ,paste0(out_path,"Index_raster.tif")
#                ,paste0(in_path2,"Index_to_ESA_clip_float.tif"),sep=" "))
# }
# 
# ESA.idx <- raster(paste(in_path2,"Index_to_ESA_clip_float.tif",sep=""))
# 
# 
# LC <- stack(c(ESA.raw,ESA.idx))
# 
# out_l <- list()
# for(dd in 1:nlayers(LC)){
#   rr <- LC[[dd]]
#   nr <- nrow(rr)
#   nc <- ncol(rr)
#   mm <- matrix(nrow = nr, ncol = nc)
#   cnt <- 1
#   for(ii in 1:nrow(rr)){
#     
#     if(!(ii %% 10000)){
#       vec_tmp <- getValues(rr, (ii-9999), 10000)
#       print(ii)
#       flush.console()
#       if(cnt == 1){
#         vec <- vec_tmp
#         cnt <- 2
#       } else{
#         vec <- c(vec, vec_tmp)
#       }
#     }
#     
#   }
#   vec <- c(vec,getValues(rr, (ii-ii %% 10000), (ii %% 10000)))
#   
#   out_l[[dd]] <- vec
#   names(out_l)[dd] <- names(rr)
# }
# 
# #rm(list=ls()[! ls() %in% c("out_l","ESA.leg")])
# gc()
# 
# out_l_red <- list()
# for(ll in 1:2){
#   out_l_red[[ll]] <- out_l[[ll]][!is.na(out_l[[1]]) & out_l[[1]] != 210]
# }
# 
# names(out_l_red) <- c("NB_LAB", "idx")
# out_df <- as_data_frame(out_l_red)
# 
# sum_out_df <- out_df %>% group_by(idx, NB_LAB) %>% summarise(count = n())
# 
# out_df_spread <- sum_out_df %>% spread(NB_LAB, count, fill = 0)
# 
# out_df_spread <- out_df_spread[out_df_spread$idx > 220,]

########
### WDPA
########
#library(foreign)

#WDPA.med <- data.frame(Value=raster(paste0(out_path,"WDPA.tif"))[])
#WDPA.df <- read.dbf("raster_other/WDPA_1km.tif.vat.dbf", as.is = T)
#WDPA.med <- join(WDPA.med, WDPA.df, by = "Value")

if(first){
  WDPA.r <- raster(paste0(out_path,"WDPA.tif"))
  levels(WDPA.r)<- WDPA.df$IUCN_CAT
  r_list <- list()
  
  for(ii in 1:7){
    tmp.r <- WDPA.r
    tmp.r[] <- NA #mask(WDPA.r, WDPA.r != ii, maskvalue=1)
    writeRaster(tmp.r, filename=paste0(out_path,"WDPA_", WDPA.df$IUCN_CAT[ii], ".tif"), format="GTiff", overwrite=TRUE)
    
    system(paste("gdalwarp -r average"
                 ,paste0(out_path,"WDPA_", WDPA.df$IUCN_CAT[ii], "_1km.tif")
                 ,paste0(out_path,"WDPA_", WDPA.df$IUCN_CAT[ii], ".tif"),sep=" "))
    rm(tmp.r)
  }
}

WDPA_Ia <- raster(paste0(out_path,"WDPA_Ia.tif"))
WDPA_Ib <- raster(paste0(out_path,"WDPA_Ib.tif"))
WDPA_II <- raster(paste0(out_path,"WDPA_II.tif"))
WDPA_III <- raster(paste0(out_path,"WDPA_III.tif"))
WDPA_IV <- raster(paste0(out_path,"WDPA_IV.tif"))
WDPA_V <- raster(paste0(out_path,"WDPA_V.tif"))
WDPA_VI <- raster(paste0(out_path,"WDPA_VI.tif"))

#rr <- raster(paste0(out_path,"WDPA_", WDPA.df$IUCN_CAT[ii], "_1km.tif"))
#zz <- raster(paste0(out_path,"Index_raster.tif"))
#rr <- try(extend(rr, zz))
#dd <- zonal(rr, zz, 'sum', na.rm = TRUE)

########
### create csv.df

csv.r <- nb.stack[[1]]
csv.r[] <- 1:length(nb.stack[[1]][])
r.pts <- rasterToPoints(csv.r, spatial=TRUE)
names(r.pts) <- "idx"

csv.df <- data.frame(idx = r.pts@data,
                     long = coordinates(r.pts)[,1],
                     lat = coordinates(r.pts)[,2],
                     diversity = nb.div,
                     median.decl = nb.median,
                     spp.decl = neg.spp,
                     n.spp = n.spp,
                     HF.delta.mean,
                     HF.delta.median,
                     HF.93 = HF1993[],
                     HF.09 = HF2009[],
                     WDPA_Ia = WDPA_Ia[],
                     WDPA_Ib = WDPA_Ib[],
                     WDPA_II = WDPA_II[],
                     WDPA_III = WDPA_III[],
                     WDPA_VI = WDPA_IV[],
                     WDPA_V =WDPA_V[],
                     WDPA_VI = WDPA_VI[]
)#,                
#ss)

csv.df <- join(csv.df, out_df_spread, by = "idx")

csv.df.red <- csv.df[!is.na(csv.df$median.decl),]

csv.df.red$flag <- 1

write.csv(csv.df.red, "csv.df.WDPA.csv", row.names = F)

idx.df <- data.frame(idx = idx.r[])
idx.df <- join(idx.df, csv.df.red, by = "idx")

idx.r2 <- idx.r
idx.r2[] <- idx.df$flag

idx.r3 <- idx.r
csv.df2 <- csv.df
names(csv.df2)[1] <- "ID"
levels(idx.r3) <- csv.df2

writeRaster(idx.r2, filename=paste0(out_path,"Index_raster_flag.tif"), format="GTiff", overwrite=TRUE)
#not run as this takes forever
#writeRaster(idx.r3, filename=paste0(out_path,"Index_raster_atr_table.tif"), format="GTiff", overwrite=TRUE)

# csv.df analysis
library(tidyverse)
csv.df <- read.csv("csv.df.csv")

str(csv.df)

#keep only rows that represent the top 10% of diversity values
cutoff <- quantile(csv.df$diversity, probs = c(0.25, 0.5, 0.75, 0.9))
csv.df.red <- csv.df[csv.df$diversity > cutoff[4],]

#not sure how useful this is, but here the setup to group by a variable and create summaries
test <- csv.df.red %>% group_by(spp.decl) %>% summarise(count = n(),
                                                        med.decl = mean(median.decl))
print(test, n = nrow(test))

####################################################################################
### Land cover change
####################################################################################
first <- 0

lc_in <- "D:/Work/LiberEro/_Project/_MS2_Scott_trends/Scott_trends/land_use/SSPs_may2017/"
out_path <- "D:/Work/LiberEro/_Project/_MS2_Scott_trends/Scott_trends/raster_other/"
setwd(lc_in)

if(first){
  fls <- list.files(pattern = ".asc")
  tt <- stack(x = fls)
  pr <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  crs(tt) <- pr
  
  writeRaster(tt, filename= ".tif", bylayer = TRUE, suffix = "names", format="GTiff")

  fls <- list.files(pattern = ".tif?")
  
  for(ii in fls){
    system(paste("gdalwarp -r near",
                 paste0(lc_in,ii), 
                 paste0(out_path,ii),sep=" "))
  }
}

setwd(out_path)
fls <- c("ssp1_year_50.tif", "ssp2_year_50.tif", "ssp3_year_50.tif", "year_2000.tif")
tt2 <- as.data.frame(stack(fls))
idx <- as.data.frame(raster("Index_raster.tif"))

lc.df <- data.frame(idx = idx, tt2)
names(lc.df)[1] <- "idx"

setwd(owd)
csv.df <- read.csv("csv.df.csv")

csv.df.lc <- join(csv.df[,1:11], lc.df, by = "idx")

write.csv(csv.df.lc, "csv.df.land_change.csv", row.names = FALSE)

####################################################################################
### Figures
####################################################################################
library(raster)
library(readxl)
library(dplyr)
library(plyr)
library(tidyverse)
library(sf)
library(raster)
library(RColorBrewer)
library(here)
library(leaflet)
library(png)
library(fields)

load("raster_setup.RData")

here("R") %>% 
  list.files(full.names = TRUE) %>% 
  walk(source)

# natural earth political boundaries
ne_land <- read_sf("data/ne-land.gpkg") %>% st_geometry()
ne_country_lines <- read_sf("data/ne-country-lines.gpkg") %>% st_geometry()
ne_state_lines <- read_sf("data/ne-state-lines.gpkg") %>% st_geometry()


idx.r <- nb.stack[[1]]
idx.r[] <- 1:length(nb.stack[[1]][])

figs <- read_xlsx("Paper 2 data for figures_updated.xlsx")

figs2 <- mutate(figs,
       Fig1a = ifelse(Fig1a.Diversity.high.priority == 1, 3, ifelse(Fig1a.Diversity.mod.priority == 1, 2, 1)),
       Fig1b = ifelse(Fig1b.Med.declines.high.priority == 1, 3, ifelse(Fig1b.Med.declines.mod.priority  == 1, 2, 1)),
       Fig1c = ifelse(Fig1c.No.declining.high.priority == 1, 3, ifelse(Fig1c.No.declining.mod.priority   == 1, 2, 1)),
       
       Fig2a = ifelse(Fig2a.Overlap.diversity.med.decline.high == 1, 3, ifelse(Fig2a.Overlap.diversity.med.decline.mod == 1, 2, 1)),
       Fig2b = ifelse(Fig2b.Overlap.diversity.no.decline.high == 1, 3, ifelse(Fig2b.Overlap.diversity.no.decline.mod == 1, 2, 1)),
       Fig2c = ifelse(Fig2c.Overlap.med.decline.no.decline.high == 1, 3, ifelse(Fig2c.Overlap.med.decline.no.decline.mod   == 1, 2, 1))
)

idx.df <- data.frame(idx = idx.r[])
idx.df <- join(idx.df, figs2, by = "idx")

fig.r <- idx.r
fig.r[] <- idx.df$Fig1a
plot(fig.r)

fig.lst <- list()

for(ii in 19:24){
  fig.r[] <- idx.df[,ii]
  fig.lst[[ii- 18]] <- fig.r
}

fig.st <- stack(fig.lst)
names(fig.st) <- c("Fig1a", "Fig1b", "Fig1c", "Fig2a", "Fig2b", "Fig2c")


# process for visualization
crs <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
abd_plot <- fig.st %>% 
  stem_to_na() %>% 
  projectRaster(crs = crs) %>% 
  sqrt() %>% 
  stem_crop()
# remove migration cells covered by abundance in another season
#abd_plot[["migration"]] <- mask(abd_plot[["migration"]],
#                                calc(abd_plot[[c("breeding", "winter")]], 
#                                     mean, na.rm = TRUE),
#                                inverse = TRUE)


# plot
e <- extent(abd_plot)
text_col <- "black"

fig <- names(abd_plot)
palette <- c("Greens", "Blues", "BuPu", "Reds")
legend_offsets <- c(0.01, 0.06, 0.11, 0.16)
# prepare vector layers
land <- st_transform(ne_land, crs = crs)
country <- st_transform(ne_country_lines, crs = crs)
state <- st_transform(ne_state_lines, crs = crs)

# plot fig
for (ii in seq_along(fig)) {
  
  title <- paste(fig[ii])
  # print map
  here("out_all_spp/fig", paste0(fig[ii], ".png")) %>% 
    png(width = 3000, height = 3000, res = 300)
  par(mar = c(0, 0, 0, 0), oma = c(0,0,0,0), bg = "white")
  plot(land, col = "grey40", border = NA, xlim = e[1:2], ylim = e[3:4])
  
  pal <- brewer.pal(3, "Set1")
  pal <- colorRampPalette(pal)
  #pal <- colorFactor("RdYlBu", domain = NULL)
  #pal <- colorQuantile(pal, values(abd_plot[[fig[ii]]]), n = 8, #probs = seq(0, 1, length.out = n + 1),
  #              na.color = "#808080", alpha = FALSE, reverse = FALSE)
  plot(abd_plot[[gsub(" ", ".",fig[ii])]], add = TRUE, col = pal(256), legend = FALSE, 
       maxpixels = ncell(abd_plot))
  add_legend(fig[ii], pal, legend_offsets[3], low_high = TRUE,
             text_col = text_col)
  
  # boundaries
  plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
  plot(country, col = "black", lwd = 1, add = TRUE)
  
  # title
  # plot bounds
  usr <- par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]
  # labels
  text(x = usr[1] + 0.05 * xwidth, y = usr[3] + 0.21 * yheight,
       labels = title, pos = 4, font = 1, cex = 1.5, col = text_col)
  
  #rasterImage(logo,usr[1] + 0.01 * xwidth, usr[3] + 0.03 * yheight,
  #            usr[1] + 0.38 * xwidth, usr[3] + 0.09 * yheight)
  dev.off()
}






