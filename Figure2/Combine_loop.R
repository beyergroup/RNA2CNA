source.dir <- "/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
source(paste0(source.dir,"Functions/Compute_Number_CNAs_per_Cell.R"))

input_path <- c("Baron2016/output/scaledata/genewise/FixedMovingWindow/","Segerstolpe2016/output/scaledata/genewise/FixedMovingWindow/","Xin2016/output/scaledata/genewise/FixedMovingWindow/",
                "Muraro2016/output/scaledata/genewise/FixedMovingWindow/","Lawlor2017/output/scaledata/genewise/FixedMovingWindow/","Enge2017/output/scaledata/genewise/FixedMovingWindow/",
                "Peng2019/output/scaledata/v2/genewise/healthy/FixedMovingWindow/")
input_path <- c("Baron2016/output/scaledata/genewise/FixedMovingWindow/","Muraro2016/output/scaledata/genewise/FixedMovingWindow/",
                "Peng2019/output/scaledata/v2/genewise/healthy/FixedMovingWindow/") # Downsampling/Prob/
output_dir <- c("FDR01/") #c("chr_filtered/Rm_Segments/threshold_10/FDR005/")
FDR <- 0.1
# size_threshold <- 10

# 0.05 0.1 0.2 0.3 

for(i in 1:length(input_path)){
  name <-  strsplit(input_path[i],"/")[[1]][1]
  counts_se <- readRDS(file=paste0("/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Pancreas/",input_path[i],'data.rds'))
  Segments <- readRDS(paste0("/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Pancreas/",input_path[i],"CBS_pvalues.rds"))
  print(paste0("Chr removed in dataset: ",name))
  print(setdiff(y=names(table(Segments$chrom)), x=as.character(1:23)))
  
  Segments$CNA<-NA
  Segments[which(Segments$seg.mean > 0),"CNA"] <- "Amp"
  Segments[which(Segments$seg.mean < 0),"CNA"] <- "Del"
  
  #rm non sig calls
  Segments[which(Segments$FDR > FDR),"CNA"] <- "Zero"
  # # rm small CNAs
  # Segments[which(Segments$num.mark < size_threshold),"CNA"] <- "Zero" ####################
  Segments_split <- split(Segments, Segments$Cell)
  
  Cells <- unique(Segments$Cell)
  
  Number_CNAs_per_cell_Sep <- Compute_Number_CNAs_per_Cell(Segments_split,counts_se@meta.data,Cells=Cells)
  Number_CNAs_per_cell_Comb <- Compute_Number_CNAs_per_Cell(Segments_split,counts_se@meta.data,Cells=Cells,Sep=FALSE)
  
  if(name=="Baron2016"){
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    Number_CNAs_per_cell_Comb$Age <- revalue(Number_CNAs_per_cell_Comb$Age, c('Donor1'=17,'Donor2'=51,'Donor3'=38,'Donor4'=59))
    Number_CNAs_per_cell_Comb$Age <- as.numeric(as.character(Number_CNAs_per_cell_Comb$Age))
    
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    Number_CNAs_per_cell_Comb$Sex <- revalue(Number_CNAs_per_cell_Comb$Sex,  c('Donor1'='m','Donor2'='f','Donor3'='m','Donor4'='f'))
    Number_CNA_metadata <-Number_CNAs_per_cell_Comb
  
    #Sep
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    Number_CNAs_per_cell_Sep$Age <- revalue(Number_CNAs_per_cell_Sep$Age, c('Donor1'=17,'Donor2'=51,'Donor3'=38,'Donor4'=59))
    Number_CNAs_per_cell_Sep$Age <- as.numeric(as.character(Number_CNAs_per_cell_Sep$Age))
    #Number_CNAs_per_cell_Sep$CNA <- revalue(Number_CNAs_per_cell_Sep$CNA, c("Amp"="Amplification","Del"="Deletion"))
    
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    Number_CNAs_per_cell_Sep$Sex <- revalue(Number_CNAs_per_cell_Sep$Sex,  c('Donor1'='m','Donor2'='f','Donor3'='m','Donor4'='f'))
  }
  if(name=="Segerstolpe2016"){
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    Number_CNAs_per_cell_Comb$Age <- revalue(Number_CNAs_per_cell_Comb$Age, c('Donor_1'=22,'Donor_2'=23,'Donor_3'=25,'Donor_4'=27,
                                                                                'Donor_5_T2D'=37,'Donor_6'=43,'Donor_7'=48,'Donor_8_T2D'=52,
                                                                                'Donor_9_T2D'=55,'Donor_10_T2D'=57))
    Number_CNAs_per_cell_Comb$Age <- as.numeric(as.character(Number_CNAs_per_cell_Comb$Age))
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    Number_CNAs_per_cell_Comb$Sex <- revalue(Number_CNAs_per_cell_Comb$Sex, c('Donor_1'='m','Donor_2'='m','Donor_3'='m','Donor_4'='m','Donor_5_T2D'='f','Donor_6'='m','Donor_7'='f','Donor_8_T2D'='m','Donor_9_T2D'='f','Donor_10_T2D'='m'))
    Number_CNA_metadata <-Number_CNAs_per_cell_Comb
    
    
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    Number_CNAs_per_cell_Sep$Age <- revalue(Number_CNAs_per_cell_Sep$Age, c('Donor_1'=22,'Donor_2'=23,'Donor_3'=25,'Donor_4'=27,
                                                                                'Donor_5_T2D'=37,'Donor_6'=43,'Donor_7'=48,'Donor_8_T2D'=52,
                                                                                'Donor_9_T2D'=55,'Donor_10_T2D'=57))
    Number_CNAs_per_cell_Sep$Age <- as.numeric(as.character(Number_CNAs_per_cell_Sep$Age))
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    Number_CNAs_per_cell_Sep$Sex <- revalue(Number_CNAs_per_cell_Sep$Sex, c('Donor_1'='m','Donor_2'='m','Donor_3'='m','Donor_4'='m','Donor_5_T2D'='f','Donor_6'='m','Donor_7'='f','Donor_8_T2D'='m','Donor_9_T2D'='f','Donor_10_T2D'='m'))

  }
  if(name=="Xin2016"){
    Donor_info <- readODS::read_ods('/data/public/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/Pancreas/Xin/data/Donor_info.ods')
    Donor_info$Donor <- c('Donor02','Donor07','Donor01','Donor14','Donor04','Donor18','Donor05','Donor17','Donor03','Donor11','Donor06','Donor15','Donor16','Donor08','Donor13','Donor09','Donor10','Donor12')
    Donor_info <- Donor_info[order(Donor_info$Donor),c('Donor','Donor ID','Age','Ethnicity','Gender','BMI','HbA1c','Weight (lbs)','Cause of Death')]
    
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    all.equal(levels(Number_CNAs_per_cell_Comb$Age),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Comb$Age) <- Donor_info$Age
    Number_CNAs_per_cell_Comb$Age <- as.numeric(as.character(Number_CNAs_per_cell_Comb$Age))
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    all.equal(levels(Number_CNAs_per_cell_Comb$Sex),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Comb$Sex) <- Donor_info$Gender
    Number_CNA_metadata <-Number_CNAs_per_cell_Comb
    
    
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    all.equal(levels(Number_CNAs_per_cell_Sep$Age),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Sep$Age) <- Donor_info$Age
    Number_CNAs_per_cell_Sep$Age <- as.numeric(as.character(Number_CNAs_per_cell_Sep$Age))
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    all.equal(levels(Number_CNAs_per_cell_Sep$Sex),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Sep$Sex) <- Donor_info$Gender
  }
  if(name=="Muraro2016"){
    
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    Number_CNAs_per_cell_Comb$Age <- revalue(Number_CNAs_per_cell_Comb$Age, c('Donor_1'=23,'Donor_2'=48,'Donor_3'=54,'Donor_4'=59))
    Number_CNAs_per_cell_Comb$Age <- as.numeric(as.character(Number_CNAs_per_cell_Comb$Age))
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    Number_CNAs_per_cell_Comb$Sex <- revalue(Number_CNAs_per_cell_Comb$Sex, c('Donor_1'='m','Donor_2'='f','Donor_3'='m','Donor_4'='m'))
    Number_CNA_metadata <-Number_CNAs_per_cell_Comb
    
    
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    Number_CNAs_per_cell_Sep$Age <- revalue(Number_CNAs_per_cell_Sep$Age, c('Donor_1'=23,'Donor_2'=48,'Donor_3'=54,'Donor_4'=59))
    Number_CNAs_per_cell_Sep$Age <- as.numeric(as.character(Number_CNAs_per_cell_Sep$Age))
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    Number_CNAs_per_cell_Sep$Sex <- revalue(Number_CNAs_per_cell_Sep$Sex, c('Donor_1'='m','Donor_2'='f','Donor_3'='m','Donor_4'='m'))
  }
  if(name=="Enge2017"){
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    Number_CNAs_per_cell_Comb$Age <- revalue(Number_CNAs_per_cell_Comb$Age, c('Donor_1'=1,'Donor_2'=5,'Donor_3'=6,'Donor_4'=21,
                                                                                'Donor_5'=22,'Donor_6'=38,'Donor_7'=44,'Donor_8'=54))
    Number_CNAs_per_cell_Comb$Age <- as.numeric(as.character(Number_CNAs_per_cell_Comb$Age))
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    Number_CNAs_per_cell_Comb$Sex <- revalue(Number_CNAs_per_cell_Comb$Sex, c('Donor_1'='m','Donor_2'='m','Donor_3'='m','Donor_4'='m','Donor_5'='m','Donor_6'='f','Donor_7'='f','Donor_8'='m'))
    Number_CNA_metadata <-Number_CNAs_per_cell_Comb
    
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    Number_CNAs_per_cell_Sep$Age <- revalue(Number_CNAs_per_cell_Sep$Age, c('Donor_1'=1,'Donor_2'=5,'Donor_3'=6,'Donor_4'=21,
                                                                              'Donor_5'=22,'Donor_6'=38,'Donor_7'=44,'Donor_8'=54))
    Number_CNAs_per_cell_Sep$Age <- as.numeric(as.character(Number_CNAs_per_cell_Sep$Age))
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    Number_CNAs_per_cell_Sep$Sex <- revalue(Number_CNAs_per_cell_Sep$Sex, c('Donor_1'='m','Donor_2'='m','Donor_3'='m','Donor_4'='m','Donor_5'='m','Donor_6'='f','Donor_7'='f','Donor_8'='m'))

  }
  if(name=="Peng2019"){
    
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    Number_CNAs_per_cell_Comb$Age <- revalue(Number_CNAs_per_cell_Comb$Age, c("N1"="30", "N2"="31","N3"="34", "N4"="41","N5"="42", "N6"="50","N7"="52", "N8"="53","N9"="55", "N10"="64","N11"="65", 
                                                                              "T1"="36","T2"="44", "T3"="51","T4"="52","T5"="54", "T6"="54","T7"="54","T8"="56", "T9"="58","T10"="58","T11"="59", "T12"="59",
                                                                              "T13"="59","T14"="61", "T15"="64","T16"="64","T17"="65", "T18"="66","T19"="67","T20"="67", "T21"="68","T22"="70","T23"="71", 
                                                                              "T24"="72"))
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    Number_CNAs_per_cell_Comb$Sex <- revalue(Number_CNAs_per_cell_Comb$Sex, c("N1"="f", "N2"="f","N3"="m", "N4"="m","N5"="f", "N6"="m","N7"="f", "N8"="m","N9"="m", "N10"="f","N11"="f", 
                                                                              "T1"="m","T2"="f", "T3"="m","T4"="m","T5"="m", "T6"="f","T7"="m","T8"="f", "T9"="f","T10"="f","T11"="f", "T12"="m",
                                                                              "T13"="m","T14"="m", "T15"="m","T16"="m","T17"="f", "T18"="f","T19"="f","T20"="f", "T21"="f","T22"="m","T23"="f", 
                                                                              "T24"="f"))
    Number_CNA_metadata<-Number_CNAs_per_cell_Comb
    
    #SEP
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    Number_CNAs_per_cell_Sep$Age <- revalue(Number_CNAs_per_cell_Sep$Age, c("N1"="30", "N2"="31","N3"="34", "N4"="41","N5"="42", "N6"="50","N7"="52", "N8"="53","N9"="55", "N10"="64","N11"="65", 
                                                                              "T1"="36","T2"="44", "T3"="51","T4"="52","T5"="54", "T6"="54","T7"="54","T8"="56", "T9"="58","T10"="58","T11"="59", "T12"="59",
                                                                              "T13"="59","T14"="61", "T15"="64","T16"="64","T17"="65", "T18"="66","T19"="67","T20"="67", "T21"="68","T22"="70","T23"="71", 
                                                                              "T24"="72"))
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    Number_CNAs_per_cell_Sep$Sex <- revalue(Number_CNAs_per_cell_Sep$Sex, c("N1"="f", "N2"="f","N3"="m", "N4"="m","N5"="f", "N6"="m","N7"="f", "N8"="m","N9"="m", "N10"="f","N11"="f", 
                                                                              "T1"="m","T2"="f", "T3"="m","T4"="m","T5"="m", "T6"="f","T7"="m","T8"="f", "T9"="f","T10"="f","T11"="f", "T12"="m",
                                                                              "T13"="m","T14"="m", "T15"="m","T16"="m","T17"="f", "T18"="f","T19"="f","T20"="f", "T21"="f","T22"="m","T23"="f", 
                                                                              "T24"="f"))
  }
  if(name=="Lawlor2017"){
    metadata_Lawlor <- readRDS('/data/public/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/Pancreas/Lawlor/data/metadata_Lawlor.rds')
    
    Donor_info <- rbind(subset(metadata_Lawlor, Donor == 'Donor01')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor02')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor03')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor04')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor05')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor06')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor07')[1,c('Donor','Age','Sex','Condition','Race','BMI')],
                        subset(metadata_Lawlor, Donor == 'Donor08')[1,c('Donor','Age','Sex','Condition','Race','BMI')])
    
    Number_CNAs_per_cell_Comb$Age <- Number_CNAs_per_cell_Comb$Donor
    levels(Number_CNAs_per_cell_Comb$Age)
    all.equal(levels(Number_CNAs_per_cell_Comb$Age),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Comb$Age) <- Donor_info$Age
    Number_CNAs_per_cell_Comb$Age <- as.numeric(as.character(Number_CNAs_per_cell_Comb$Age))
    Number_CNAs_per_cell_Comb$Sex <- Number_CNAs_per_cell_Comb$Donor
    all.equal(levels(Number_CNAs_per_cell_Comb$Sex),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Comb$Sex) <- Donor_info$Sex
    Number_CNA_metadata <-Number_CNAs_per_cell_Comb
    
    Number_CNAs_per_cell_Sep$Age <- Number_CNAs_per_cell_Sep$Donor
    levels(Number_CNAs_per_cell_Sep$Age)
    all.equal(levels(Number_CNAs_per_cell_Sep$Age),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Sep$Age) <- Donor_info$Age
    Number_CNAs_per_cell_Sep$Age <- as.numeric(as.character(Number_CNAs_per_cell_Sep$Age))
    Number_CNAs_per_cell_Sep$Sex <- Number_CNAs_per_cell_Sep$Donor
    all.equal(levels(Number_CNAs_per_cell_Sep$Sex),as.character(Donor_info$Donor))
    levels(Number_CNAs_per_cell_Sep$Sex) <- Donor_info$Sex
  }
  
  
  saveRDS(Number_CNA_metadata,file=paste0('output/',output_dir,name,'_Number_CNA_metadata.rds')) # Downsampling/Prob/
  saveRDS(Number_CNAs_per_cell_Sep,file=paste0('output/',output_dir,name,'_Number_CNA_metadata_Amp_Del.rds')) # Downsampling/Prob/
  
}


# [1] "Chr removed in dataset: Baron2016"
# [1] "8"  "13" "15" "18" "20" "21" "22" "23"
# [1] "Chr removed in dataset: Segerstolpe2016"
# [1] "18" "21" "23"
# [1] "Chr removed in dataset: Xin2016"
# [1] "13" "18" "20" "21" "22" "23"
# [1] "Chr removed in dataset: Muraro2016"
# [1] "13" "18" "21" "23"
# [1] "Chr removed in dataset: Lawlor2017"
# [1] "13" "21" "23"
# [1] "Chr removed in dataset: Enge2017"
# [1] "13" "18" "21" "23"
# [1] "Chr removed in dataset: Peng2019"
# [1] "13" "18" "21" "23"