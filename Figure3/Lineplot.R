library(dplyr)
library(ggpubr)

## Function to sort data frames according to the genomic position of the genes
## Input:   data = data matrix; rows -> genes_hgnc
##          gencode = gencode object
genome_sort <- function(data){
  data <- data[order(data$chr,data$start),]
  chromosomeorder <- data.frame(chr=data$chr,
                                start=data$start,
                                end=data$end,
                                ensembl=data$gene_id,
                                hgnc=data$gene_name)
  assign("chromosomeorder", chromosomeorder, .GlobalEnv)
  return(data)
}

lineplot_CNAs <- function(Del_Amp_profile, file_name, ylimit = NULL , chromosome = NULL, per_celltype = FALSE, pdf=FALSE,lab_title=NULL,path="",centromer_index=NULL){
  if(per_celltype == TRUE){
    message("Plot lineplots per celltype")

    Del_Amp_profile_plot <- lapply(Del_Amp_profile, genome_sort)
    Del_Amp_profile_plot <- ldply (Del_Amp_profile, data.frame)
    number_celltypes <- length(Del_Amp_profile)
    number_genes <- nrow(Del_Amp_profile[[1]])
    height_plot <- 1000
    print(names(Del_Amp_profile))
  }else{
    message("Plot combined lineplots")
    Del_Amp_profile <- genome_sort(Del_Amp_profile)
    Del_Amp_profile_plot <- Del_Amp_profile
    number_celltypes <- 1
    number_genes <- nrow(Del_Amp_profile)
    height_plot <-500
  }

  if ( is.null(chromosome) ){
    message('Plot all chromosomes.')
    chromosome_boundaries <- which(duplicated(chromosomeorder$chr)==FALSE)[-1]
    if(per_celltype == TRUE){
      p <- plot_lineplot(Del_Amp_profile_plot, ylimit, start=1, end=number_genes, number_celltypes, type='Del',
                         chromosome_boundaries, ylab_title="Percent of Deletions",title = NULL, per_celltype)
      p<-ggarrange(plotlist = p,
                ncol = 1, nrow = length(p)+1,
                heights = c(rep(1, length(p)),0.15))
      if(pdf){
        message('Save deletion lineplots in pdf file:')
        pdf(paste0("plots/Deletion_Profile_", file_name,".pdf"), width=15, height = height_plot/100+5)
        print(annotate_figure(p,
                              left = text_grob("Percent of Deletions",  rot = 90,size = 30),
                              bottom = text_grob("Chromosome",size = 30)#,right = text_grob(Anot,  rot = -90,size = 20)
        ))
        dev.off()
      }

      message('Save deletion lineplots in png file:')
      png(paste0("plots/Deletion_Profile_", file_name,".png"), width=1000, height = height_plot)
      print(annotate_figure(p,
                            left = text_grob("Percent of Deletions",  rot = 90,size = 30),
                            bottom = text_grob("Chromosome",size = 30)#,right = text_grob(Anot,  rot = -90,size = 20)
                            ))
      dev.off()

      p <- plot_lineplot(Del_Amp_profile_plot, ylimit, start=1, end=number_genes, number_celltypes, type='Amp',
                         chromosome_boundaries, ylab_title="Percent of Amplifications",title = NULL, per_celltype)
      p<-ggarrange(plotlist = p,
                   ncol = 1, nrow = length(p)+1,
                   heights = c(rep(1, length(p)),0.15))
      if(pdf){
        message('Save amplification lineplots in pdf file.')
        pdf(paste0("plots/Amplification_Profile_", file_name,".pdf"), width=15, height = height_plot/100+5)
        print(annotate_figure(p,
                              left = text_grob("Percent of Amplifications",  rot = 90,size = 30),
                              bottom = text_grob("Chromosome",size = 30)#,right = text_grob(Anot,  rot = -90,size = 25)
        ))
        dev.off()
      }
      message('Save amplification lineplots in png file.')
      png(paste0("plots/",path,"Amplification_Profile_", file_name,".png"), width=1000, height = height_plot)
      print(annotate_figure(p,
                            left = text_grob("Percent of Amplifications",  rot = 90,size = 30),
                            bottom = text_grob("Chromosome",size = 30)#,right = text_grob(Anot,  rot = -90,size = 25)
      ))
      dev.off()

    }else{
        if( is.null(lab_title) ){lab_title="Percent of "}

      p <- plot_lineplot(Del_Amp_profile_plot, ylimit, start=1, end=number_genes, number_celltypes, type='Del',
                         chromosome_boundaries, ylab_title=paste0(lab_title,"Deletions"),title = NULL, per_celltype,centromer_index)
      message('Save deletion lineplots in png file.')
      png(paste0("plots/",path,"Deletion_Profile_", file_name,".png"), width=1000, height = height_plot)
      print(p)
      dev.off()
      
      if(pdf){
        message('Save deletion lineplots in pdf file:')
        pdf(paste0("plots/Deletion_Profile_", file_name,".pdf"), width=15, height = height_plot/100-2)
        print(p)
        dev.off()
      }
      
      p <- plot_lineplot(Del_Amp_profile_plot, ylimit, start=1, end=number_genes, number_celltypes, type='Amp',
                         chromosome_boundaries, ylab_title=paste0(lab_title,"Amplifications"),title = NULL, per_celltype,centromer_index)
      message('Save amplification lineplots in png file.')
      png(paste0("plots/",path,"Amplification_Profile_", file_name,".png"), width=1000, height = height_plot)
      print(p)
      dev.off()
      
      if(pdf){
        message('Save amplification lineplots in pdf file.')
        pdf(paste0("plots/Amplification_Profile_", file_name,".pdf"), width=15, height = height_plot/100-2)
        print(p)
        dev.off()
      }
      
    }
  }else{
    message(paste('Plot chromosome',chromosome,'.'))

    if(per_celltype == TRUE){
      Del_Amp_profile_chr <- lapply(Del_Amp_profile, function(x) x[chromosomeorder$chr == chromosome,])
      Del_Amp_profile_chr <- ldply(Del_Amp_profile_chr, data.frame)

    }else{Del_Amp_profile_chr <- Del_Amp_profile[chromosomeorder$chr == chromosome,]}

    chr_index <- which(chromosomeorder$chr == chromosome)
    start <- chr_index[1]
    end <- chr_index[length(chr_index)]


    if(per_celltype == TRUE){
############################
    }else{
      message('Save deletion lineplots in png file.')
      p <- plot_lineplot(Del_Amp_profile_chr, ylimit, start=start, end=end, number_celltypes, type='Del',
                         ylab_title="Percent of Deletions",title = paste("Chromosome",chromosome),
                         per_celltype = per_celltype)
      png(paste0("plots/",path,"Deletion_Profile_",file_name,"_chr_",chromosome,".png"),width=1000, height = height_plot)
      print(p)
      dev.off()
      
      if(pdf){
        message('Save deletion lineplots in pdf file:')
        pdf(paste0("plots/Deletion_Profile_", file_name,".pdf"), width=15, height = height_plot/100+5)
        print(annotate_figure(p,
                              left = text_grob("Percent of Deletions",  rot = 90,size = 30),
                              bottom = text_grob("Chromosome",size = 30)#,right = text_grob(Anot,  rot = -90,size = 20)
        ))
        dev.off()
      }
      
      message('Save amplification lineplots in png file.')
      p <- plot_lineplot(Del_Amp_profile_chr, ylimit, start=start, end=end, number_celltypes, type='Amp',
                         ylab_title="Percent of Amplifications",title = paste("Chromosome",chromosome),
                         per_celltype = per_celltype)
      png(paste0("plots/",path,"Amplification_Profile_",file_name,"_chr_",chromosome,".png"),width=1000, height = height_plot)
      print(p)
      dev.off()
      
      if(pdf){
        message('Save amplification lineplots in pdf file.')
        pdf(paste0("plots/Amplification_Profile_", file_name,".pdf"), width=15, height = height_plot/100+5)
        print(annotate_figure(p,
                              left = text_grob("Percent of Amplifications",  rot = 90,size = 30),
                              bottom = text_grob("Chromosome",size = 30)#,right = text_grob(Anot,  rot = -90,size = 25)
        ))
        dev.off()
      }
    }

  }
}


plot_lineplot <- function(Del_Amp_profile_plot,ylimit,start,end,number_celltypes,type,chromosome_boundaries=NULL,
                          ylab_title,title = NULL, per_celltype,centromer_index=NULL){
  if ( is.null(ylimit) ){
    ylimit <- max(Del_Amp_profile_plot[,type], na.rm = TRUE)
  }
  #browser()
  x_label<-c(1,chromosome_boundaries)+c(diff(c(1,chromosome_boundaries))/2,600)
  y_label<-rep(c(-ylimit/10,-(ylimit/10)-(ylimit/20)),round(length(c(1,chromosome_boundaries))/2))
  y_label<-y_label[1:length(x_label)]
  if(type=="Amp"){color<-"#D16103"}else{color<-"#005AB5"}
  #print(chromosome_boundaries)
  #print(colnames(Del_Amp_profile_plot))
  p <- ggplot(Del_Amp_profile_plot,  aes_string(x = rep(start:end,number_celltypes), y =  type)) +
    geom_line(size = 1,color=color)+
    xlim(start,end+500) +
    geom_vline(xintercept = chromosome_boundaries,col="black") +
    geom_vline(xintercept = centromer_index,col="darkgrey",linetype = "dashed") +
    xlab("\n\nChromosomes") +
    ylab(ylab_title) +
    ggtitle(title) +
    theme(axis.title = element_text(size = 25), axis.text.y = element_text(size = 20),axis.text.x = element_blank(),
          title= element_text(size = 25), legend.title = element_text(size = 30), axis.ticks.x = element_blank(),
          legend.text = element_text(size = 26)) +
    scale_color_viridis_d()+
    annotate("text", x =x_label, y = y_label,
             label = unique(chromosomeorder$chr), size = 8)+
    coord_cartesian(ylim = c(0, ylimit), clip = "off")


   if(per_celltype == TRUE){
     Tissue <- unique(Del_Amp_profile_plot$tissue)
     p<-list()
     for(tissue in 1:(length(Tissue)-1) ){
   
     df<-Del_Amp_profile_plot[which(Del_Amp_profile_plot$tissue == Tissue[tissue]),]
     p[[tissue]] <- ggplot(df,  aes_string(x = rep(start:end,1), y =  type)) +
       geom_line(size = 1,color=color)+
       geom_vline(xintercept = chromosome_boundaries,col="darkgrey") +
       xlab("") +
       ylab("") +
       #ggtitle(Tissue[tissue]) +
       theme(axis.title = element_blank(), axis.text.y = element_text(size = 18),axis.text.x = element_blank(),
             title= element_text(size = 25), legend.title = element_text(size = 30), axis.ticks.x = element_blank(),
             legend.text = element_text(size = 26)) +
       scale_color_viridis_d()+
       annotate("text", x =end+1500, y = ylimit-ylimit/2,
                label = Tissue[tissue], size = 7, angle = -90)+
       coord_cartesian(ylim = c(0, ylimit),xlim = c(start,end+500), clip = "off")
     }
     #browser()
    # last tissue
     x_label<-c(1,chromosome_boundaries)+c(diff(c(1,chromosome_boundaries))/2,600)
     y_label<-rep(c(-ylimit/5),length(unique(chromosomeorder$chr)))

      tissue <-length(Tissue)
       df<-Del_Amp_profile_plot[which(Del_Amp_profile_plot$tissue == Tissue[tissue]),]
       p[[tissue]] <- ggplot(df,  aes_string(x = rep(start:end,1), y =  type)) +
         geom_line(size = 1,color=color)+
         geom_vline(xintercept = chromosome_boundaries,col="darkgrey") +
         xlab("") +
         ylab("") +
         #ggtitle(Tissue[tissue]) +
         theme(axis.title = element_blank(), axis.text.y = element_text(size = 18),axis.text.x = element_blank(),
               title= element_text(size = 25), legend.title = element_text(size = 30), axis.ticks.x = element_blank(),
               legend.text = element_text(size = 26)) +
         scale_color_viridis_d()+
         annotate("text", x =x_label, y = y_label,
                  label = unique(chromosomeorder$chr), size = 8)+
         annotate("text", x =end+1500, y = ylimit-ylimit/2,
                  label = Tissue[tissue], size = 7, angle = -90)+
         coord_cartesian(ylim = c(0, ylimit), xlim = c(start,end+500),clip = "off") #

   }
  return(p)
}
