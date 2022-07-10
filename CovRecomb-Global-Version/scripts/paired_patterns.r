# Libraries
library(optparse)
library(ggplot2)

option_list <- list(
  make_option(c("-d", "--file_address"), type = "character", default = "/home/soniali/Desktop/CovRecomb-Global-Version/data/2022_02_12/3_recom_pattern/",
              action = "store", help = "The file address of the lineage_paired.csv and vairant_paired.csv"
  ),
  make_option(c("-l", "--lineage_paired_filename"), type = "character", default = "lineage_paired.csv",
              action = "store", help = "The filename of lineage_paired.csv"
  ),
  make_option(c("-v", "--variant_paired_filename"), type = "character", default = "variant_paired.csv",
              action = "store", help = "The filename of vairant_paired.csv"
  )
)

opt = parse_args(OptionParser(option_list = option_list, usage = "The script to plot the heatmap for lineage-paired and variant-paired patterns for independent recombination events."))
out_dir_pdf  = opt$file_address
# Read the file lineage-paired
lineage_paired_file = paste0(opt$file_address,opt$lineage_paired_filename)

if (file.exists(lineage_paired_file)){
  df <- read.csv(lineage_paired_file,
                 header = FALSE)
  colnames(df) <- c("A","B","Recom")
  size_group = c()
  for (i in df$Recom){
    i = as.integer(i)
    if (i >= 8){
      size_group = c(size_group,"9")
    } else if (i >= 4){
      size_group = c(size_group,"4 ~ 5")
    } else if (i >= 1) {
      size_group = c(size_group,"1 ~ 3")
    }
  }
  df$size_group = size_group
  
  # PDF
  title_main = "Lineage-paired recombination events"
  p1 <- ggplot(df,aes(x=A, y=B, size = Recom, color = size_group)) +#
    scale_fill_manual(values=c("1~3" = "darkgreen", "4~5" = "blue", "9" = "red"))+
    geom_point(alpha=0.5) +#,color = "#8696a7"
    scale_size(range = c(3,35), name="Number of recombination event",breaks = c(1,4,7)) +
    # scale_color_viridis(discrete=TRUE, guide="none") +
    theme_light()+ #+
    theme(legend.position="none")
  
  p2 <- p1 + labs(x = "Lineage X", y = "Lineage Y", title = title_main)
  p3 <- p2 + theme(axis.text.x=element_text(size = 8, vjust = 0.5, hjust = 0.5, angle = 45, face = "plain", color='black'),
                   axis.text.y=element_text(size = 10, vjust = 0.5, face = "plain", hjust = 0.5, angle = 0, color='black'))
  p4 <- p3 + theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold"))  #将title居中
  p4 <- p4 + theme(axis.title.y = element_text(size = 14, face = "bold"))
  p5 <- p4 + theme(axis.title.x = element_text(size = 14, face = "bold"))
  
  ggsave(filename = paste(paste(out_dir_pdf,title_main),".pdf"),plot = p5,device = "pdf",width = 16, height =16)
} else {
  print("There is no lineage_paired file, please cheack the raw folder.")
}


# Read the file variant-paired
variant_paired_file = paste0(opt$file_address, opt$variant_paired_filename)
if (file.exists(variant_paired_file)){
  df <- read.csv(variant_paired_file,
                 header = TRUE)
  colnames(df) <- c("A","B","Recom")
  size_group = c()
  for (i in df$Recom){
    i = as.integer(i)
    if (i >= 17){
      size_group = c(size_group,"17-22")
    } else if (i >= 10){
      size_group = c(size_group,"10-16")
    } else if (i >= 4) {
      size_group = c(size_group,"4-9")
    }else if (i >= 1) {
      size_group = c(size_group,"2-3")
    }
  } 
  df$size_group = size_group
  
  # PDF
  title_main = "Variant-paired recombination events"
  p1 <- ggplot(df,aes(x=X, y=Y, size = Recom, color = "c0c0c0")) +
    geom_point(alpha=0.5) +
    scale_size(range = c(3,10), name="Number of recombination event",breaks = c(1,4,7)) +
    theme_light()+ 
    theme(legend.position="bottom")
  
  p2 <- p1 + labs(x = "Variant X", y = "Variant Y", title = "")
  p3 <- p2 + theme(axis.text.x=element_text(size = 6, vjust = 0.5, hjust = 0.5, angle = 45, face = "bold", color='black'),
                   axis.text.y=element_text(size = 6, vjust = 0.5, face = "bold", hjust = 0.5, angle = 0, color='black'))
  p4 <- p3 + theme(plot.title = element_text(hjust = 0.5,size=12,face = "bold"))  #将title居中
  p4 <- p4 + theme(axis.title.y = element_text(size = 10, face = "bold"))
  p5 <- p4 + theme(axis.title.x = element_text(size = 10, face = "bold"))
  
  ggsave(filename = paste(paste(out_dir_pdf,title_main),".pdf"),plot = p5,device = "pdf",width = 4, height =4.5)
} else {
  print("There is no variant_paired file, please cheack the raw folder.")
}

