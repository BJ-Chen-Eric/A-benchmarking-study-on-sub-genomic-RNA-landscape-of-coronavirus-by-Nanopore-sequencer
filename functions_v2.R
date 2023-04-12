library(seqinr)
library(stringr)
library(Biostrings)
library(bioseq)
library(dplyr)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(msa)
library(stringdist)
library(cowplot)
library(gggenes)
library(ggtranscript)
library(ComplexHeatmap)
library(circlize)
library(ComplexUpset)
library(paletteer)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(data.table) 
library(ComplexHeatmap)
library(circlize)
library(ComplexUpset)
library(ggtranscript)
library(gggenes)
library(paletteer)
## functions
filter_raw <- function(list, short_frag=20, score=60)  {
  for(i in seq_along(list))  {
    A <- list[[i]][,1:6]
    A <- A[A$V5 >=score, ]
    A$frag <- A$V3-A$V2+1
    read <- A$V4[A$frag < short_frag]
    list[[i]] <- A[!(A$V4 %in% read), 1:6]
  }
  return(list)
}


align_specy_table <- function(path, replace_pattern2='', ref_name, prefix='')  {
  read_as_list <- function(path, seq='', prefix='')  {
    ipak <- function(pkg){
      new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
      if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
      sapply(pkg, require, character.only = TRUE)
    }
    packages <- c('dplyr', 'purrr', 'fs')
    ipak(packages)
    mini_pre.bed <- list()
    bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
    bed_path <- as_fs_path(bed_path)
    mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
    names(mini_pre.bed) <- bed_path
    pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
    for(i in 1:2)
    {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
    return(mini_pre.bed)
    sapply(packages, require, character.only = TRUE)
    lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
  }
  A <- read_as_list(path = path, prefix = prefix)
  class <- as.data.frame(matrix(nrow = length(A), ncol = 4))
  for(i in seq_along(A)) {
    B <- as.data.frame(A[[i]])
    B <- as.data.frame(B[,3])
    class[i,'libs'] <- names(A)[i]
    unmap <- length(B[,1][B[,1] %in% '*'])
    class[i,1] <- unmap
    virus <- length(B[,1][B[,1] %in% ref_name])
    class[i,2] <- virus
    contamination <- length(B[,1][grep(B[,1], pattern = 'NC_001140.6')])
    class[i,3] <- contamination
    class[i,4] <- nrow(B)-unmap-virus-contamination
    # class[i,1:3] <- round(class[i,1:3]/rowSums(class[i,1:3])*100, digits = 2)
  }
  colnames(class)[1:4] <- c('Unmap', 'Virus', 'RCS', 'Host')
  return(class)
}



barplot_align_proportion <- function(df
                                     , x_title_size=36, x_test_size=28
                                     , y_title_size=36, y_text_size=28
                                     , title_size=36, legend_size=28
                                     , caption_size=24, label_size = 7
                                     , line_size=2, point_size=2
                                     , title='Library reads species proportion'
                                     , write=F, write_path=getwd(), out_prefix='')
{
  if(!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)
  gg_theme= theme(
    axis.title.x = element_text(size = x_title_size),
    axis.text.x = element_text(size =  x_test_size),
    axis.title.y = element_text(size = y_title_size),
    axis.text.y = element_text(size = y_text_size),
    title = element_text(size = caption_size),
    legend.text = element_text(size = legend_size), 
    legend.key.size = unit(2, 'lines'))
  A <- df %>% 
    dplyr::mutate(ID = 1:n()) %>%  
    tidyr::gather(key, value, c(1:3)) 
  
  df <- A
  df$libs <- sub(df$libs, pattern = '[0-9]', replacement = '')
  df[seq(1, nrow(df), 2),'ID'] <- seq(1, nrow(df)/2)
  df[seq(2, nrow(df), 2),'ID'] <- seq(1, nrow(df)/2)
  
  
  library(dplyr)
  df.summary <- df %>%
    dplyr::group_by(ID) %>%
    summarise(
      sd = sd(value, na.rm = TRUE),
      value = mean(value)
    )
  df.summary <- as.data.frame(df.summary)
  df.summary[,4:5] <- df[seq(1, nrow(df), 2), c(1,3)]
  
  df.summary$ID <- rep(c(1,2,3,4),3)
  df.summary$value <- as.numeric(df.summary$value)
  df.summary$sd <- as.numeric(df.summary$sd)
  df.summary$key <- factor(df.summary$key, levels = c('Host', 'Virus', 'Unmap'))
  
  # color <- colorRampPalette(c(bar_color2, bar_color1))
  library(ggplot2)
  plot <- 
    ggplot(df.summary, aes(key, value)) +
    geom_line(size=1, aes(group=ID, colour=libs))+
    geom_point(size=7, aes(shape=libs, colour=libs)) +
    geom_errorbar(aes(ymax=value+sd, ymin= value-sd, colour=libs), width = 0.15) +
    labs(title = title
         , x='species'
         , y='Percentage')+
    scale_color_discrete(name = "")+
    theme_linedraw()+gg_theme
  
  if(isTRUE(write)) {
    if(isFALSE(grepl(write_path, pattern = '&/')))
    {write_path <- paste(write_path,'/', sep = '')}
    png(paste(write_path, out_prefix, '_Alignment species proportion.png', sep = ''), width = 1600, height = 900)
    print(plot)
    dev.off()
  }
  if(isFALSE(write)) {
    print(plot)
  }
}



fragment_table <- function(list, max_frag=4, replicate=T)  {
  library(dplyr)
  frag <- as.data.frame(matrix('', nrow = 0,ncol = max_frag+2))
  colnames(frag) <- c(seq(1,max_frag), paste('≥', (max_frag+1), sep = ''), 'libs')
  index <- c(seq(1,max_frag), paste('≥', (max_frag+1), sep = ''))
  for(k in seq_along(list)) {
    BED <- list[[k]]
    BED[,2] <- BED[,2]+1
    BED <- BED[order(BED$V4), ]
    BED$inedx <- seq(1, nrow(BED))
    name <- as.data.frame(table(BED$V4))
    A <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    A$Var1 <- as.numeric(A$Var1)
    row <- nrow(A[A$Var1 <= max_frag,])
    A[(row+1), 'Freq'] <- colSums(as.matrix(A[(row+1):nrow(A), 'Freq']))
    A[(row+1),"Var1"] <- paste('≥', (max_frag+1),sep = '')
    A <- A[-c((row+2):nrow(A)),]
    for(j in seq(1,max_frag+1)){
      value <- A$Freq[A$Var1 %in% index[j]]
      frag[k,j] <- ifelse(identical(value, numeric(0)), yes = 0, no = value)
      frag[k,'libs'] <- names(list)[k]
    }
    frag[k,1:(max_frag+1)] <- round(as.numeric(frag[k,1:(max_frag+1)])/colSums(as.matrix(as.numeric(frag[k,1:(max_frag+1)])))*100, digits = 2)
  }
  A <- frag %>% 
    dplyr::mutate(ID = 1:n()) %>%  
    tidyr::gather(key, value, c(seq(1,(max_frag+1))))
  A$value <- as.numeric(A$value)
  if(isTRUE(replicate))  {
    df <- A
    df$libs <- sub(df$libs, pattern = '[0-9]', replacement = '')
    df[seq(1, nrow(df), 2),'ID'] <- seq(1, nrow(df)/2)
    df[seq(2, nrow(df), 2),'ID'] <- seq(1, nrow(df)/2)
    df.summary <- df %>%
      group_by(ID) %>%
      summarise(
        sd = sd(value, na.rm = TRUE),
        value = mean(value)
      )
    df.summary <- as.data.frame(df.summary)
    df.summary[,4:5] <- df[seq(1, nrow(df), 2), c(1,3)]
    
    df.summary$ID <- seq(1, 4)
    df.summary$value <- as.numeric(df.summary$value)
    df.summary$sd <- as.numeric(df.summary$sd)
    A <- df.summary
  }
  A$key <- factor(A$key, levels =  c(seq(1,max_frag), paste('≥', (max_frag+1), sep = '')))
  return(A)
}


read_as_list <- function(path, sep='auto', prefix='', header='auto', cols=NA, file_type='') {
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('dplyr', 'purrr', 'fs', 'seqinr', 'stringr', 'data.table')
  ipak(packages)
  out_list <- list()
  file_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  # file_type <- sub(str_extract(file_path[1], pattern = '\\..*'), pattern = '.', replacement = '')
  file_path <- as_fs_path(file_path)
  if(file_type == 'fasta')  {
    out_list <- file_path %>% map(.f = function(path){read.fasta(path, as.string = T, whole.header = T)})
  }
  # else if(file_type == 'csv')  {
  #   out_list <- file_path %>% map(.f = function(path){read.csv(path, header = header, sep = sep)})
  # }
  else{
    if(is.na(cols))  {
      out_list <- file_path %>% map(.f = function(path){fread(path, header = header
                                                              , stringsAsFactors = F, sep = sep)})
    }
    else  {
      out_list <- file_path %>% map(.f = function(path){fread(path, header = header, select = cols
                                                              , stringsAsFactors = F, sep = sep)})
    }
    out_list <- lapply(out_list, as.data.frame)
  }
  names(out_list) <- file_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)  {
    names(out_list) <- sub(names(out_list), pattern = pattern[i], replacement = '')
  }
  
  return(out_list)
  sapply(packages, require, character.only = TRUE)
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}


combine_rep <- function(list)  {
  x=0
  com_bed <- list()
  for(i in seq(1, length(list), 2))  {
    x=x+1
    rep1 <- rbind(list[[i]], list[[i+1]])
    com_bed[[x]] <- rep1
    names(com_bed)[x] <- sub(names(list)[i], pattern = '[0-9]', replacement = '')
  }
  return(com_bed)
}


# fragment_table_filter <- function(bed_list, gap)  {
#   frag_ori <- as.data.frame(matrix(data = '', nrow = 5, ncol = 8))
#   frag_filter <- as.data.frame(matrix(data = '', nrow = 5, ncol = 8))
#   colnames(frag_filter) <- colnames(frag_ori) <- names(bed)
#   for(i in seq(1,5))  {
#     A <- (bed, frag = i)
#     frag_ori[(i),] <- unlist(lapply(A, nrow))
#     if(i>=2)  {
#       for(j in seq_along(A))  {
#         B <- A[[j]]
#         dfgap <- as.data.frame(matrix(data = '', nrow = nrow(B), ncol = 1))
#         for(k in 2:i) {
#           dfgap[,(k-1)] <- B[,paste('start', k, sep = '')]-B[,paste('end', (k-1), sep = '')]
#         }
#         dfgap$min <- apply(dfgap, MARGIN = 1,FUN = min)
#         B$gap <- dfgap$min
#         B <- B[B$gap > gap, ]
#         A[[j]] <- B
#       }
#     }
#     frag_filter[(i),] <- unlist(lapply(A, nrow))
#   }
#   
#   frag_ori <- sapply(frag_ori, as.numeric)
#   # frag_ori[5,] <- 100-colSums(frag_ori[1:4,])+frag_ori[5,]
#   frag_filter <- sapply(frag_filter, as.numeric)
#   # frag_filter[5,] <- 100-colSums(frag_filter[1:4,])+frag_filter[5,]
#   for(i in seq_along(bed))  {
#     row <- length(unique(bed[[i]]$V4))
#     frag_ori[,i] <- frag_ori[,i]/row*100
#     frag_filter[,i] <- frag_filter[,i]/row*100
#   }
#   
#   diff <- (frag_ori-frag_filter)/frag_ori
#   
#   result <- list(frag_ori=frag_ori, frag_filter=frag_filter, diff=diff)
#   result <- lapply(result, as.data.frame)
#   for(i in seq_along(result))  {
#     A <- result[[i]]
#     rownames(A) <- c(seq(1,4), '≥5')
#     A['libs',] <- colnames(A)
#     result[[i]] <- A
#   }
#   # result <- lapply(result, t)
#   
#   return(result)
# }

plot_table <- function(recom_table, from, to)  {
  df <- recom_table %>% 
    dplyr::mutate(ID = 1:n()) %>%  
    tidyr::gather(key, value, c(seq(from ,to)))
  df$value <- as.numeric(df$value)
  df$libs <- sub(df[,1], pattern = '[0-9]', replacement = '')
  df[seq(1, nrow(df), 2),'ID'] <- seq(1, nrow(df)/2)
  df[seq(2, nrow(df), 2),'ID'] <- seq(1, nrow(df)/2)
  
  library(dplyr)
  df.summary <- df %>%
    dplyr::group_by(ID) %>%
    summarise(
      sd = sd(value, na.rm = TRUE),
      value = mean(value)) %>% as.data.frame()
  df.summary[,4:5] <- df[seq(1, nrow(df), 2), c(1,3)]
  df <- df.summary
  df$ID <- seq(1, 4)
  df$value <- as.numeric(df$value)
  df$sd <- as.numeric(df$sd)
  return(df)
}


align_length <- function(bed_list)  {
  out <- as.data.frame(matrix(data = '', nrow = 8, ncol = 6))
  colnames(out) <- c(1,2,3,4,'≥5','libs')
  for(j in 1:5) {
    library(dplyr)
    frag1 <- bed %>%
      specific_frag_format_organ(frag = j) 
    
    for(i in seq_along(frag1)) {
      A <- frag1[[i]]
      length <- as.data.frame(matrix(data = '', nrow = nrow(A), ncol = 1))
      for(k in 1:j) {
        length[,k] <- A[,paste('end', k, sep = '')]-A[,paste('start', k, sep = '')]
      }
      A$length <- rowSums(as.matrix(length))
      out[i,j] <- round(as.numeric(mean(A$length)), digits = 0)
      out[i,'libs'] <- names(frag1)[i]
    }
  }
  return(out)
}

alignment_qscore <- function(seq_path, bed_list)  {
  seq_sum <- read_as_list(path = seq_path, header = T)
  
  for(i in seq_along(seq_sum))  {
    A <- seq_sum[[i]]
    A <- A[,c(2,10,14,15,17,18)]
    seq_sum[[i]] <- A
  }
  
  out <- as.data.frame(matrix(data = '', nrow = 8, ncol = 5))
  colnames(out) <- c(1:4, '≥5')
  for(j in 1:5) {
    library(dplyr)
    frag1 <- bed_list %>%
      specific_frag_format_organ(frag = j) 
    
    qscore <- list()
    for(i in seq_along(frag1))  {
      A <- frag1[[i]]
      Q <- seq_sum[[i]]
      out[i,j] <- mean(Q$mean_qscore_template[Q$read_id %in% A$V4])
      out[i,'libs'] <- names(frag1)[i]
    }
  }
  return(out)
}

qctable <- function(qc_table_path)  {
  qc_table <- read.csv(file = qc_table_path)
  rownames(qc_table) <- qc_table[,1]; qc_table <- qc_table[,-c(1)]
  out <- qc_table[c(5,7,9,10),]
  out[,1:4] <- apply(out[,1:4], MARGIN = 2, as.numeric)
  out <- out*100
  out <- as.data.frame(t(out)) %>% tibble::rownames_to_column()
  colnames(out)[1] <- "libs"
  out <- out[,c(2:5,1)] 
  return(out)
}


all_frag_terminal <- function(list)  {
  output_list <- list()
  for(k in seq_along(list))  {
    BED <- list[[k]][,1:6]
    BED[,2] <- BED[,2]+1
    BED <- BED[order(BED$V4), ]
    name <- as.data.frame(table(BED$V4))
    split <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    fragment <- as.numeric(split[2:nrow(split),1])
    split <- BED[BED$V4 %in% name$Var1[name$Freq %in% 1],]
    split$frag <- 1
    split$s2 <- split$e2 <- ''
    split <- split[, c(1:3, 8, 9, 4:6, 7)]
    B <- split[order(split$V4),]
    colnames(B) <- c('chr', 'start1', 'end1', 'start2', 'end2', 'name', 'score', 'strand', 'frag')
    for(i in fragment) {
      split <- BED[BED$V4 %in% name$Var1[name$Freq %in% i],]
      split$frag <- i
      split <- split[order(split$V4),]
      D <- split[seq(i,nrow(split), i),2:7]
      C <- split[seq(1,nrow(split), i),1:3]
      C <- cbind(C,D)
      C[,9] <- i
      colnames(C) <- c('chr', 'start1', 'end1', 'start2', 'end2', 'name', 'score', 'strand', 'frag')
      B <- rbind(B,C)
    }
    # colnames(B) <- c('chr', 'start1', 'end1', 'start2', 'end2', 'name', 'score', 'strand', 'frag')
    output_list[[k]] <- B
    names(output_list)[k] <- names(list)[k]
  }
  return(output_list)
}


terminal_identification <- function(list, blurry, length, UTR5, table=T)  {
  out <- as.data.frame(matrix('', nrow = 8,ncol = 5))
  colnames(out) <- c("+5'UTR_+3'UTR","-5'UTR_+3'UTR","+5'UTR_-3'UTR","-5'UTR_-3'UTR",'libs')
  for(j in seq_along(list))  {
    A <- list[[j]]
    A$end3 <- A$end2 > (length-blurry)
    end3 <- A[A$end3 %in% T, -c(ncol(A))]
    if(nrow(end3)>0)  {
      end3$leader <- end3$start1 < UTR5
      end3$class <- ifelse(end3$leader==T, yes = "+5'UTR_+3'UTR", no = "-5'UTR_+3'UTR")
      end3 <- end3[,-c(8)]
    }
    noend3 <- A[A$end3 %in% F, -c(ncol(A))]
    if(nrow(noend3)>0) {
      noend3$leader <- noend3$start1 < UTR5
      noend3$class <- ifelse(noend3$leader==T, yes = "+5'UTR_-3'UTR", no = "-5'UTR_-3'UTR")
      noend3 <- noend3[,-c(8)]
    }
    B <- rbind(end3, noend3)
    list[[j]] <- B
    A <- as.data.frame(prop.table(table(B$class)), stringsAsFactors = F)
    for(i in seq(1:4)) {
      value <- A$Freq[A$Var1 %in% colnames(out)[i]]
      out[j,i] <- ifelse(identical(value, numeric(0)), yes = 0, no = value)
      out[j,'libs'] <- names(list)[j]
    }
  }
  if(isTRUE(table))  {
    library(dplyr)
    A <- out %>% 
      dplyr::mutate(ID = 1:n()) %>%  
      tidyr::gather(key, value, c(1:4)) 
    
    A$key <- factor(A$key, levels = c("+5'UTR_+3'UTR","-5'UTR_+3'UTR","+5'UTR_-3'UTR","-5'UTR_-3'UTR"))
    A$value <- as.numeric(A$value)
    return(A)
  }else(return(list))
  
}


# function 5: format organization depending on split frag
specific_frag_format_organ <- function(list, frag, min_gap=50)  {
  specific_fragment <- list()
  for(j in seq_along(list))  {
    BED <- list[[j]]
    BED[,2] <- BED[,2]+1
    BED$V4 <- paste(BED$V4, BED$V6, sep = '_')
    name <- as.data.frame(table(BED[,4]), stringsAsFactors =F)
    split <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    name <- name$Var1[name$Freq %in% frag]
    if(length(name) >0)  { 
      BED <- BED[BED[,4] %in% name,]
      BED <- BED[order(BED[,4], decreasing = F),]
      fragment <- list()
      if(frag==1) {
        B <- BED
      }
      else if(frag!=1) {
        for(i in 1:frag) {
          if(i ==1)
          {A <- BED[seq(i,nrow(BED), frag),][,1:3]}
          if(i ==frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:6]}
          if(i>1 & i<frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:3]}
          fragment[[i]] <- A
        }
        B <- as.data.frame(matrix(data = 0, nrow = nrow(fragment[[1]]), ncol = 0))
        for(i in rev(seq_along(fragment))) {
          A <- fragment[[i]]
          B <- cbind(A,B)
        }
      }
      
      colnames(B)[seq(2,(ncol(B)-3), 2)] <- 
        paste('start', seq(1,frag,1), sep = '')
      colnames(B)[seq(3,(ncol(B)-2), 2)] <- 
        paste('end', seq(1,frag,1), sep = '')
      specific_fragment[[j]] <- B
      names(specific_fragment)[j] <- names(list)[j]
    }
  }
  
  return(specific_fragment)
}


# function 6: format blast
all_frag_terminal_blast <- function(list)  {
  output_list <- list()
  for(k in seq_along(list))  {
    BED <- list[[k]][, c(1,3,4,7,8,9,10)]
    colnames(BED) <- c('read', 'Q', 'len', 'qs', 'qe', 'rs', 're')
    name <- as.data.frame(table(BED$read))
    split <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    fragment <- as.numeric(split[2:nrow(split),1])
    split <- BED[BED$read %in% name$Var1[name$Freq %in% 1],]
    split$frag <- 1
    split$s2 <- split$e2 <- ''
    split <- split[, c(1:3, 4:5, 9, 10, 6:7)]
    # B <- split[order(split$V4),]
    # colnames(B) <- c('chr', 'start1', 'end1', 'start2', 'end2', 'name', 'score', 'strand', 'frag')
    
    for(i in fragment) {
      split <- BED[BED$read %in% name$Var1[name$Freq %in% i],]
      split$frag <- i
      # split <- split[order(split$read),]
      D <- split[seq(i,nrow(split), i),1:7]
      C <- split[seq(1,nrow(split), i),c(1, 6:8)]
      C <- cbind(D,C)
      C[,9] <- i
      colnames(C) <- c('chr', 'start1', 'end1', 'start2', 'end2', 'name', 'score', 'strand', 'frag')
      B <- rbind(B,C)
    }
    # colnames(B) <- c('chr', 'start1', 'end1', 'start2', 'end2', 'name', 'score', 'strand', 'frag')
    output_list[[k]] <- B
    names(output_list)[k] <- names(list)[k]
  }
  return(output_list)
}

# Function3: identify the sgRNA (for files from Function5)
## gff_file have to be the required format
## sgmRNA_BA <- sgmRNA_boundary_allow
specific_DVG_class <- function(list, sgmRNA_BA=10, leader, UTR5, UTR3, frag, out_prefix=''
                               , length, write=F, write_path=getwd())
{ 
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  # if(frag==1)  {
  #   output_list <- list()
  #   for(i in seq_along(list))  {
  #     f1 <- list[[i]]
  #     f1$length <- f1[,3]-f1[,2]
  #     f1$complete <- f1$length > (length-whole_genome_tolerance)
  #     complete <- f1[f1$complete %in% TRUE,-c(ncol(f1))]
  #     if(nrow(complete)>0)  {
  #       complete$class <- c('complete genome')
  #     }
  #     if(nrow(complete)==0)  {
  #       complete <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
  #     
  #     f1 <- f1[f1$complete %in% FALSE,-c(ncol(f1))]
  #     f1$utr3 <- f1$end> UTR3
  #     delete_utr3 <- f1[f1$utr3 %in% FALSE,-c(ncol(f1))]
  #     if(nrow(delete_utr3)>0)  {
  #       delete_utr3$utr5 <- delete_utr3$start < UTR5
  #       delete_utr3_utr5 <- delete_utr3[delete_utr3$utr5 %in% FALSE, -c(ncol(delete_utr3))]
  #       if(nrow(delete_utr3_utr5)>0)  {
  #         delete_utr3_utr5$class <- c('Δ 3\' and 5\' genome')
  #       }
  #       if(nrow(delete_utr3_utr5)==0)
  #       {delete_utr3_utr5 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
  #       
  #       delete_utr3 <- delete_utr3[delete_utr3$utr5 %in% TRUE, -c(ncol(delete_utr3))]
  #       if(nrow(delete_utr3)>0)  {
  #         delete_utr3$class <- c('Δ 3\' genome')
  #       }
  #       if(nrow(delete_utr3)==0)
  #       {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
  #     }
  #     if(nrow(delete_utr3)==0)
  #     {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
  #     
  #     f1 <- f1[f1$utr3 %in% TRUE,-c(ncol(f1))]
  #     f1$utr5 <- f1$start < UTR5
  #     delete_utr5 <- f1[f1$utr5 %in% FALSE,-c(ncol(f1))]
  #     if(nrow(delete_utr5)>0)  {
  #       delete_utr5$class <- c('Δ leader sgmRNA')
  #     }
  #     if(nrow(delete_utr5)==0)
  #     {delete_utr5 <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
  #     
  #     near_complete <- f1[f1$utr5 %in% T,-c(ncol(f1))]
  #     if(nrow(near_complete)>0)  {
  #       near_complete$class <- c('near complete genome')
  #     }
  #     if(nrow(near_complete)==0)
  #     {near_complete <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
  #     
  #     f1 <- rbind(delete_utr3, delete_utr3_utr5, delete_utr5, near_complete)
  #     output_list[[i]] <- f1
  #     names(output_list)[[i]] <- names(list)[[i]]
  #   }
  #   
  #   for(j in 1:length(output_list))  {
  #     B <- output_list[[j]]
  #     rownames(B) <- seq(1,nrow(B))
  #     if(nrow(B[B$class %in% 'Δ leader sgmRNA', ])>0)  {
  #       A <- B[-c(as.numeric(rownames(B[B$class %in% 'Δ leader sgmRNA', ]))),]
  #       sgmRNA <- list()
  #       A$sub_class <- c('')
  #       D <- B[B$class %in% c('Δ leader sgmRNA'), ]
  #       if(nrow(D)>0)  {
  #         for(i in seq_len(nrow(gff)))  {
  #           D$sub_class_bulin <- (gff[i,4]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(gff[i,4]+sgmRNA_BA)
  #           E <- D[D$sub_class_bulin==T,]
  #           E <- E[,-c(ncol(E))]
  #           if(nrow(E)>0)  {
  #             E$sub_class <- gff[i,9]
  #           }
  #           sgmRNA[[i]] <- E
  #           D <- D[D$sub_class_bulin ==F,]
  #         }
  #         D <- D[D$sub_class_bulin ==F,]
  #         D <- D[,-c(ncol(D))]
  #         if(nrow(D)>0)
  #         {D$sub_class <- c('uncanonical sgmRNAs')}
  #         sgmRNA[[(length(sgmRNA)+1)]] <- D
  #         D <- ldply(sgmRNA, rbind)
  #         D <- rbind(A,D)
  #         output_list[[j]] <- D 
  #       }
  #     }
  #     else if(nrow(B[B$class %in% 'sgmRNA', ])==0)  {
  #       B$sub_class <- c('')
  #       B <- B[order(B$freq, decreasing = T), ]
  #       output_list[[j]] <- B
  #     }
  #   }
  # }
  
  output_list <- list()
  if(frag==2)  {
    for(i in seq_along(list))  {
      A <- list[[i]]
      B <- A
      
      B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
      UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
      if(nrow(UTR3_F)>0)  {
        C <- UTR3_F
        C$UTR5 <- C[,'start1'] <= UTR5
        D <- C[C$UTR5 ==T,-c(ncol(C))]
        if(nrow(D)>0)  {
          D$class <- c('DVG')
          D$sub_class <- c('Δ 3\' DVG')
        }
        if(nrow(D)==0)  {
          D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        E <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(E)>0)  {
          E$class <- c('DVG')
          E$sub_class <- c('Δ 5\' Δ3\' DVG')
        }
        if(nrow(E)==0) {
          E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        delete_UTR_Di <- rbind(D,E)
      }
      
      B <- B[B$UTR3 ==T,-c(ncol(B))]
      if(nrow(B)>0)  {
        C <- B
        C$UTR5 <- C[,'start1'] < UTR5
        delete_5UTR_Di <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(delete_5UTR_Di)>0)  {
          delete_5UTR_Di$class <- c('DVG')
          delete_5UTR_Di$sub_class <- c('Δ 5\' DVG')
        }
        if(nrow(delete_5UTR_Di)==0)  {
          delete_5UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
      }
      
      class <- rbind(delete_UTR_Di, delete_5UTR_Di)
      
      B <- A[!(A$V4 %in% class$V4), ]
      B$class <- 'DVG'
      B$sub_class <- "5' 3' DVG"
      
      output_list[[i]] <- rbind(B, class)
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  # if(frag>2)  {
  #   output_list <- list()
  #   for(i in seq_along(list))  {
  #     frag_other <- list[[i]]
  #     if(nrow(frag_other)>0)  {
  #       B <- frag_other
  #       dfgap <- as.data.frame(matrix(data = '', nrow = nrow(B), ncol = 1))
  #       for(j in 2:frag) {
  #         dfgap[,(j-1)] <- frag_other[,paste('end', j, sep = '')]-frag_other[,paste('start', (j-1), sep = '')]
  #       }
  #       dfgap$min <- apply(dfgap, MARGIN = 1,FUN = min)
  #       B$gap <- dfgap$min
  #       B$small_gap <- B[,'gap'] <gap
  #       small_gap <- B[B$small_gap ==T,-c(ncol(B))]
  #       if(nrow(small_gap)>0)  {
  #         small_gap$class <- c('small_gap')
  #         small_gap$sub_class <- c('')
  #       }
  #       if(nrow(small_gap)==0)  {
  #         small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #       }
  #       
  #       B <- B[B$small_gap ==F,-c(ncol(B))]
  #       B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
  #       UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
  #       if(nrow(UTR3_F)>0)  {
  #         C <- UTR3_F
  #         C$UTR5 <- C[,'start1'] <= UTR5
  #         D <- C[C$UTR5 ==T,-c(ncol(C))]
  #         if(nrow(D)>0) {
  #           D$class <- c('DVG')
  #           D$sub_class <- c('Δ 3\' DVG')
  #         }
  #         if(nrow(D)==0)  {
  #           D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #         }
  #         E <- C[C$UTR5 ==F,-c(ncol(C))]
  #         if(nrow(E)>0)  {
  #           E$class <- c('DVG')
  #           E$sub_class <- c('Δ 3\' 5\' DVG')
  #         }
  #         if(nrow(E)==0)  {
  #           E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #         }
  #         delete_UTR_Di <- rbind(D,E)
  #       }
  #       if(nrow(UTR3_F)==0)  {
  #         delete_UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #       }
  #       
  #       UTR3_T <- B[B$UTR3 ==T,-c(ncol(B))]
  #       if(nrow(UTR3_T)>0)  {
  #         C <- UTR3_T
  #         C$UTR5 <- C[,'start1'] > UTR5
  #         D <- C[C$UTR5 ==T,-c(ncol(C))]
  #         if(nrow(D)>0)  {
  #           D$class <- c('DVG')
  #           D$sub_class <- c('Δ 5\' DVG')
  #         }
  #         if(nrow(D)==0)  { 
  #           D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #         }
  #         E <- C[C$UTR5 ==F,-c(ncol(C))]
  #         if(nrow(E)>0)  {
  #           E$class <- c('DVG')
  #           E$sub_class <- c('')
  #         }
  #         if(nrow(E)==0)  {
  #           E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #         }
  #         delete_UTR5_Di <- rbind(D,E)
  #       }
  #       if(nrow(UTR3_T)==0) {
  #         delete_UTR5_Di <-  as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
  #       }
  #       class <- rbind(delete_UTR5_Di, delete_UTR_Di, small_gap)
  #     }
  #     output_list[[i]] <- class
  #     names(output_list)[[i]] <- names(list)[[i]]
  #   }
  # }
  
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
    return(output_list)
  }
  if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
    if(isFALSE(grepl(write_path, pattern = '&/'))){
      write_path <- paste(write_path,'/', sep = '')
    }
    lapply(seq_along(output_list)
           , function(i) 
             write.csv(output_list[[i]]
                       , file = paste0(write_path
                                       ,names(output_list[i]),'_', frag, 'fragment_', out_prefix, ".csv"),
                       row.names = FALSE, quote = F))
  }
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}

# Function3: identify the sgRNA (for files from Function5)
## gff_file have to be the required format
## sgmRNA_BA <- sgmRNA_boundary_allow
# specific_reads_sgRNA_identify <- function(list, gap=50, gff_path, sgmRNA_BA=10
#                                           , leader, UTR5, UTR3, frag, out_prefix=''
#                                           , whole_genome_tolerance=20, length
#                                           , isrecom=F, write=F, write_path=getwd())
# { 
#   ipak <- function(pkg){
#     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#     if (length(new.pkg)) 
#       install.packages(new.pkg, dependencies = TRUE)
#     sapply(pkg, require, character.only = TRUE)
#   }
#   packages <- c('plyr', 'dplyr')
#   ipak(packages)
#   
#   gff <- read.csv(file = gff_path, header = F)
#   
#   if(frag==1)  {
#     output_list <- list()
#     for(i in seq_along(list))  {
#       f1 <- list[[i]]
#       f1$length <- f1[,3]-f1[,2]
#       f1$complete <- f1$length > (length-whole_genome_tolerance)
#       complete <- f1[f1$complete %in% TRUE,-c(ncol(f1))]
#       if(nrow(complete)>0)  {
#         complete$class <- c('complete genome')
#       }
#       if(nrow(complete)==0)  {
#         complete <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
#       
#       f1 <- f1[f1$complete %in% FALSE,-c(ncol(f1))]
#       f1$utr3 <- f1$end> UTR3
#       delete_utr3 <- f1[f1$utr3 %in% FALSE,-c(ncol(f1))]
#       if(nrow(delete_utr3)>0)  {
#         delete_utr3$utr5 <- delete_utr3$start < UTR5
#         delete_utr3_utr5 <- delete_utr3[delete_utr3$utr5 %in% FALSE, -c(ncol(delete_utr3))]
#         if(nrow(delete_utr3_utr5)>0)  {
#           delete_utr3_utr5$class <- c('Δ 3\' and 5\' genome')
#         }
#         if(nrow(delete_utr3_utr5)==0)
#         {delete_utr3_utr5 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
#         
#         delete_utr3 <- delete_utr3[delete_utr3$utr5 %in% TRUE, -c(ncol(delete_utr3))]
#         if(nrow(delete_utr3)>0)  {
#           delete_utr3$class <- c('Δ 3\' genome')
#         }
#         if(nrow(delete_utr3)==0)
#         {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
#       }
#       if(nrow(delete_utr3)==0)
#       {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
#       
#       f1 <- f1[f1$utr3 %in% TRUE,-c(ncol(f1))]
#       f1$utr5 <- f1$start < UTR5
#       delete_utr5 <- f1[f1$utr5 %in% FALSE,-c(ncol(f1))]
#       if(nrow(delete_utr5)>0)  {
#         delete_utr5$class <- c('Δ leader sgmRNA')
#       }
#       if(nrow(delete_utr5)==0)
#       {delete_utr5 <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
#       
#       near_complete <- f1[f1$utr5 %in% T,-c(ncol(f1))]
#       if(nrow(near_complete)>0)  {
#         near_complete$class <- c('near complete genome')
#       }
#       if(nrow(near_complete)==0)
#       {near_complete <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
#       
#       f1 <- rbind(delete_utr3, delete_utr3_utr5, delete_utr5, near_complete)
#       output_list[[i]] <- f1
#       names(output_list)[[i]] <- names(list)[[i]]
#     }
#     
#     for(j in 1:length(output_list))  {
#       B <- output_list[[j]]
#       rownames(B) <- seq(1,nrow(B))
#       if(nrow(B[B$class %in% 'Δ leader sgmRNA', ])>0)  {
#         A <- B[-c(as.numeric(rownames(B[B$class %in% 'Δ leader sgmRNA', ]))),]
#         sgmRNA <- list()
#         A$sub_class <- c('')
#         D <- B[B$class %in% c('Δ leader sgmRNA'), ]
#         if(nrow(D)>0)  {
#           for(i in seq_len(nrow(gff)))  {
#             D$sub_class_bulin <- (gff[i,4]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(gff[i,4]+sgmRNA_BA)
#             E <- D[D$sub_class_bulin==T,]
#             E <- E[,-c(ncol(E))]
#             if(nrow(E)>0)  {
#               E$sub_class <- gff[i,9]
#             }
#             sgmRNA[[i]] <- E
#             D <- D[D$sub_class_bulin ==F,]
#           }
#           D <- D[D$sub_class_bulin ==F,]
#           D <- D[,-c(ncol(D))]
#           if(nrow(D)>0)
#           {D$sub_class <- c('uncanonical sgmRNAs')}
#           sgmRNA[[(length(sgmRNA)+1)]] <- D
#           D <- ldply(sgmRNA, rbind)
#           D <- rbind(A,D)
#           output_list[[j]] <- D 
#         }
#       }
#       else if(nrow(B[B$class %in% 'sgmRNA', ])==0)  {
#         B$sub_class <- c('')
#         B <- B[order(B$freq, decreasing = T), ]
#         output_list[[j]] <- B
#       }
#     }
#   }
#   
#   
#   if(frag==2)  {
#     output_list <- list()
#     for(i in seq_along(list))  {
#       A <- list[[i]]
#       B <- A
#       
#       B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
#       UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
#       if(nrow(UTR3_F)>0)  {
#         C <- UTR3_F
#         C$UTR5 <- C[,'start1'] <= UTR5
#         D <- C[C$UTR5 ==T,-c(ncol(C))]
#         if(nrow(D)>0)  {
#           D$class <- c('DVG')
#           D$sub_class <- c('Δ 3\' DVG')
#         }
#         if(nrow(D)==0)  {
#           D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
#         }
#         E <- C[C$UTR5 ==F,-c(ncol(C))]
#         if(nrow(E)>0)  {
#           E$class <- c('DVG')
#           E$sub_class <- c('Δ 5\' Δ3\' DVG')
#         }
#         if(nrow(E)==0) {
#           E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
#         }
#         delete_UTR_Di <- rbind(D,E)
#       }
#       
#       B <- B[B$UTR3 ==T,-c(ncol(B))]
#       if(nrow(B)>0)  {
#         C <- B
#         C$UTR5 <- C[,'start1'] < UTR5
#         delete_5UTR_Di <- C[C$UTR5 ==F,-c(ncol(C))]
#         if(nrow(delete_5UTR_Di)>0)  {
#           delete_5UTR_Di$class <- c('DVG')
#           delete_5UTR_Di$sub_class <- c('Δ 5\' DVG')
#         }
#         if(nrow(delete_5UTR_Di)==0)  {
#           delete_5UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
#         }
#       }
#       
#       C <- C[C$UTR5 ==T, -c(ncol(C))]
#       C$leader <- C$end1 <= leader+sgmRNA_BA
#       sgmRNA <- C[C$leader == T, -c(ncol(C))]
#       if(nrow(sgmRNA)>0)  {
#         sgmRNA$class <- c('sgmRNA')
#         sgmRNA$sub_class <- c('')
#       }
#       if(nrow(sgmRNA)==0) {sgmRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
#       }
#       
#       DiRNA <- C[C$leader ==F,-c(ncol(C))]
#       if(nrow(DiRNA)>0) {
#         DiRNA$class <- c('DVG')
#         DiRNA$sub_class <- c('5\' 3\' DVG')
#       }
#       if(nrow(DiRNA)==0) {
#         DiRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
#       }
#       class <- rbind(delete_UTR_Di, delete_5UTR_Di, DiRNA, sgmRNA)
#       
#       B <- class
#       rownames(B) <- seq(1,nrow(B))
#       if(nrow(B[B$class %in% 'sgmRNA', ])>0)  {
#         A <- B[!(B$class %in% 'sgmRNA'),]
#         sgmRNA <- list()
#         D <- B[B$class %in% c('sgmRNA'), ]
#         if(nrow(D)>0) {
#           for(j in seq_len(nrow(gff)))  {
#             D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,paste('start', 2, sep = '')] & D[,paste('start', 2, sep = '')]<(gff[j,4]+sgmRNA_BA)
#             E <- D[D$sub_class_bulin==T,]
#             E <- E[,-c(ncol(E))]
#             if(nrow(E)>0)  {
#               E$sub_class <- gff[j,9]
#             }
#             sgmRNA[[j]] <- E
#             D <- D[D$sub_class_bulin ==F,]
#           }
#           D <- D[D$sub_class_bulin ==F,]
#           D <- D[,-c(ncol(D))]
#           if(nrow(D)>0)  {
#             D$sub_class <- c('noncan sgmRNAs')
#           }
#           sgmRNA[[(length(sgmRNA)+1)]] <- D
#           D <- do.call(rbind, sgmRNA)
#           B <- rbind(A,D)
#           if(isTRUE(isrecom)) {
#             B <- B[order(B$freq, decreasing = T),]
#           }
#         }
#       }
#       if(nrow(B[B$sub_class %in% 'noncan sgmRNAs', ])>0)  {
#         A <- B[!(B$sub_class %in% 'noncan sgmRNAs'),]
#         sgmRNA <- list()
#         D <- B[B$sub_class %in% c('noncan sgmRNAs'), ]
#         if(nrow(D)>0) {
#           for(j in seq_len(nrow(gff)))  {
#             D$sub_class_bulin <- (gff[j,4])<D[,paste('start', 2, sep = '')] & D[,paste('start', 2, sep = '')]<(gff[j,5])
#             E <- D[D$sub_class_bulin==T,]
#             E <- E[,-c(ncol(E))]
#             if(nrow(E)>0)  {
#               E$sub_class_nonsgm <- gff[j,9]
#             }
#             sgmRNA[[j]] <- E
#             D <- D[D$sub_class_bulin ==F,]
#           }
#           D <- D[D$sub_class_bulin ==F,]
#           D <- D[,-c(ncol(D))]
#           if(nrow(D)>0)  {
#             D$sub_class_nonsgm <- c('out of frame')
#           }
#           sgmRNA[[(length(sgmRNA)+1)]] <- D
#           D <- do.call(rbind, sgmRNA)
#           A$sub_class_nonsgm <- ''
#           B <- rbind(A,D)
#           if(isTRUE(isrecom)) {
#             B <- B[order(B$freq, decreasing = T),]
#           }
#         }
#       }
#       output_list[[i]] <- B
#       names(output_list)[[i]] <- names(list)[[i]]
#     }
#   }
#   
#   if(frag>2)  {
#     output_list <- list()
#     for(i in seq_along(list))  {
#       frag_other <- list[[i]]
#       if(nrow(frag_other)>0)  {
#         B <- frag_other
#         dfgap <- as.data.frame(matrix(data = '', nrow = nrow(B), ncol = 1))
#         for(j in 2:frag) {
#           dfgap[,(j-1)] <- frag_other[,paste('end', j, sep = '')]-frag_other[,paste('start', (j-1), sep = '')]
#         }
#         dfgap$min <- apply(dfgap, MARGIN = 1,FUN = min)
#         B$gap <- dfgap$min
#         B$small_gap <- B[,'gap'] <gap
#         small_gap <- B[B$small_gap ==T,-c(ncol(B))]
#         if(nrow(small_gap)>0)  {
#           small_gap$class <- c('small_gap')
#           small_gap$sub_class <- c('')
#         }
#         if(nrow(small_gap)==0)  {
#           small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#         }
#         
#         B <- B[B$small_gap ==F,-c(ncol(B))]
#         B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
#         UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
#         if(nrow(UTR3_F)>0)  {
#           C <- UTR3_F
#           C$UTR5 <- C[,'start1'] <= UTR5
#           D <- C[C$UTR5 ==T,-c(ncol(C))]
#           if(nrow(D)>0) {
#             D$class <- c('DVG')
#             D$sub_class <- c('Δ 3\' DVG')
#           }
#           if(nrow(D)==0)  {
#             D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#           }
#           E <- C[C$UTR5 ==F,-c(ncol(C))]
#           if(nrow(E)>0)  {
#             E$class <- c('DVG')
#             E$sub_class <- c('Δ 3\' 5\' DVG')
#           }
#           if(nrow(E)==0)  {
#             E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#           }
#           delete_UTR_Di <- rbind(D,E)
#         }
#         if(nrow(UTR3_F)==0)  {
#           delete_UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#         }
#         
#         UTR3_T <- B[B$UTR3 ==T,-c(ncol(B))]
#         if(nrow(UTR3_T)>0)  {
#           C <- UTR3_T
#           C$UTR5 <- C[,'start1'] > UTR5
#           D <- C[C$UTR5 ==T,-c(ncol(C))]
#           if(nrow(D)>0)  {
#             D$class <- c('DVG')
#             D$sub_class <- c('Δ 5\' DVG')
#           }
#           if(nrow(D)==0)  { 
#             D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#           }
#           E <- C[C$UTR5 ==F,-c(ncol(C))]
#           if(nrow(E)>0)  {
#             E$class <- c('DVG')
#             E$sub_class <- c('')
#           }
#           if(nrow(E)==0)  {
#             E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#           }
#           delete_UTR5_Di <- rbind(D,E)
#         }
#         if(nrow(UTR3_T)==0) {
#           delete_UTR5_Di <-  as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
#         }
#         class <- rbind(delete_UTR5_Di, delete_UTR_Di, small_gap)
#         if(isTRUE(isrecom)) {class <- class[order(class$freq, decreasing = T),]
#         }
#       }
#       output_list[[i]] <- class
#       names(output_list)[[i]] <- names(list)[[i]]
#     }
#   }
#   
#   if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
#     return(output_list)
#   }
#   if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
#     if(isFALSE(grepl(write_path, pattern = '&/'))){
#       write_path <- paste(write_path,'/', sep = '')
#     }
#     lapply(seq_along(output_list)
#            , function(i) 
#              write.csv(output_list[[i]]
#                        , file = paste0(write_path
#                                        ,names(output_list[i]),'_', frag, 'fragment_', out_prefix, ".csv"),
#                        row.names = FALSE, quote = F))
#   }
#   lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
# }


# sgRNA identify class table
class_table <- function(list, libs, prop=T)  {
  class <- as.data.frame(matrix('', nrow = libs, ncol = 2))
  colnames(class) <- c('sgmRNA',  'DVG')
  index <- c('sgmRNA',  'DVG')
  for(i in seq_along(list)) {
    A <- list[[i]] 
    print(nrow(A))
    if(prop==T)  {
      A <- as.data.frame(prop.table(table(A$class)), stringsAsFactors = F)
      A[,2] <- sapply(A[,2], as.numeric)*100
    }
    if(prop==F)  {
      A <- as.data.frame(table(A$class), stringsAsFactors = F)
      A[,2] <- sapply(A[,2], as.numeric)
    }
    
    
    for(j in seq(1,2)) {
      value <- A$Freq[A$Var1 %in% index[j]]
      class[i,j] <- ifelse(identical(value, numeric(0)), yes = 0, no = value)
      rownames(class)[i] <- names(list)[i]
    }
  }
  return(class)
}


sgmRNA_subclass_table <- function(list, gff_path, prop=T, filtered=F
                                  , colN='sub_class') {
  gff <- read.csv(file = gff_path, header = F)
  class <- as.data.frame(matrix('', nrow = length(list),ncol = length(unique(gff$V9))+1))
  colnames(class) <- c(unique(gff$V9), 'noncan sgmRNAs')
  for(i in seq_along(list)) {
    A <- list[[i]] 
    if(isFALSE(filtered))  {
      A <- A[A$class %in% 'sgmRNA', ]
    }
    
    if(prop==T)  {
      A <- as.data.frame(prop.table(table(A[, colN])), stringsAsFactors = F)
      A[,2] <- sapply(A[,2], as.numeric)*100
    }
    if(prop==F)  {
      A <- as.data.frame(table(A[, colN]), stringsAsFactors = F)
      A[,2] <- sapply(A[,2], as.numeric)
    }
    for(j in seq_len(ncol(class))){
      value <- A$Freq[A$Var1 %in% colnames(class)[j]]
      class[i,j] <- ifelse(identical(value, numeric(0)), yes = 0, no = value)
      if(is.na(class[i,j])) {
        class[i,j] <- 0
      }
      rownames(class)[i] <- names(list)[i]
    }
  }
  return(class)
}


sgRNA_table <- function(list, gff_path, colname_subset
                        , sgRNA, colname_table) {
  gff <- read.csv(file = gff_path, header = F)
  class <- as.data.frame(matrix('', nrow = length(list),ncol = length(unique(gff$V9))+1))
  colnames(class) <- c(unique(gff$V9), 'noncan sgmRNAs')
  for(i in seq_along(list)) {
    A <- list[[i]] 
    A <- A[A[, colname_subset] %in% sgRNA, ]
    A <- as.data.frame(prop.table(table(A[, colname_table])), stringsAsFactors = F)
    for(j in seq_len(ncol(class))){
      value <- A$Freq[A$Var1 %in% colnames(class)[j]]
      class[i,j] <- ifelse(identical(value, numeric(0)), yes = 0, no = value*100)
      if(is.na(class[i,j])) {
        class[i,j] <- 0
      }
      rownames(class)[i] <- names(list)[i]
    }
  }
  return(class)
}

# diRNA subclass percentage
# DiRNA_subclass_table <- function(list) {
#   class <- as.data.frame(matrix('', nrow = 8,ncol = 4))
#   colnames(class) <- c("5' 3' Di", "delete 3'UTR Di", "delete 5' 3' Di", "delete 5' Di")
#   for(i in seq_along(list)) {
#     A <- list[[i]]
#     A <- A[A$class %in% 'DiRNA', ]
#     A <- as.data.frame(prop.table(table(A$sub_class)), stringsAsFactors = F)
#     for(j in seq_len(ncol(class))){
#       value <- A$Freq[A$Var1 %in% colnames(class)[j]]
#       class[i,j] <- ifelse(identical(value, numeric(0)), yes = 0, no = value*100)
#       if(is.na(class[i,j])) {
#         class[i,j] <- 0
#       }
#       rownames(class)[i] <- names(list)[i]
#     }
#   }
#   return(class)
# }
# 


# diRNA subclass percentage/ N containing proportion
# sum 2 replications
DiRNA_subclass_table <- function(list, prop=T, position=F) {
  class <- as.data.frame(matrix('', nrow = 0,ncol =4 ))
  colnames(class) <- c('libs', 'sub_class', 'percent', 'N_proportion')
  # c("5' 3' Di", "Δ 3'UTR Di", "Δ 5' 3' Di", "Δ 5' Di")
  for(i in seq_along(list)) {
    A <- list[[i]]
    A <- A[A$class %in% c('DVG', 'sgmRNA'), ]
    B <- as.list(A %>% dplyr::group_split(sub_class))
    for(j in seq_along(B))  {
      C <- as.data.frame(B[[j]])
      if(isTRUE(position))  {
        C[, 'N_containing'] <- ifelse(C[, 'start2'] >30740, yes = 0, no = 1)
      }
      if(isFALSE(position)){
        C[, 'N_containing'] <- ifelse(C[, 'Peptide'] =='', yes = 0, no = 1)
      }
      N <- colSums(as.matrix(C$N_containing))/nrow(C)
      names(B)[j] <- unique(C$sub_class)
      B[[j]] <- N
    }
    if(prop==T)  {
      A <- as.data.frame(prop.table(table(A$sub_class)), stringsAsFactors = F)
      A[,2] <- sapply(A[,2], as.numeric)*100
    }
    if(prop==F)  {
      A <- as.data.frame(table(A$sub_class), stringsAsFactors = F)
      A[,2] <- sapply(A[,2], as.numeric)
    }
    D <- as.data.frame(matrix('', nrow = 0,ncol = 4))
    colnames(D) <- c('libs', 'sub_class', 'percent', 'N_proportion')
    x=0
    for(j in c(names(B)))  {
      value <- A$Freq[A$Var1 %in% j ]
      x=x+1
      D[x, 'libs'] <- names(list)[i]
      D[x, 'sub_class'] <- j
      D[x, 'percent'] <- ifelse(identical(value, numeric(0)), yes = 0, no = value)
      D[x, 'N_proportion'] <- B[[j]]
    }
    class <- rbind(class, D)
  }
  colnames(class) <- c('libs', 'sub_class', 'Percentage', 'N_proportion')
  return(class)
}



cut_by_n_bp <- function(list, genome_length, bin=10, frag='', cunstom_index=F
                        , index, cols) {
  if(isTRUE(cunstom_index)){
    index <- index
  }
  if(isFALSE(cunstom_index)){
    index <- as.data.frame(matrix(data = 0, nrow = length(seq(0, 31032,bin)), ncol = 0))
    index[,1] <- seq(0, 31032, bin)
    index[nrow(index),1] <- 31032
  }
  
  out_list <- list()
  for(j in seq_along(list)) {
    A <- list[[j]]
    A$index <- seq(1, nrow(A))
    if(frag!='')  {
      colN <- c(paste('start', seq_len(frag), sep = ''), paste('end', seq_len(frag), sep = ''))
    }else (colN <- colnames(A)[cols])
    for(k in colN) {
      C <- as.data.frame(matrix(data = '', nrow = 0, ncol = 0))
      for(i in 2:nrow(index)) {
        A[,k] <- as.numeric(A[,k])
        A[,'clu'] <- A[,k] > index[i-1,1] & A[,k] <= index[i,1]
        row <- A$index[A$clu %in% T]
        if(length(row)>0) {
          B <- A[A$index %in% row,-c(ncol(A))]
          B[rownames(A[A$clu %in% T,]),k] <- index[i,1]
          C <- rbind(B,C)
          A <- A[!(A$index %in% row),-c(ncol(A))]
        }
      }
      A <- C
    }
    
    out_list[[j]] <- A
    names(out_list)[j] <- names(list)[j]
  }
  return(out_list)
}

tobin_specific_column <- function(df, genome_length, bin=10, column) {
  index <- as.data.frame(matrix(data = 0, nrow = length(seq(0, 31032, 10)), ncol = 0))
  index[,1] <- seq(0, 31032, 10)
  index[nrow(index),1] <- 31032
  
  A <- df
  A$index <- seq(1, nrow(A))
  colN <- colnames(A)[column]
  for(k in colN) {
    C <- as.data.frame(matrix(data = '', nrow = 0, ncol = 0))
    for(i in 2:nrow(index)) {
      A[,k] <- as.numeric(A[,k])
      A[,'clu'] <- A[,k] > index[i-1,1] & A[,k] <= index[i,1]
      row <- A$index[A$clu %in% T]
      if(length(row)>0) {
        B <- A[A$index %in% row,-c(ncol(A))]
        B[rownames(A[A$clu %in% T,]),k] <- index[i,1]
        C <- rbind(B,C)
        A <- A[!(A$index %in% row),-c(ncol(A))]
      }
    }
    A <- C
  }
  
  A <- A[,-c(ncol(A))]
  return(A)
}

top_lib_recom_matrix <- function(lib_list, frag, class, all=F, top_n=100, position=2) {
  recom_table_list <- list()
  recom_list <- list()
  inter_count <- ''
  x=0
  for(i in seq_along(lib_list)) {
    lib1 <- lib_list[[i]]
    if(identical(class, 'sgmRNA'))  {
      lib1 <- lib1[lib1$class %in% 'sgmRNA', ]
      lib1 <- lib1[!(lib1$sub_class %in% 'uncanonical sgmRNAs'), ]
    }
    if(identical(class, 'DVG'))  {
      lib1 <- lib1[lib1$class %in% 'DVG', ]
    }
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    if(position==4)  {colNames <-  colnames(lib1)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))]}
    if(position==2)  {colNames <-  colnames(lib1)[c(3:(frag*2),(frag*2+6), (frag*2+7))]}
    lib1$recom <-  do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    lib1 <- lib1[lib1$Freq >=2, ]
    if(isFALSE(all))  {
      lib1 <- lib1 %>% slice_max(order_by = Freq, n = top_n)
      recom_list[[i]] <- lib1$Var1
      recom_table_list[[i]] <- lib1
    }
    
    if(isTRUE(all))  {
      recom_list[[i]] <- lib1$Var1
      recom_table_list[[i]] <- lib1
    }
  }
  names(recom_list) <- names(lib_list)
  names(recom_table_list) <- names(lib_list)
  
  all_recom <- Reduce(union, recom_list)
  lib_matrix <- as.data.frame(matrix(data = 0, nrow = length(all_recom), ncol = 8))
  lib_matrix$index <- seq(1, nrow(lib_matrix))
  rownames(lib_matrix) <- all_recom; colnames(lib_matrix) <- names(lib_list)
  
  for(i in seq_along(recom_table_list)) {
    B <- recom_table_list[[i]]
    B <- B[order(B$Var1, decreasing = T), ]
    index <- lib_matrix[rownames(lib_matrix) %in% B$Var1, 9]
    C <- lib_matrix[lib_matrix[,9] %in% index, ]
    C <- C[order(rownames(C), decreasing = T), ]
    lib_matrix <- lib_matrix[!(lib_matrix[,9] %in% index), ]
    print(table(rownames(C)== B$Var1))
    C[,i] <- B$Freq
    lib_matrix <- rbind(lib_matrix,C)
  }
  lib_matrix <- lib_matrix[, -c(ncol(lib_matrix))]
}



shared_recom <- function(list, sgRNA, min, bi=T, coordinate=2, sgmRNA=F)  {
  ## for DiRNA
  recom_list <- list()
  recom_name <- list()
  for(i in seq_along(list))  {
    lib <- list[[i]]
    lib <- lib[lib$class %in% sgRNA,]
    if(coordinate==2)  {
      if(sgmRNA==T)  {
        lib$index <- paste(lib$sgmRNA, lib$end1, lib$start2
                           , sep = '_')
      }
      if(sgmRNA==F)  {
        lib$index <- paste(lib$sub_class, lib$end1, lib$start2,sep = '_')
      }
      recom <- as.data.frame(table(lib$index), stringsAsFactors = F)
      recom_list[[i]] <- recom[recom$Freq >=min, ]
      recom_name[[i]] <- recom$Var1[recom$Freq >=min]
    }
    if(coordinate==4)  {
      if(sgmRNA==T)  {
        lib$index <- paste(lib$sgmRNA, lib$start1, lib$end1, lib$start2, lib$end2
                           , sep = '_')
      }
      if(sgmRNA==F)  {
        lib$index <- paste(lib$sub_class,lib$start1, lib$end1, lib$start2, lib$end2,  sep = '_')
      }
      recom <- as.data.frame(table(lib$index), stringsAsFactors = F)
      recom_list[[i]] <- recom[recom$Freq >=min, ]
      recom_name[[i]] <- recom$Var1[recom$Freq >=min]
    }
  }
  names(recom_list) <- names(list)
  
  recom <- as.data.frame(Reduce(x = recom_name, f = union))
  # inter <- as.data.frame(Reduce(x = recom_list, f = intersect))
  
  for(i in seq_along(recom_list))  {
    lib <- recom_list[[i]]
    lib <- lib[order(lib[,1], decreasing = T), ]
    recom[,(i+1)] <- 0
    A <- recom[recom[,1] %in% lib[,1],]; A <- A[order(A[,1], decreasing = T), ]
    recom <- recom[!(recom[,1] %in% lib[,1]),]
    A[,(i+1)] <- lib[,2]
    recom <- rbind(A, recom)
  }
  colnames(recom)[2:(length(recom_list)+1)] <- names(recom_list)
  if(bi==TRUE)  {
    colnames(recom)[1] <- 'bi'
    recom[,2:5] <- ifelse(recom[,2:5]>0, yes = 1, no = 0)
  }
  if(bi==FALSE)  {
    colnames(recom)[1] <- 'Freq'
  }
  return(recom)
}


shared_recom_expression <- function(list, prop=T, sgRNA, min=1, bi=T
                                    , coordinate=2, sgmRNA=F, repro='')  {
  recom_df <- shared_recom(list = list, sgRNA = sgRNA, sgmRNA=sgmRNA
                           , min = min, bi = bi, coordinate = coordinate)
  if(colnames(recom_df)[1] == 'Freq')  {
    recom_df[,2:ncol(recom_df)] <- ifelse(recom_df[,2:ncol(recom_df)]>0, yes = 1, no = 0)
  }
  recom_df$sum <- rowSums(as.matrix(recom_df[,2:ncol(recom_df)]))
  if(identical('', repro))  {
    repro <- length(list)
  }
  shared <- recom_df[,1][recom_df$sum >= repro]
  
  recom_express <- data.frame(shared)
  recom_express$index <- seq(1, nrow(recom_express))
  recom_express <- recom_express[order(recom_express$shared, decreasing = T), ]
  recom_df <- shared_recom(list = list, sgRNA = sgRNA, sgmRNA=sgmRNA
                           , min = min, bi = F, coordinate = coordinate)
  # for(i in seq_along(list))  {
  #   A <- list[[i]]
  #   if(coordinate==2)  {
  #     if(sgmRNA==T)  {
  #       A$recom <- paste( A$sgmRNA,A$end1, A$start2, sep = '_')
  #     }
  #     if(sgmRNA==F)  {
  #       A$recom <- paste(A$sub_class,A$end1, A$start2,  sep = '_')
  #     }
  #   }
  #   if(coordinate==4)  {
  #     if(sgmRNA==T)  {
  #       A$recom <- paste( A$sgmRNA,A$start1, A$end1, A$start2, A$end2, sep = '_')
  #     }
  #     if(sgmRNA==F)  {
  #       A$recom <- paste(A$sub_class,A$start1, A$end1, A$start2, A$end2,  sep = '_')
  #     }
  #   }
  #   A <- A[A$recom %in% shared, ]
  #   if(prop==TRUE)  {
  #     A <- as.data.frame(prop.table(table(A$recom)))
  #     A <- A[order(A$Var1, decreasing = T), ]
  #     recom_express[, (i+2)] <- A$Freq*100
  #   }
  #   if(prop==FALSE)  {
  #     A <- as.data.frame(table(A$recom))
  #     A <- A[order(A$Var1, decreasing = T), ]
  #     recom_express[, (i+2)] <- A$Freq
  #   }
  # }
  
  for(i in 2:ncol(recom_df))  {
    A <- recom_df[recom_df$Freq %in% recom_express[,1], ]
    A <- A[order(A$Freq, decreasing = T), ]
    recom_express[,(i+1)] <- A[, i]
  }
  recom_express <- recom_express[,-c(2)]
  colnames(recom_express)[2:5] <- colnames(recom_df)[2:5]

  rownames(recom_express) <- recom_express[,1]
  recom_express <- recom_express[,-c(1)]
  if(identical(4, coordinate))  {
    recom_express[, 5:9] <- do.call(rbind, strsplit(rownames(recom_express), split = '_'))
    recom_express[, 6:9] <- apply(recom_express[, 6:9], MARGIN = 2, as.numeric) 
  }
  else if(identical(2, coordinate))  {
    recom_express[, 5:7] <- do.call(rbind, strsplit(rownames(recom_express), split = '_'))
    recom_express[, 6:7] <- apply(recom_express[, 6:7], MARGIN = 2, as.numeric)
  }
  recom_express[, 1:4] <- apply(recom_express[, 1:4], 2, as.numeric)
  if(isTRUE(prop))  {
    for(i in seq_along(list))  {
      recom_express[, i] <- recom_express[, i]/nrow(list[[i]])
    }
  }
  return(recom_express)
}


ORF_filter <- function(orf_list, bed_list, sgRNA)  {
  out_list <- list()
  for(i in seq_along(orf_list))  {
    lib <- names(orf_list)[i]
    if(sgRNA=='all')  {
      # A <- do.call(rbind, orf_list[[lib]])
      A <- orf_list[[lib]]
      A <- A[!(A$Peptide1 == ''& A$Peptide2 == ''), ]
      A$recom <- paste(A$start1, A$end1, A$start2, A$end2, sep = '_')
    }
    else if(sgRNA!='all') {
      # A <- orf_list[[lib]][[sgRNA]]
      A <- orf_list[[lib]]
      A <- A[!(A$Peptide1 == ''& A$Peptide2 == ''), ]
      A$recom <- paste(A$start1, A$end1, A$start2, A$end2, sep = '_')
    }
    B <- bed_list[[lib]]
    print(nrow(B))
    B$recom <- paste(B$start1, B$end1, B$start2, B$end2, sep = '_')
    B <- B[B$recom %in% A$recom, ]
    print(nrow(B))
    out_list[[i]] <- B
    print(nrow(lib))
    names(out_list)[i] <- lib
  }
  return(out_list)
}



recom_AT <- function(list, AT_region, thread=50, write=F, frag
                     , write_path, prefix='', string_length=15)  {
  AT_region[,c(3:5)] <- t(apply(AT_region[,c(3:5)], 1, as.numeric))
  AT_region <- AT_region[AT_region[,5] >=thread, ]
  
  seq_region <- function(x)  {seq(x, x+(string_length-1), 1)}
  region <- list()
  for(i in seq_len(nrow(AT_region)))  {
    region[[i]] <- seq_region(AT_region[i,'start'])
  }
  region <- Reduce(union, region)
  
  result <- as.data.frame(matrix(data = 0, nrow = length(list),ncol = 6))
  colnames(result) <- c('within_AT', 'total', 'percentage', 'libs', 'frag'
                        , 'AT')
  for(i in seq_along(list))  {
    A <- list[[i]]; A <- A[A$class %in% 'DVG', ]
    colN <- colnames(A)[3:(frag*2)]
    A <- unlist(A[,colN])
    result[i,'total'] <- length(A)
    result[i,'within_AT'] <- length(A[A %in% region])
    result[i,'libs'] <- names(list)[i]
    result[i,'frag'] <- paste('frag', frag, sep = '')
    result[i,'AT'] <- paste(thread, '%', sep = '')
  }
  result[,'percentage'] <- result[,'within_AT']/result[,'total']*100
  if(isTRUE(write))  {
    write.csv(result, file = paste(write_path, prefix, 'recom_AT.csv', sep = '')
              , quote = F, )
  }
  return(result)
}

depth_table <- function(depth_list, prop=T, filter=F, region_from, region_to)  {
  depth <- list()
  x=0
  for(i in seq(1, length(depth_list), 2))  {
    B <- depth_list[[i]]
    C <- depth_list[[i+1]]
    B$V3 <- as.numeric(B$V3)+as.numeric(C$V3)
    if(isTRUE(prop))  {
      B$V3 <- log2((B$V3)+1)
      # B$V3 <- (B$V3+abs(min(B$V3))+1)
    }
    if(isTRUE(filter))  {
      B <- B[region_from:region_to, ]
    }
    B$con <- sub(names(depth_list)[i], pattern = '[0-9]', replacement = '')
    x=x+1
    depth[[x]] <- B
    names(depth)[x] <- sub(names(depth_list)[i], pattern = '[0-9]', replacement = '')
  }
  
  return(do.call(rbind, depth))
}


# for all fragment 
all_frag_ident <- function(list)  {
  all <- all_frag_terminal(list)
  # all_1 <- list()
  # for(i in seq_along(all))  {
  #   all_1[[i]] <- all[[i]][all[[i]]$frag %in% 1,]
  #   all[[i]] <- all[[i]][!(all[[i]]$frag %in% 1),]
  # }
  # 
  # for(i in seq_along(all_1))  {
  #   df <- all_1[[1]]
  #   
  # }
  
  for(i in seq_along(all))  {
    df <- all[[i]]
    # df_2 <- df[df$frag %in% 2,]
    # df_3 <- df[!(df$frag %in% 2),]
    # df_2 <- specific_reads_sgRNA_identify(list = list(df_2), gap = 50, length = length
    #                                       , gff_path = gff_path, frag = 2
    #                                       , leader = leader, UTR5 = UTR5, UTR3 = UTR3
    #                                       , isrecom = F, write = F) %>% as.data.frame()
    # df_2$recom <-  do.call(paste, c(df_2[, c(2:5)], sep='_'))
    
    
    # df_3 <- specific_reads_sgRNA_identify(list = list(df_3), gap = 50
    #                                       , frag = 3, length = length
    #                                       , gff_path = gff_path
    #                                       , leader = leader, UTR5 = UTR5, UTR3 = UTR3
    #                                       , isrecom = T, write = F) %>% as.data.frame()
    # df_3$recom <-  do.call(paste, c(df_3[, c(9, 2:5)], sep='_'))
    df_allpos <- list()
    x=0
    for(j in unique(df$frag))  {
      x=x+1
      A <- specific_frag_format_organ(list = list[i], frag = j, min_gap = 201) %>% 
        as.data.frame()
      if(nrow(A) > 0)  {
        A$recom <- do.call(paste, c(A[, c(2:(2*j+1))], sep='_'))
        A <- A[, c((2*j+2), ncol(A))]
        colnames(A) <- c('reads', 'recom')
        df_allpos[[x]] <- A
      }
    }
    df_allpos <- do.call(rbind, df_allpos)
    df <- df[df$name %in% df_allpos[,1], ]
    df <- merge(df, df_allpos, by.x=c("name"), by.y=colnames(df_allpos)[1])
    
    # A <- grep(df_23$class, pattern = 'Di')
    # df_23[A, 'recom'] <- paste('D', df_23[A, 'recom'], sep='_')
    # df_23[A, 'RNA_recom'] <- paste('D', df_23[A, 'RNA_recom'], sep='_')
    # 
    # A <- grep(df_23$class, pattern = 'sg')
    # B <- grep(df_23$sub_class, pattern = 'non')
    # A <- A[!(A %in% B)]
    # df_23[A, 'recom'] <- paste('c_sg', df_23[A, 'recom'], sep='_')
    # df_23[A, 'RNA_recom'] <- paste('c_sg', df_23[A, 'RNA_recom'], sep='_')
    # df_23[B, 'recom'] <- paste('nc_sg', df_23[B, 'recom'], sep='_')
    # df_23[B, 'RNA_recom'] <- paste('nc_sg', df_23[B, 'RNA_recom'], sep='_')
    
    
    # B <- as.data.frame(table(df$recom), stringsAsFactors = F)
    # B <- B[order(B$Var1, decreasing = T), ]
    # df <- df %>% dplyr::distinct(recom, .keep_all = T)
    # df <- merge(df, B, by.x = c('recom'), by.y = c('Var1'))
    # df_23 <- df_23[, c(2:12, 1, 14, 13)]
    all[[i]] <- df
  }
  
  # for(i in seq_along(all_1))  {
  #   # colnames(all[[i]]) <- colnames(all_1[[i]])
  #   all[[i]] <- rbind(all[[i]], all_1[[i]])
  # }
  # all <- as.data.frame(do.call(rbind, all))
  return(all)
}


express_by_fc <- function(recom, number=30)  {
  # if(fc==T)  {
  #   n_count <- list()
  #   for(i in seq_along(list))  {
  #     A <- list[[i]]
  #     n_count[[i]] <- nrow(A[A$class %in% sgRNA & A$sub_class %in% sub_class, ])
  #   }
  #   n_count <- unlist(n_count)
  #   out <- sweep(recom ,2 , n_count,'/')
  #   out$mean <- rowMeans(out)
  #   out <- out[order(out$mean, decreasing = T), ]
  # }
  
  recom$mean <- rowMeans(recom)
  out <- recom %>% top_n(n = number, wt = mean)
  out <- out[order(out$mean, decreasing = T), ]
  dvg <- rownames(out)
  if(number<=nrow(out))  {
    out <- do.call(rbind, strsplit(rownames(out)[1:number], split = '_'))
  }
  if(number>nrow(out))  {
    out <- do.call(rbind, strsplit(rownames(out)[1:nrow(out)], split = '_'))
  }
  out <- as.data.frame(out)
  # out$V5 <- dvg
  # out$V5 <- paste(out$V5, '_', seq(1, nrow(out), 1), sep = '')
  return(out)
}

DVG_segment <- function(ref, recom, ORF_name='DVG')  {
  out <- as.data.frame(matrix('', nrow = 0, ncol = 9))
  colnames(out) <- colnames(ref)
  if(ORF_name=='DVG')  {
    recom$recom <- paste(recom[, 1], paste(recom[, 2], recom[, 3],recom[, 4]
                               , recom[, 5], sep = '_'), sep = ' ')
    out[seq(1, nrow(recom)*2, 2), 3] <- recom[,1]
    out[seq(2, nrow(recom)*2, 2), 3] <- recom[,3]
    out[seq(1, nrow(recom)*2, 2), 4] <- recom[,2]
    out[seq(2, nrow(recom)*2, 2), 4] <- recom[,4]
    out[seq(1, nrow(recom)*2, 2), 1] <- recom$recom 
    out[seq(2, nrow(recom)*2, 2), 1] <- recom$recom 
    out$ORF <- ORF_name
  }
  if(ORF_name=='sgm')  {
    recom$recom <- paste(recom[, 1], paste(recom[, 2], recom[, 3],recom[, 4]
                                           , recom[, 5], sep = '_'), sep = ' ')
    out[seq(1, nrow(recom)*2, 2), 3] <- recom[,2]
    out[seq(2, nrow(recom)*2, 2), 3] <- recom[,4]
    out[seq(1, nrow(recom)*2, 2), 4] <- recom[,3]
    out[seq(2, nrow(recom)*2, 2), 4] <- recom[,5]
    out[seq(1, nrow(recom)*2, 2), 1] <- recom$recom
    out[seq(2, nrow(recom)*2, 2), 1] <- recom$recom
    out$ORF <- ORF_name
  }
  
  # out$molecule <- factor(out$molecule, levels = unique(out$molecule))
  out$strand <- '+'
  out[,c(3,4)] <- apply(out[,c(3,4)], 2, as.integer)
  return(out)
}

library(ggplot2)
library(dplyr)
gg_theme= theme(
  title = element_text(size = 30),
  plot.subtitle = element_text(size = 26),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 24), #, angle=90,hjust=0.95,vjust=0.2
  axis.title.y = element_text(size = 28),
  axis.text.y = element_text(size = 20),
  plot.caption = element_text(size = 30), 
  legend.title=element_text(size=20),
  legend.text = element_text(size = 16), 
  legend.key.size = unit(2, 'lines'), 
  # legend.key.height = unit(5, "cm"),
  strip.text.x = element_text(size = 20), 
  strip.background = element_blank(), 
  text = element_text(size = 12))


gff <- read.csv(file = '~/Analysis/reference/bcov_edit_gff.csv', header = F)
gff <- gff[rownames(gff[!(gff$V9 %in% 'ORF1a'),]), ]
gff[(nrow(gff)+1):(nrow(gff)+3),] <- gff[1,]
rownames(gff) <- seq(1,nrow(gff))
gff[(nrow(gff)-2),c(4,5,9)] <- c(1, 70, 'leader')
gff[(nrow(gff)-1),c(4,5,9)] <- c(71, 210, "5'UTR")
gff[(nrow(gff)),c(4,5,9)] <- c(30741, 31032, "3'UTR")
gff$V4 <- as.numeric(gff$V4)
gff$V5 <- as.numeric(gff$V5)
gff <- gff[order(gff$V4, decreasing = F),]

gene_table <- gff[,c(1,9,4,5,6,7)]
gene_table[,5] <- 'forward'
colnames(gene_table) <- c('molecule', 'Gene', 'start', 'end', 'strand', 'direction')
gene_table$ORF <- factor(x = gene_table$Gene, levels = gene_table$Gene)
gene_table$gap <- gene_table$end-gene_table$start
gene_table[,'label'] <- gene_table[ ,'Gene']
colnames(gene_table)
# cut_frag2 <- read_as_list(path = '~/Analysis/bovine/viral_v2/', prefix = '.bed$') %>%
#   combine_rep() %>%
#   specific_frag_format_organ(frag = 2) %>%
#   specific_reads_sgRNA_identify(gap = 50, length = 31032, sgmRNA_BA = 10
#                                 , gff_path = '~/Analysis/reference/bovine_edit_gff.csv'
#                                 , leader = 64, UTR5 = 210, UTR3 = 30741, frag = 2) %>%
#   cut_by_n_bp(genome_length = 31032, bin = 10, frag = 2)


# frag2 <- read_as_list(path = '~/Analysis/bovine/viral_v2/', prefix = '.bed$') %>%
#   combine_rep() %>% 
#   specific_frag_format_organ(frag = 2) %>% 
#   specific_reads_sgRNA_identify(gap = 50, length = 31032, sgmRNA_BA = 10
#                                 , gff_path = '~/Analysis/reference/bovine_edit_gff.csv'
#                                 , leader = 64, UTR5 = 210, UTR3 = 30741, frag = 2)
