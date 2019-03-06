performance_comparison <- function(gene_only_SCRB,
                                   gene_filter,
                                   threshold,
                                   data_de){
  
  # Parameter in the function
  # gene_only_SCRB: the gene list without dropouts in SCRB
  # gene_filter: the gene list of the data
  # threshold: the threshold of genes
  # data_de: the combined data
  
  
  tmp_dropseq <- data_sc_Dropseq[match(gene_only_SCRB,gene_filter),]
  index <- rowSums(tmp_dropseq == 0) > threshold
  colnames(index) <- NULL
  common_diff <- gene_only_SCRB[which(index)]
  
  data_ks <- data_de[match(common_diff,gene_filter),]
  data_ks_drop <- data_ks[,dataType == 2]
  data_ks_SCRB <- data_ks[,dataType == 1]
  data_ks_drimpute <- data_ks[,dataType == 3]
  data_ks_scimpute <- data_ks[,dataType == 4]
  data_ks_magic <- data_ks[,dataType == 5]
  data_ks_viper <- data_ks[,dataType == 7]
  data_ks_scrabble <- data_ks[,dataType == 8]
  
  ks_value <- matrix(0,nrow = dim(data_ks)[1], ncol = 6)
  
  for (i in c(1:dim(data_ks)[1])){
    
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_drop[i,]))
    ks_value[i,1] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_drimpute[i,]))
    ks_value[i,2] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_scimpute[i,]))
    ks_value[i,3] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_magic[i,]))
    ks_value[i,4] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_viper[i,]))
    ks_value[i,5] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_scrabble[i,]))
    ks_value[i,6] <- tmp$statistic
  }
  
  tmp <- melt(ks_value)
  
  ksData <- data.frame(ks=tmp$value,group=tmp$Var2)
  my_comparisons <- list( c("1", "6"),c("2", "6"), c("3", "6"), c("4", "6"), c("5", "6"))
  pval1 <- compare_means(ks ~ group,data = ksData, method = "t.test", paired = T)
  pval <- c(paste0('p1 = ',formatC(pval1$p[5], format = "e", digits = 2)),
            paste0('p1 = ',formatC(pval1$p[9], format = "e", digits = 2)),
            paste0('p1 = ',formatC(pval1$p[12], format = "e", digits = 2)),
            paste0('p1 = ',formatC(pval1$p[14], format = "e", digits = 2)),
            paste0('p1 = ',formatC(pval1$p[15], format = "e", digits = 2)))
  pl <- ggboxplot(ksData, x = "group", y = "ks", fill = "group",
                  palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#ef8a62", "#6ebb00")) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    ylim(c(0,1.6)) + 
    theme_bw() +
    ggtitle(paste0("Genes used: ", length(common_diff), "  Dropout cutoff: ", formatC(threshold*100/76, digits = 0, format = "f"),"%")) +
    geom_signif(comparisons = my_comparisons, 
                annotations = pval,
                tip_length = 0.03,
                y_position = c(1.5,1.4,1.3, 1.2, 1.1)) +
    theme(text=element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Data Sets") + 
    ylab("K-S Statistics") +
    scale_fill_manual( values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07",  "#ef8a62", "#6ebb00"),
                       name="Method",
                       breaks=c("1", "2", "3", "4", "5", "6"),
                       labels=c("Dropout", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  
  return(pl)
  
}

performance_distribution_comparison <- function(gene_only_SCRB,
                                                gene_filter,
                                                threshold,
                                                data_de){
  
  # Parameter in the function
  # gene_only_SCRB: the gene list without dropouts in SCRB
  # gene_filter: the gene list of the data
  # threshold: the threshold of genes
  # data_de: the combined data
  
  tmp_dropseq <- data_sc_Dropseq[match(gene_only_SCRB,gene_filter),]
  index <- rowSums(tmp_dropseq == 0) > threshold
  colnames(index) <- NULL
  common_diff <- gene_only_SCRB[which(index)]
  
  data_ks <- data_de[match(common_diff,gene_filter),]
  data_ks_drop <- data_ks[,dataType == 2]
  data_ks_SCRB <- data_ks[,dataType == 1]
  data_ks_drimpute <- data_ks[,dataType == 3]
  data_ks_scimpute <- data_ks[,dataType == 4]
  data_ks_magic <- data_ks[,dataType == 5]
  data_ks_viper <- data_ks[,dataType == 7]
  data_ks_scrabble <- data_ks[,dataType == 8]
  
  p <- list()
  for (i in c(1:dim(data_ks)[1])){
    tmp1 <- data.frame( value = t(data_ks_SCRB[i,]))
    tmp1$e <- 7
    tmp2 <- data.frame( value = t(data_ks_drop[i,]))
    tmp2$e <- 1
    tmp3 <- data.frame( value = t(data_ks_drimpute[i,]))
    tmp3$e <- 2
    tmp4 <- data.frame( value = t(data_ks_scimpute[i,]))
    tmp4$e <- 3
    tmp5 <- data.frame( value = t(data_ks_magic[i,]))
    tmp5$e <- 4
    tmp6 <- data.frame( value = t(data_ks_viper[i,]))
    tmp6$e <- 5
    tmp7 <- data.frame( value = t(data_ks_scrabble[i,]))
    tmp7$e <- 6
    tmp <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)
    colnames(tmp) <- c("value","e")
    data <- tmp
    mu <- ddply(data, "e", summarise, grp.mean=mean(value))
    p[[i]] <- ggplot(data, aes(x = value,color = as.factor(e), fill=as.factor(e))) +
      geom_density( alpha = 0.2,adjust = 10) +
      theme_bw() + 
      scale_fill_manual(labels = c("Dropseq","DrImpute","scImpute","MAGIC","VIPER","SCRABBLE","SCRBseq"),
                        values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#ef8a62", "#6ebb00","#8B008B")) + 
      scale_color_manual(labels = c("Dropseq","DrImpute","scImpute","MAGIC","VIPER","SCRABBLE","SCRBseq"),
                         values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#ef8a62", "#6ebb00","#8B008B")) +
      ggtitle(paste0(common_diff[i])) +
      theme(legend.position = "bottom",
            legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x="log10(Normalized Expression levels + 1)", y="Density")+
      scale_x_continuous(limits=c(-1,5)) +
      scale_y_continuous(limits = c(0,1))
  }
  
  # return the handler of the plot
  return(p)

}


performance_mean <- function(gene_only_SCRB,gene_filter,threshold, data_sc_Dropseq, data_de){
  
  tmp_dropseq <- data_sc_Dropseq[match(gene_only_SCRB,gene_filter),]
  index <- rowSums(tmp_dropseq == 0) > threshold
  colnames(index) <- NULL
  common_diff <- gene_only_SCRB[which(index)]
  
  data_ks <- data_de[match(common_diff,gene_filter),]
  data_ks_drop <- data_ks[,dataType == 2]
  data_ks_SCRB <- data_ks[,dataType == 1]
  data_ks_drimpute <- data_ks[,dataType == 3]
  data_ks_scimpute <- data_ks[,dataType == 4]
  data_ks_magic <- data_ks[,dataType == 5]
  data_ks_viper <- data_ks[,dataType == 7]
  data_ks_scrabble <- data_ks[,dataType == 8]
  
  ks_value <- matrix(0,nrow = dim(data_ks)[1], ncol = 7)
  
  for (i in c(1:dim(data_ks)[1])){
    
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_drop[i,]))
    ks_value[i,1] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_drimpute[i,]))
    ks_value[i,2] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_scimpute[i,]))
    ks_value[i,3] <- tmp$statistic
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_magic[i,]))
    ks_value[i,4] <- tmp$statistic
    
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_viper[i,]))
    ks_value[i,5] <- tmp$statistic
    
    tmp <- ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_scrabble[i,]))
    ks_value[i,6] <- tmp$statistic
    
    ks_value[i,7] <- mean(as.matrix(data_ks_SCRB[i,]))
  }
  
  ks_value <- ks_value[order(ks_value[,7]),]
  
  tmp <- melt(ks_value[,1:6])
  
  print(c(ks_value[1,7],ks_value[20,7], ks_value[40,7], ks_value[56,7]))
  
  ksData <- data.frame(ks=tmp$value,group=tmp$Var2, xvalue = rep(1:length(ks_value[,7]),6))
  
  
  
  pl <- ggplot(ksData, aes(x=xvalue, y=ks)) +
    geom_vline(xintercept= 1:length(ks_value[,7]),  colour="grey",  alpha=0.4) +
    geom_jitter(aes(colour = as.factor(group)), width = 0) +
    coord_fixed(ratio = 15) +
    scale_color_manual( values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#ef8a62", "#6ebb00"),
                        name="Method",
                        breaks=c("1", "2", "3", "4", "5", "6"),
                        labels=c("Dropout", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")) +
    ylim(c(0,1)) + 
    xlim(c(1,56)) +
    theme_bw() +
    ggtitle(paste0("Genes used: ", length(common_diff), "  Dropout cutoff: ", formatC(threshold*100/76, digits = 0, format = "f"),"%")) +
    theme(text=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position="bottom")+
    xlab("Mean of log10(Normalized Data + 1)") + 
    ylab("K-S Statistics")
  
  return(pl)
}
