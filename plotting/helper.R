library(grid)
library(gridExtra)
library(ggplot2)
library(tools)
library(tidyr)
library(hash)

# Borrowed from: https://rpubs.com/sjackman/grid_arrange_shared_legend
# Thanks to Shaun Jackman
grid_arrange_shared_legend <- function(show_legend, num_rows, legend_plot_id, ...) {
  plots <- list(...)
  if (show_legend) {
    g <- ggplotGrob(plots[[legend_plot_id]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    plots_arg <- lapply(plots, function(x)  x + theme(legend.position="none"))
    plots_arg$nrow <- num_rows
    grid.arrange(
      do.call(arrangeGrob, plots_arg),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight))
  }
  else {
    g <- ggplotGrob(plots[[legend_plot_id]] + theme(legend.position="none"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    plots_arg <- lapply(plots, function(x)  x + theme(legend.position="none"))
    plots_arg$nrow <- num_rows
    grid.arrange(
      do.call(arrangeGrob, plots_arg),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc")))
    # grid.arrange(
    #   do.call(arrangeGrob, lapply(plots, function(x)
    #     x + theme(legend.position="none"))),
    #   ncol = 1,
    #   heights = unit.c(unit(1, "npc")))
    
    # g <- ggplotGrob(plots[[legend_plot_id]] + theme(legend.position="bottom"))$grobs
    # legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    # grid.arrange(
    #   do.call(arrangeGrob, lapply(plots, function(x)
    #     x + theme(legend.position="none"))),
    #   ncol = 1,
    #   heights = unit.c(unit(1, "npc")))
  }
}

# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# Plot hap.py ROCs and PASS points
#
# This script requires the ggplot2 package.
# To install it, run this command in R:
#
#   install.packages(c("ggplot2"))
#
# Usage:
#
# This script runs on a set of hap.py results. Result n will have been run
# with hap.py -o prefix_n. The names for each result are optional, they can
# be used to specify a custom label for the ROCs in the plot.
#
# run Rscript rocplot.Rscript [-pr] output_name prefix_1:name_1 ... prefix_n:name_n
#
# Use the -pr switch to produce a precision-recall curve rather than a TPR/FPR curve
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#
# Read single hap.py / xcmp dataset
read_single = function(x) {
  nx = strsplit(x, "\\:")[[1]]
  
  if(length(nx) == 1) {
    name = basename(file_path_sans_ext(x))
  } else {
    x = nx[1]
    name = nx[2]
  }
  # cat(sprintf("Reading %s as %s\n", x, name))
  
  all_results = list()
  
  # all_results$roc_data_snp_all = read.csv(paste(x, "roc.Locations.SNP", "csv", "gz", sep="."))
  # all_results$roc_data_indel_all = read.csv(paste(x, "roc.Locations.INDEL", "csv", "gz", sep="."))
  
  all_results$roc_data_snp_pass = read.csv(paste(x, "roc.Locations.SNP.PASS", "csv", "gz", sep="."))
  all_results$roc_data_indel_pass = read.csv(paste(x, "roc.Locations.INDEL.PASS", "csv", "gz", sep="."))
  
  sel_snp_file = paste(x, "roc.Locations.SNP.SEL", "csv", "gz", sep=".")
  
  # if we have a selectively-filtered ROC, don't show "PASS" ROC
  if(file.exists(sel_snp_file)) {
    all_results$roc_data_snp_sel = read.csv(sel_snp_file)
  } else {
    all_results$roc_data_snp_sel = all_results$roc_data_snp_pass
  }
  
  # all_results$roc_data_snp_all = head(subset(all_results$roc_data_snp_all,
  #                                            QQ == min(all_results$roc_data_snp_all["QQ"])),
  #                                     n=1)
  # all_results$roc_data_snp_pass = head(subset(all_results$roc_data_snp_pass,
  #                                            QQ == min(all_results$roc_data_snp_pass["QQ"])),
  #                                      n=1)
  all_results$roc_data_snp_pass = head(subset(all_results$roc_data_snp_pass,
                                       METRIC.F1_Score == max(all_results$roc_data_snp_pass["METRIC.F1_Score"])),
                                       n=1)
  sel_min = head(subset(all_results$roc_data_snp_sel,
                        QQ == min(all_results$roc_data_snp_sel["QQ"])),
                 n=1)
  
  # all_results$roc_connector_snp =
  #     rbind(all_results$roc_data_snp_all, sel_min)
  all_results$roc_connector_snp = sel_min
  all_results$roc_connector_snp$Filter = "CONN"
  all_results$roc_data_snp_sel$Filter = "ROC"
  
  sel_indel_file = paste(x, "roc.Locations.INDEL.SEL", "csv", "gz", sep=".")
  if(file.exists(sel_indel_file)) {
    all_results$roc_data_indel_sel = read.csv(sel_indel_file)
  } else {
    # use PASS ROC if no SEL ROC present
    all_results$roc_data_indel_sel = all_results$roc_data_indel_pass
  }
  
  # just keep single ALL and PASS point
  # all_results$roc_data_indel_all = head(subset(all_results$roc_data_indel_all,
  #                                              QQ == min(all_results$roc_data_indel_all["QQ"])),
  #                                       n=1)
  all_results$roc_data_indel_pass = head(subset(all_results$roc_data_indel_pass,
                                                QQ == min(all_results$roc_data_indel_pass["QQ"])),
                                         n=1)
  
  sel_min = head(subset(all_results$roc_data_indel_sel,
                        QQ == min(all_results$roc_data_indel_sel["QQ"])),
                 n=1)
  
  # all_results$roc_connector_indel =
  #     rbind(all_results$roc_data_indel_all, sel_min)
  all_results$roc_connector_indel = sel_min
  all_results$roc_connector_indel$Filter = "CONN"
  all_results$roc_data_indel_sel$Filter = "ROC"
  
  result = do.call(rbind, all_results)
  row.names(result) = NULL
  result$filename = x
  # result$name = name
  result$name = strsplit(name, "-")[[1]][1]
  # print(name)
  # prefix <- data.frame(x=name) %>% separate(x, c("A"), sep="-")
  # print(prefix)
  # result$name = c(prefix)
  
  result$igroup = paste(result$name,
                        result$Filter,
                        result$Type)
  return(result)
}

read_single_with_map = function(x, map_caller_name) {
  nx = strsplit(x, "\\:")[[1]]
  
  if(length(nx) == 1) {
    name = basename(file_path_sans_ext(x))
  } else {
    x = nx[1]
    name = nx[2]
  }
  # cat(sprintf("Reading %s as %s\n", x, name))
  
  all_results = list()
  
  # all_results$roc_data_snp_all = read.csv(paste(x, "roc.Locations.SNP", "csv", "gz", sep="."))
  # all_results$roc_data_indel_all = read.csv(paste(x, "roc.Locations.INDEL", "csv", "gz", sep="."))
  
  all_results$roc_data_snp_pass = read.csv(paste(x, "roc.Locations.SNP.PASS", "csv", "gz", sep="."))
  all_results$roc_data_indel_pass = read.csv(paste(x, "roc.Locations.INDEL.PASS", "csv", "gz", sep="."))
  
  sel_snp_file = paste(x, "roc.Locations.SNP.SEL", "csv", "gz", sep=".")
  
  # if we have a selectively-filtered ROC, don't show "PASS" ROC
  if(file.exists(sel_snp_file)) {
    all_results$roc_data_snp_sel = read.csv(sel_snp_file)
  } else {
    all_results$roc_data_snp_sel = all_results$roc_data_snp_pass
  }
  
  # all_results$roc_data_snp_all = head(subset(all_results$roc_data_snp_all,
  #                                            QQ == min(all_results$roc_data_snp_all["QQ"])),
  #                                     n=1)
  all_results$roc_data_snp_pass = head(subset(all_results$roc_data_snp_pass,
                                              QQ == min(all_results$roc_data_snp_pass["QQ"])),
                                       n=1)
  # all_results$roc_data_snp_pass = head(subset(all_results$roc_data_snp_pass,
  #                                      METRIC.F1_Score == max(all_results$roc_data_snp_pass["METRIC.F1_Score"])),
  #                                      n=1)
  # sel_min = head(subset(all_results$roc_data_snp_sel,
  #                       QQ == min(all_results$roc_data_snp_sel["QQ"])),
  #                n=1)
  sel_min = head(subset(all_results$roc_data_snp_sel,
                        METRIC.F1_Score == max(all_results$roc_data_snp_sel["METRIC.F1_Score"])),
                 n=1)
  
  # all_results$roc_connector_snp =
  #     rbind(all_results$roc_data_snp_all, sel_min)
  all_results$roc_connector_snp = sel_min
  all_results$roc_connector_snp$Filter = "CONN"
  all_results$roc_data_snp_sel$Filter = "ROC"
  
  sel_indel_file = paste(x, "roc.Locations.INDEL.SEL", "csv", "gz", sep=".")
  if(file.exists(sel_indel_file)) {
    all_results$roc_data_indel_sel = read.csv(sel_indel_file)
  } else {
    # use PASS ROC if no SEL ROC present
    all_results$roc_data_indel_sel = all_results$roc_data_indel_pass
  }
  
  # just keep single ALL and PASS point
  # all_results$roc_data_indel_all = head(subset(all_results$roc_data_indel_all,
  #                                              QQ == min(all_results$roc_data_indel_all["QQ"])),
  #                                       n=1)
  all_results$roc_data_indel_pass = head(subset(all_results$roc_data_indel_pass,
                                                QQ == min(all_results$roc_data_indel_pass["QQ"])),
                                         n=1)
  
  # sel_min = head(subset(all_results$roc_data_indel_sel,
  #                       QQ == min(all_results$roc_data_indel_sel["QQ"])),
  #                n=1)
  sel_min = head(subset(all_results$roc_data_indel_sel,
                        METRIC.F1_Score == max(all_results$roc_data_indel_sel["METRIC.F1_Score"])),
                 n=1)
  
  # all_results$roc_connector_indel =
  #     rbind(all_results$roc_data_indel_all, sel_min)
  all_results$roc_connector_indel = sel_min
  all_results$roc_connector_indel$Filter = "CONN"
  all_results$roc_data_indel_sel$Filter = "ROC"
  
  result = do.call(rbind, all_results)
  row.names(result) = NULL
  result$filename = x
  # result$name = name
  result$name = strsplit(name, "-")[[1]][1]
  # print(name)
  # prefix <- data.frame(x=name) %>% separate(x, c("A"), sep="-")
  # print(prefix)
  # result$name = c(prefix)
  
  # Map names
  result$name <- values(map_caller_name, keys=result$name)
  
  result$igroup = paste(result$name,
                        result$Filter,
                        result$Type)
  return(result)
}

# Plot P/R curves
plot_data = function(pdata, highlight=NULL, is.PR=FALSE) {
  # precision / recall curve
  if(is.PR) {
    xaxis = "METRIC.Recall"
    yaxis = "METRIC.Precision"
  } else {
    # approximate ROC-style curve (FPR is not correct)
    # xaxis = "FPR"
    # yaxis = "TPR"
    # pdata$FPR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    # pdata$TPR = pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    xaxis = "FNR"
    yaxis = "FDR"
    pdata$FDR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    pdata$FNR = 1 - pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    cc = complete.cases(pdata[, c(xaxis, yaxis)])
    pdata = pdata[cc, ]
  }
  
  plt = ggplot(pdata, aes_string(x=xaxis, y=yaxis, color="name"))
  facet_wrap(~Type)
  
  pdata <- pdata %>% arrange(METRIC.Recall)
  if (is.null(highlight)) {
    # ROC lines
    plt = plt +
      geom_path(data = subset(pdata, Filter == "ROC"),
                mapping=aes(group=igroup),
                size=1.5,
                linetype=1)
  } else {
    plt = plt +
      geom_line(data = subset(pdata, Filter == "ROC"),
                mapping=aes(group=igroup),
                size=0,
                linetype="dashed")
    plt = plt +
      geom_line(data = pdata %>% filter(Filter == "ROC", name %in% highlight),
                mapping=aes(group=igroup),
                size=1,
                linetype=1)
    plt = plt +
      geom_line(data = pdata %>% filter(Filter == "ROC", !(name %in% highlight)),
                # data = subset(pdata, Filter == "ROC" & !(name %in% highlight)) %>% arrange(name_rev),
                mapping=aes(group=igroup),
                size=0.5,
                linetype="dashed")
  }
  
  # Connector between ALL and start of ROC
  plt = plt +
    geom_line(data = subset(pdata, Filter == "CONN"),
              mapping=aes(group=igroup),
              size=1.5,
              linetype=4)
  
  plt = plt +
    geom_point(data = subset(pdata, Filter %in% c("CONN")),
               mapping=aes(group=igroup),
               size=4)
  
  # plt = plt +
  #     geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
  #                mapping=aes(shape=Filter, group=igroup),
  #                size=8)
  # plt = plt +
  #     geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
  #                mapping=aes(group=igroup),
  #                size=8)
  
  xl_min = max(0,
               min(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) - 0.02)
  xl_max = min(1.0,
               max(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) + 0.02)
  yl_min = max(0,
               min(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) - 0.01)
  yl_max = min(1.0,
               max(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) + 0.01)
  
  plt = plt +
    scale_color_brewer("", palette="Set1") +
    xlim(c(xl_min, xl_max)) +
    ylim(c(yl_min, yl_max)) +
    theme_bw(base_size=18)
  
  return(plt)
}

plot_data_combined = function(
  pdata, highlight=NULL, is.PR=FALSE, xl_min=NULL, xl_max=NULL, yl_min=NULL, yl_max=NULL) {
  # precision / recall curve
  if(is.PR) {
    xaxis = "METRIC.Recall"
    yaxis = "METRIC.Precision"
  } else {
    # approximate ROC-style curve (FPR is not correct)
    # xaxis = "FPR"
    # yaxis = "TPR"
    # pdata$FPR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    # pdata$TPR = pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    xaxis = "FNR"
    yaxis = "FDR"
    pdata$FDR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    pdata$FNR = 1 - pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    cc = complete.cases(pdata[, c(xaxis, yaxis)])
    pdata = pdata[cc, ]
  }
  
  plt = ggplot(pdata, aes_string(x=xaxis, y=yaxis, color="name"))
  facet_wrap(~Type)
  
  if (is.null(highlight)) {
    # ROC lines
    plt = plt +
      geom_line(data = subset(pdata, Filter == "ROC"),
                mapping=aes(group=igroup),
                size=1,
                linetype=1)
  } else {
    plt = plt +
      geom_line(data = subset(pdata, Filter == "ROC"),
                mapping=aes(group=igroup),
                size=0,
                linetype="dashed")
    plt = plt +
      geom_line(data = pdata %>% filter(Filter == "ROC", name %in% highlight),
                mapping=aes(group=igroup),
                size=1,
                linetype=1)
    plt = plt +
      geom_line(data = pdata %>% filter(Filter == "ROC", !(name %in% highlight)),
                mapping=aes(group=igroup),
                size=0.5,
                linetype="dashed")
    # plt = plt +
    #   geom_line(data = subset(pdata, Filter == "ROC" & name %in% highlight),
    #             mapping=aes(group=igroup),
    #             size=1,
    #             linetype=1)
  }
  
  # Connector between ALL and start of ROC
  plt = plt +
    geom_line(data = subset(pdata, Filter == "CONN"),
              mapping=aes(group=igroup),
              size=1,
              linetype=4)
  
  plt = plt +
    geom_point(data = subset(pdata, Filter %in% c("CONN")),
               mapping=aes(group=igroup),
               size=2)
  
  # plt = plt +
  #     geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
  #                mapping=aes(shape=Filter, group=igroup),
  #                size=8)
  # plt = plt +
  #     geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
  #                mapping=aes(group=igroup),
  #                size=8)
  
  xl_min_raw <- min(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]])
  xl_max_raw <- max(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]])
  xl_extend <- (xl_max_raw - xl_min_raw) * 0.05
  yl_min_raw <- min(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]])
  yl_max_raw <- max(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]])
  yl_extend <- (yl_max_raw - yl_min_raw) * 0.1

  if (is.null(xl_min)){
    xl_min = max(0,
                 min(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) - xl_extend)
  }
  if (is.null(xl_max)){
    xl_max = min(1.0,
                 max(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) + xl_extend)
  }
  if (is.null(yl_min)){
    yl_min = max(0,
                 min(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) - yl_extend)
  }
  if (is.null(yl_max)){
    yl_max = min(1.0,
                 max(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) + yl_extend)
  }
  
  plt = plt +
    scale_color_brewer("", palette="Set1") +
    xlim(c(xl_min, xl_max)) +
    ylim(c(yl_min, yl_max)) +
    theme_bw(base_size=18)
  
  # plt = plt + facet_wrap(Type ~ Coverage)#, scales = "free")
  
  return(plt)
}

plot_data_commonness = function(pdata, highlight=NULL, is.PR=FALSE) {
  # precision / recall curve
  if(is.PR) {
    xaxis = "METRIC.Recall"
    yaxis = "METRIC.Precision"
  } else {
    # approximate ROC-style curve (FPR is not correct)
    # xaxis = "FPR"
    # yaxis = "TPR"
    # pdata$FPR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    # pdata$TPR = pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    xaxis = "FNR"
    yaxis = "FDR"
    pdata$FDR = pdata$QUERY.FP / (pdata$QUERY.TOTAL - pdata$QUERY.UNK)
    pdata$FNR = 1 - pdata$TRUTH.TP / (pdata$TRUTH.TP + pdata$TRUTH.FN)
    cc = complete.cases(pdata[, c(xaxis, yaxis)])
    pdata = pdata[cc, ]
  }
  
  plt = ggplot(pdata, aes_string(x=xaxis, y=yaxis, color="name"))
  facet_wrap(~Type)
  
  pdata <- pdata %>% arrange(METRIC.Recall)
  # ROC lines
  # plt = plt +
  #   geom_path(data = subset(pdata, Filter == "ROC"),
  #             mapping=aes(group=igroup),
  #             size=1.5,
  #             linetype=1)
  if (is.null(highlight)) {
    plt = plt +
      geom_point(data = subset(pdata, Filter == "ROC"),
                mapping=aes(group=igroup),
                size=1.5)
  } else {
    plt = plt +
      geom_point(data = subset(pdata, Filter == "ROC"),
                 mapping=aes(group=igroup),
                 size=0)
    plt = plt +
      geom_point(data = subset(pdata, Filter == "ROC" & !(name %in% highlight)),
                 mapping=aes(group=igroup),
                 size=0.1)
    plt = plt +
      geom_point(data = subset(pdata, Filter == "ROC" & (name %in% highlight)),
                 mapping=aes(group=igroup),
                 size=1.5)
  }

  # Connector between ALL and start of ROC
  # plt = plt +
  #   geom_line(data = subset(pdata, Filter == "CONN"),
  #             mapping=aes(group=igroup),
  #             size=0.5,
  #             linetype=1)
  
  plt = plt +
    geom_point(data = subset(pdata, Filter %in% c("CONN")),
               mapping=aes(group=igroup),
               size=4)
  
  # plt = plt +
  #     geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
  #                mapping=aes(shape=Filter, group=igroup),
  #                size=8)
  # plt = plt +
  #     geom_point(data = subset(pdata, Filter %in% c("PASS", "ALL")),
  #                mapping=aes(group=igroup),
  #                size=8)
  
  xl_min = max(0,
               min(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) - 0.02)
  xl_max = min(1.0,
               max(subset(pdata, Filter %in% c("PASS", "ALL"))[[xaxis]]) + 0.02)
  yl_min = max(0,
               min(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) - 0.01)
  yl_max = min(1.0,
               max(subset(pdata, Filter %in% c("PASS", "ALL"))[[yaxis]]) + 0.01)
  
  plt = plt +
    scale_color_brewer("", palette="Set1") +
    xlim(c(xl_min, xl_max)) +
    ylim(c(yl_min, yl_max)) +
    theme_bw(base_size=18)
  
  return(plt)
}

# Update METRIC.F1_Score for a happy roc DataFrames
# We use this function for DataFrames gone through special processing (e.g. Recall and Precision from different DataFrames)
update_F1 <- function(df){
  df$METRIC.F1_Score = (2 / ((1 / df$METRIC.Recall) + (1 / df$METRIC.Precision) ))
  return(df)
}

plot_af_hist <- function(df, type){
  if (type %in% c('SNP', 'INDEL')){
    df <- df %>% filter(TYPE == type)
  }
  else if (type != 'ALL'){
    print('Incorrect type')
    return()
  }
  hs <- hist(df$AF, breaks=50, plot=F)
  total_counts <- sum(hs$counts)
  # print(total_counts)
  p <- ggplot() + geom_bar(mapping=aes(x=hs$mids, y=hs$density), stat = "identity", fill="lightblue", color="white") +
    scale_y_continuous(sec.axis = sec_axis(~(. * total_counts/100), name = "Count")) +
    theme_bw() +
    xlab("Population allele frequency") + ylab("Fraction (%)")
  return(p)
}

plot_af_hist_pair <- function(df, type){
  if (type %in% c('SNP', 'INDEL')){
    df <- df %>% filter(TYPE == type)
  }
  else if (type != 'ALL'){
    print('Incorrect type')
    return()
  }
  p <- ggplot(df, aes(AF, fill=label)) +
    geom_histogram(position="dodge", bins=50) +
    theme_bw() +
    xlab("Population allele frequency") + ylab("Count") + labs(fill="") +
    theme(legend.position="bottom") +
    # scale_fill_manual(values=c("#00AFBB", "#E7B800"))
    scale_fill_brewer(palette="Paired")
  return(p)
}
