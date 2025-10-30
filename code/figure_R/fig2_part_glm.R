rm(list = ls())
#setwd("E:/sunyingying/1111LS-GNSHF-project/1111data/01.IGC_humanswiss/01.GNHSF/09.figure/Figure4_permanova_maaslin/glm_20250828")

if (!require(pscl, quietly = TRUE)) {
  install.packages("pscl")}
library(pscl)## need to install first to calculate peudo R2 of GLM model
library("parallel")
library(pbapply)

######### generate 1385 matrix
cinfo <- read.csv("GNHSF_sample_inform_normalized_1385_2_yes_all_metadata.csv")
type <- "species"
func.1385 <- function(type){
  df <- read.delim2(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample.tsv"))
  colnames(df)[1] <- "sample"
  df.1385 <- df[, colnames(df) %in% cinfo$sample]
  rownames(df.1385) <- df$sample
  df.1385 <- df.1385[!apply(df.1385, 1, function(x) all(is.na(x))), ]
  write.table(df.1385, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_rmallNA.tsv"), sep = "\t", col.names=NA)
}

func.1385("cog")
func.1385("kegg")
func.1385("species")
func.1385("genus")
func.1385("phylum")
func.1385("microprotein")
func.1385("humanprotein")


### impute 0 and log2 matrix
func.abun <- function(type){
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_rmallNA.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  df[is.na(df)] <- 0
  write.table(df, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_rmallNA_0.tsv"), sep = "\t", row.names = F)
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_rmallNA.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  df[is.na(df)] <- 1
  dflog <- log2(df[,2:ncol(df)])
  dflog <- cbind(df[,1], dflog)
  colnames(dflog)[1] <- colnames(df)[1]
  write.table(dflog, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_rmallNA_log2.tsv"), sep = "\t", row.names = F)
}
func.abun("cog")
func.abun("kegg")
func.abun("species")
func.abun("genus")
func.abun("phylum")
func.abun("microprotein")
func.abun("humanprotein")
## NA cutoff 0.9 and impute 0

func.prev <- function(x){
  sum(is.na(x))
}

func.NA <- function(type){
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_rmallNA.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  df$prevalence <- ((ncol(df)-1) - as.numeric(apply(df[, 2:ncol(df)], 1, func.prev))) / (ncol(df)-1)
  df_NA90 <- df[df$prevalence > 0.1,]
  write.table(df_NA90, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_prevelance.tsv"), sep = "\t", row.names = F)
  df_NA90 <- df_NA90[,-ncol(df_NA90)]
  df_NA90[is.na(df_NA90)] <- 0
  write.table(df_NA90, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90.tsv"), sep = "\t", row.names = F)
}
func.NA("cog")
func.NA("kegg")
func.NA("species")
func.NA("genus")
func.NA("phylum")
func.NA("microprotein")
func.NA("humanprotein")

## INT transform
func.INT <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}
type <- "genus"
func.INT.abun <- function(type){
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  df.int <- as.data.frame(t(apply(df[,2:ncol(df)], 1, func.INT)), check.name = F)
  df.int <- cbind(df[,1], df.int)
  colnames(df.int)[1] <- colnames(df)[1]
  write.table(df.int, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90_INT.tsv"), sep = "\t", row.names = F)
}
func.INT.abun("cog")
func.INT.abun("kegg")
func.INT.abun("species")
func.INT.abun("genus")
func.INT.abun("phylum")
func.INT.abun("microprotein")
func.INT.abun("humanprotein")

######################## GLM all features ###########################
sampledf <- as.data.frame(read.csv(file = "GNHSF_sample_inform_normalized_1385_2_yes_all_metadata.csv", header = TRUE, row.names = 1, stringsAsFactors = TRUE))
cat <- sort.int(unique(sampledf$edu3))
sampledf$edu3 = factor(sampledf$edu3, levels = cat)
cat <- sort.int(unique(sampledf$income4))
sampledf$income4 = factor(sampledf$income4, levels = cat)
sampledf$bristol_scale <- gsub("Normal", "0", sampledf$bristol_scale)
sampledf$bristol_scale <- gsub("Diarrhea", "1", sampledf$bristol_scale)
sampledf$bristol_scale <- gsub("Constipation", "2", sampledf$bristol_scale)
cat <- sort.int(unique(sampledf$bristol_scale))
sampledf$bristol_scale = factor(sampledf$bristol_scale, levels = cat)
sampledf$marrige4 <- gsub("unmarried", "3", sampledf$marrige4)
sampledf$marrige4 <- gsub("married", "0", sampledf$marrige4)
sampledf$marrige4 <- gsub("divorced", "1", sampledf$marrige4)
sampledf$marrige4 <- gsub("widowed", "2", sampledf$marrige4)

cat <- sort.int(unique(sampledf$marrige4))
sampledf$marrige4 = factor(sampledf$marrige4, levels = cat)

summary(sampledf)
sampledf_tmp <- sampledf
sampledf_tmp$sample <- rownames(sampledf)
sambf <- as.data.frame(read.csv(file = "GNHSF_sample_ident_20220307.csv"))
sambf.1385 <- sambf[sambf$label3 %in% rownames(sampledf), c(4, 7)]
sampledf_tmp <- merge(sampledf_tmp, sambf.1385, by.x = "sample", by.y = "label3")
sampledf_tmp <- sampledf_tmp[, c(1, ncol(sampledf_tmp), 2:(ncol(sampledf_tmp)-1))]

sampledf <- sampledf_tmp[, c(2:ncol(sampledf_tmp))]
rownames(sampledf) <- sampledf_tmp$sample
rownames(sampledf_tmp) <- sampledf_tmp$sample

sampledf$batch <- as.factor(sampledf$batch)
sampledf_tmp$batch <- as.factor(sampledf_tmp$batch)

colnames(sampledf_tmp)
colnames(sampledf)

table(sampledf$income4)


func.glm.common <- function(type){
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90_INT.tsv"), header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE))
  dft <- as.data.frame(t(df))
  colnames(dft) <- as.character(lapply(strsplit(as.character(colnames(dft)), split=" "), head, n=1))
  if (type == "microprotein"){
    colnames(dft) <- paste0("X", colnames(dft))
  }
  colnames(dft) <- gsub("\\|", ".", colnames(dft))
  colnames(dft) <- gsub("\\-", ".", colnames(dft))
  colnames(dft) <- gsub("\\.$", "", colnames(dft))
  colnames(dft)
  dft$sample <- rownames(dft)
  glm.res.final <- data.frame()
  for (i in 1:(ncol(dft)-1)){
    print(i)
    df.glm <- merge(dft[,c(i, ncol(dft))], sampledf_tmp, by.x = "sample", by.y = "sample")
    for (j in 2:ncol(sampledf)){
      print(paste0(i, "\t", j))
      fomular <- ""
      if ((colnames(sampledf)[j] == "age") || (colnames(sampledf)[j] == "sex") || (colnames(sampledf)[j] == "bristol_scale")){
        fomular <- paste0(colnames(df.glm)[2], " ~ sex + age + bristol_scale")
      }else{
        fomular <- paste0(colnames(df.glm)[2], " ~ sex + age + bristol_scale + ", colnames(sampledf)[j])
      }
      glm.run <- glm(fomular, data = df.glm, family = gaussian)
      glm.sum <- summary(glm.run)
      glm.res <- as.data.frame(glm.sum$coefficients)
      glm.res$feature <- rep(colnames(df.glm)[2], nrow(glm.res))
      glm.res$metadata <- rownames(glm.res)
      glm.res$AIC <- glm.sum$aic
      glm.res$BIC <- BIC(glm.run)
      r2 <- pR2(glm.run)
      glm.res$mcfadden.R2 <- r2["McFadden"]
      glm.res$r2ml.R2 <- r2["r2ML"]

      glm.res.row <- c()
      glm.res.add <- data.frame()
      if (length(c(grep(colnames(sampledf)[j], rownames(glm.res)))) > 0){
        glm.res.row <- c(grep(colnames(sampledf)[j], rownames(glm.res)))
        glm.res.add <- glm.res[glm.res.row,]
        glm.res.final <- rbind(glm.res.final, glm.res.add)
      }
    }
  }
  colnames(glm.res.final)
  glm.res.final$q <- p.adjust(glm.res.final$`Pr(>|t|)`, method = "BH", n = length(glm.res.final$`Pr(>|t|)`))
  q.sig <- c()
  for (i in 1:nrow(glm.res.final)) {
    if (glm.res.final$q[i] < 0.001){
      q.sig[i] <- "***"
    }else if (glm.res.final$q[i] < 0.01){
      q.sig[i] <- "**"
    }else if (glm.res.final$q[i] < 0.05){
      q.sig[i] <- "*"
    }else{
      q.sig[i] <- "nsig"
    }
  }
  glm.res.final$q.sig <- q.sig
  write.csv(glm.res.final, file = paste0("glm/", type, "_NA90_all_meta_glm_res.csv"))
}


func.glm.common.parallel <- function(type, n_cores = NULL) {
  closeAllConnections()
  # 确定使用的核心数
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1  # 留一个核心给系统
  }
  cat("使用", n_cores, "个核心进行并行计算\n")
  
  # 数据预处理（和原来一样）
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90_INT.tsv"), 
                                 header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE))
  dft <- as.data.frame(t(df))
  colnames(dft) <- as.character(lapply(strsplit(as.character(colnames(dft)), split=" "), head, n=1))
  
  if (type == "microprotein") {
    colnames(dft) <- paste0("X", colnames(dft))
  }
  colnames(dft) <- gsub("\\|", ".", colnames(dft))
  colnames(dft) <- gsub("\\-", ".", colnames(dft))
  colnames(dft) <- gsub("\\.$", "", colnames(dft))
  dft$sample <- rownames(dft)
  
  # 创建任务列表（i,j 的所有组合）
  tasks <- expand.grid(i = 1:(ncol(dft)-1), j = 2:ncol(sampledf))
  cat("总共", nrow(tasks), "个GLM任务\n")
  
  # 单个GLM任务的函数
  run_single_glm <- function(task_row, dft, sampledf, sampledf_tmp) {
    i <- task_row$i
    j <- task_row$j
    
    # 合并数据
    df.glm <- merge(dft[,c(i, ncol(dft))], sampledf_tmp, by.x = "sample", by.y = "sample")
    
    # 构建公式
    if ((colnames(sampledf)[j] %in% c("age", "sex", "bristol_scale"))) {
      fomular <- paste0(colnames(df.glm)[2], " ~ sex + age + bristol_scale")
    } else {
      fomular <- paste0(colnames(df.glm)[2], " ~ sex + age + bristol_scale + ", colnames(sampledf)[j])
    }
    
    # 运行GLM
    tryCatch({
      glm.run <- glm(fomular, data = df.glm, family = gaussian)
      glm.sum <- summary(glm.run)
      glm.res <- as.data.frame(glm.sum$coefficients)
      glm.res$feature <- rep(colnames(df.glm)[2], nrow(glm.res))
      glm.res$metadata <- rownames(glm.res)
      glm.res$AIC <- glm.sum$aic
      glm.res$BIC <- BIC(glm.run)
      r2 <- pR2(glm.run)
      glm.res$mcfadden.R2 <- r2["McFadden"]
      glm.res$r2ml.R2 <- r2["r2ML"]
      
      # 筛选相关行
      target_rows <- grep(colnames(sampledf)[j], rownames(glm.res))
      if (length(target_rows) > 0) {
        return(glm.res[target_rows, ])
      } else {
        return(NULL)
      }
    }, error = function(e) {
      cat("Error in task i=", i, " j=", j, ": ", e$message, "\n")
      return(NULL)
    })
  }
  
  # 并行执行并显示进度条
  cat("开始并行计算并显示进度条...\n")
  start_time <- Sys.time()
  
  # 使用pbapply进行并行计算并显示进度条
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl)) 
  # 将必要的对象和包导出到工作节点
  clusterExport(cl, c("dft", "sampledf", "sampledf_tmp", "pR2", "BIC", "run_single_glm"), envir=environment())
  clusterEvalQ(cl, {
    library(pscl)
  })
  
  # 使用pbapply::pblapply替代parApply，显示进度条
  results_list <- pblapply(1:nrow(tasks), function(idx) run_single_glm(tasks[idx, ], dft, sampledf, sampledf_tmp), cl = cl)
  
  end_time <- Sys.time()
  cat("并行计算完成，用时:", difftime(end_time, start_time, units = "mins"), "分钟\n")
  
  # 合并结果
  cat("合并结果...\n")
  glm.res.final <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  
  # 统一计算调整后的p值
  cat("计算调整后的p值...\n")
  glm.res.final$q <- p.adjust(glm.res.final$`Pr(>|t|)`, method = "BH", n = length(glm.res.final$`Pr(>|t|)`))
  
  # 添加显著性标记
  q.sig <- ifelse(glm.res.final$q < 0.001, "***",
                  ifelse(glm.res.final$q < 0.01, "**",
                         ifelse(glm.res.final$q < 0.05, "*", "nsig")))
  glm.res.final$q.sig <- q.sig
  
  # 保存结果
  output_file <- paste0("glm/", type, "_NA90_all_meta_glm_res.csv")
  write.csv(glm.res.final, file = output_file, row.names = FALSE)
  
  cat("分析完成！结果保存至:", output_file, "\n")
  cat("总共分析了", nrow(glm.res.final), "个有效结果\n")
  
  return(glm.res.final)
}



func.glm.common.parallel.addbatch <- function(type, n_cores = NULL) {
  closeAllConnections()
   # 确定使用的核心数
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1  # 留一个核心给系统
  }
  cat("使用", n_cores, "个核心进行并行计算\n")
  
  # 数据预处理（和原来一样）
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90_INT.tsv"), 
                                 header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE))
  dft <- as.data.frame(t(df))
  colnames(dft) <- as.character(lapply(strsplit(as.character(colnames(dft)), split=" "), head, n=1))
  
  if (type == "microprotein") {
    colnames(dft) <- paste0("X", colnames(dft))
  }
  colnames(dft) <- gsub("\\|", ".", colnames(dft))
  colnames(dft) <- gsub("\\-", ".", colnames(dft))
  colnames(dft) <- gsub("\\.$", "", colnames(dft))
  dft$sample <- rownames(dft)
  
  # 创建任务列表（i,j 的所有组合）
  tasks <- expand.grid(i = 1:(ncol(dft)-1), j = 2:ncol(sampledf))
  cat("总共", nrow(tasks), "个GLM任务\n")
  
  # 单个GLM任务的函数
  run_single_glm.batch <- function(task_row, dft, sampledf, sampledf_tmp) {
    i <- task_row$i
    j <- task_row$j
    
    # 合并数据
    df.glm <- merge(dft[,c(i, ncol(dft))], sampledf_tmp, by.x = "sample", by.y = "sample")
    
    # 构建公式
    if ((colnames(sampledf)[j] %in% c("age", "sex", "bristol_scale"))) {
      fomular <- paste0(colnames(df.glm)[2], " ~ sex + age + bristol_scale + batch")
    } else {
      fomular <- paste0(colnames(df.glm)[2], " ~ sex + age + bristol_scale + batch + ", colnames(sampledf)[j])
    }
    
    # 运行GLM
    tryCatch({
      glm.run <- glm(fomular, data = df.glm, family = gaussian)
      glm.sum <- summary(glm.run)
      glm.res <- as.data.frame(glm.sum$coefficients)
      glm.res$feature <- rep(colnames(df.glm)[2], nrow(glm.res))
      glm.res$metadata <- rownames(glm.res)
      glm.res$AIC <- glm.sum$aic
      glm.res$BIC <- BIC(glm.run)
      r2 <- pR2(glm.run)
      glm.res$mcfadden.R2 <- r2["McFadden"]
      glm.res$r2ml.R2 <- r2["r2ML"]
      
      # 筛选相关行
      target_rows <- grep(colnames(sampledf)[j], rownames(glm.res))
      if (length(target_rows) > 0) {
        return(glm.res[target_rows, ])
      } else {
        return(NULL)
      }
    }, error = function(e) {
      cat("Error in task i=", i, " j=", j, ": ", e$message, "\n")
      return(NULL)
    })
  }
  
  # 并行执行并显示进度条
  cat("开始并行计算并显示进度条...\n")
  start_time <- Sys.time()
  
  # 使用pbapply进行并行计算并显示进度条
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl)) 
  # 将必要的对象和包导出到工作节点
  clusterExport(cl, c("dft", "sampledf", "sampledf_tmp", "pR2", "BIC", "run_single_glm.batch"), envir=environment())
  clusterEvalQ(cl, {
    library(pscl)
  })
  
  # 使用pbapply::pblapply替代parApply，显示进度条
  results_list <- pblapply(1:nrow(tasks), function(idx) run_single_glm.batch(tasks[idx, ], dft, sampledf, sampledf_tmp), cl = cl)  

  
  end_time <- Sys.time()
  cat("并行计算完成，用时:", difftime(end_time, start_time, units = "mins"), "分钟\n")
  
  # 合并结果
  cat("合并结果...\n")
  glm.res.final <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  
  # 统一计算调整后的p值
  cat("计算调整后的p值...\n")
  glm.res.final$q <- p.adjust(glm.res.final$`Pr(>|t|)`, method = "BH", n = length(glm.res.final$`Pr(>|t|)`))
  
  # 添加显著性标记
  q.sig <- ifelse(glm.res.final$q < 0.001, "***",
                  ifelse(glm.res.final$q < 0.01, "**",
                         ifelse(glm.res.final$q < 0.05, "*", "nsig")))
  glm.res.final$q.sig <- q.sig
  
  # 保存结果
  output_file <- paste0("glm.batch/", type, "_NA90_all_meta_glm_res.csv")
  write.csv(glm.res.final, file = output_file, row.names = FALSE)
  
  cat("分析完成！结果保存至:", output_file, "\n")
  cat("总共分析了", nrow(glm.res.final), "个有效结果\n")  
  
  return(glm.res.final)
}

dir.create("glm")
dir.create("glm.batch")

func.glm.common.parallel("microprotein", n_cores = NULL)
func.glm.common.parallel.addbatch("microprotein", n_cores = NULL)


func.glm.common("cog")
func.glm.common("kegg")
func.glm.common("phylum")
func.glm.common("genus")
func.glm.common("species")
func.glm.common("humanprotein")
func.glm.common("microprotein")

############ sep glm results ###############
process_glm_files <- function(type, cinfo) {
  # 构造文件路径
  file_path <- paste0("glm/", type, "_NA90_all_meta_glm_res.csv")
  
  if (!file.exists(file_path)) {
    warning("File not found:", file_path)
    return()
  }
  
  # 读取数据
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  header <- colnames(data)
  new_header <- c(header[2:7], "metadata_value", header[8:length(header)])
  
  # 初始化 processed_data
  processed_data <- data.frame()
  
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    
    if (row[1] == "" || is.na(row[1])) next
    
    metadata <- as.character(row[7])
    
    # 匹配 metadata
    match_result <- sapply(names(cinfo), function(k) {
      if (grepl(paste0("^", k), metadata)) {
        return(list(metadatanew = k, metavalue = sub(paste0("^", k), "", metadata)))
      } else {
        return(NULL)
      }
    })
    
    match_result <- match_result[!sapply(match_result, is.null)]
    
    if (length(match_result) > 0) {
      metadatanew <- as.character(match_result[[1]]$metadatanew)
      metavalue   <- as.character(match_result[[1]]$metavalue)
    } else {
      metadatanew <- ""
      metavalue   <- ""
    }
    
    # ---- 核心改动：直接构建向量，再转data.frame ----
    if (ncol(data) >= 8) {
      new_row <- c(row[2:6], metadatanew, metavalue, row[8:ncol(data)])
    } else {
      new_row <- c(row[2:min(6, ncol(data))], metadatanew, metavalue)
      if (length(new_row) < length(new_header)) {
        new_row <- c(new_row, rep("", length(new_header) - length(new_row)))
      }
    }
    
    # 转为一行 data.frame，并加上列名
    new_row <- as.data.frame(as.list(new_row), stringsAsFactors = FALSE)
    colnames(new_row) <- new_header
    
    # 追加
    processed_data <- rbind(processed_data, new_row)
  }
  
  # 写出结果
  output_file <- paste0("glm/", sub("\\.csv$", "", basename(file_path)), "_spmeta.csv")
  write.csv(processed_data, output_file, row.names = FALSE, quote = TRUE)
  
  cat("Processed:", file_path, "->", output_file, "\n")
}


# 主函数
main_process <- function() {
  # 设置metadata文件路径
  cinfo_file <- "metadata_classify_add4.csv"
  
  if (!file.exists(cinfo_file)) {
    stop("Metadata classification file not found:", cinfo_file)
  }
  
  # 读取metadata分类文件
  cinfo_data <- read.csv(cinfo_file, stringsAsFactors = FALSE)
  cinfo <- setNames(rep(1, nrow(cinfo_data)), cinfo_data[, 1])
  
  # 要处理的类型列表
  types <- c("genus", "species")  # 根据需要修改
  
  # 处理每个类型
  for (type in types) {
    process_glm_files(type, cinfo)
  }
  
  cat("All processing complete!\n")
}

# 调用主函数
main_process()


########### combine res ###################

metatype <- read.csv(file = "glm/metadata_classify_add4.csv")
func.qsig <- function(x){
  q.sig <- c()
  for (i in 1:nrow(x)) {
    if (x$q[i] < 0.001){
      q.sig[i] <- "***"
    }else if (x$q[i] < 0.01){
      q.sig[i] <- "**"
    }else if (x$q[i] < 0.05){
      q.sig[i] <- "*"
    }else{
      q.sig[i] <- "nsig"
    }
  }
  return(q.sig)
}

type <- "genus"

func.adj <- function(type){
  df <- read.csv(file = paste0("glm/", type, "_NA90_all_meta_glm_res_spmeta.csv"))
  ## adjustp by meta except tg_cl/hdl_cl/ldl_cl/tc_cl/dm_med/dys_med/hyper_med
  df_excpt <- df[(df$metadata != "tg_cl" & df$metadata != "hdl_cl" & df$metadata != "ldl_cl" & df$metadata != "tc_cl"),]
  df_excpt$q <- p.adjust(df_excpt$`Pr...t..`, method = "BH", n = length(df_excpt$`Pr...t..`))
  df_excpt$q.sig <- func.qsig(df_excpt)
  write.csv(df_excpt, file = paste0("glm/", type, "_NA90_all_meta_glm_res_adj.csv"), row.names = F)
  ## adjustp by category except tg_cl/hdl_cl/ldl_cl/tc_cl/dm_med/dys_med/hyper_med
  df_excpt_type <- merge(df_excpt, metatype, by.x = "metadata", by.y = "Cinfo")
  cinfo_type <- unique(metatype$Cinfo_type)
  cinfo_type
  
  p.df <- data.frame()
  allmeta.p <- 0
  allmeta.padj <- 0
  allmeta.meta <- 0
  for (i in 1:length(cinfo_type)){
    cinfonow <- cinfo_type[i]
    df_excpt_typenow <- df_excpt_type[df_excpt_type$Cinfo_type == cinfonow,]
    df_excpt_typenow$q <- p.adjust(df_excpt_typenow$`Pr...t..`, method = "BH", n = length(df_excpt_typenow$`Pr...t..`))
    df_excpt_typenow$q.sig <- func.qsig(df_excpt_typenow)
    write.csv(df_excpt_typenow, file = paste0("glm/", type, "_NA90_", cinfonow, "_glm_res_adj_category.csv"), row.names = F)
    adj <- 0
    p <- 0
    metadata.adj <- c()
    for (i in 1:nrow(df_excpt_typenow)){
      if (df_excpt_typenow[i, 5] < 0.05){
        p <- p + 1
      }
      if (df_excpt_typenow[i, 9] != "nsig"){
        adj <- adj + 1
        metadata.adj <- append(metadata.adj, df_excpt_typenow$metadata[i])
      }
    }
    allmeta.p <- allmeta.p + p
    allmeta.padj <- allmeta.padj + adj
    
    adjv <- c(type, "p.adj", adj, cinfonow)
    pv <- c(type, "p", p, cinfonow)
    metadata.adj <- unique(metadata.adj)
    allmeta.meta <- allmeta.meta + length(metadata.adj)
    metadata.adj.num <- c(type, "metadata", length(metadata.adj), cinfonow)
    p.df.add <- data.frame()
    p.df.add <- rbind(adjv, pv, metadata.adj.num)
    p.df <- rbind(p.df, p.df.add)
  }
  total <- c(type, allmeta.p, allmeta.padj, allmeta.meta)
  p.df <- rbind(p.df, total)
  write.csv(p.df, file = paste0("glm/", type, "_NA90_glm_res_adj_statistic.csv"))
}

func.adj("cog")
func.adj("kegg")
func.adj("phylum")
func.adj("genus")
func.adj("species")
func.adj("humanprotein")
func.adj("microprotein")

