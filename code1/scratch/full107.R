#!/usr/bin/env Rscript
# =============================================================================
# Univariable Mendelian Randomization Analysis
# Exposure: Gene Expression (GTEx eQTL data)
# Outcome: Alzheimer's Disease (GWAS data)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 加载必要的R包
# -----------------------------------------------------------------------------
library(data.table)      # 高效的数据读取和处理
library(tidyverse)       # 数据清洗和转换
library(pbapply)         # 进度条显示
library(MendelianRandomization)  # MR分析的标准包
library(mvtnorm)         # 多元正态分布
library(LaplacesDemon)   # 贝叶斯统计
library(dirmult)         # Dirichlet-multinomial分布
library(invgamma)        # 逆Gamma分布

# 加载自定义函数和C++编译的Gibbs采样器
source("/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code/init_setup_seso.R")
Rcpp::sourceCpp("/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code/gibbs_seso_uhp_only.cpp", 
                cacheDir = "/gpfs/data/linchen-lab/Bowei/")

# -----------------------------------------------------------------------------
# 2. LD Clumping函数
# 目的：选择独立的遗传工具变量，去除连锁不平衡(LD)造成的冗余SNPs
# -----------------------------------------------------------------------------
clump <- function(dat, 
                  SNP_col = "rsid",           # SNP ID列名
                  pval_col = "p",             # p值列名
                  clump_kb = 250,             # LD窗口大小(kb)
                  clump_r2 = 0.1,             # LD r²阈值
                  clump_p = 0.999,            # p值阈值
                  bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR",  # 参考基因型数据
                  plink_bin = genetics.binaRies::get_plink_binary(),  # PLINK路径
                  pop = "EUR") {              # 人群
  
  # 准备clumping所需的数据格式
  df <- data.frame(rsid = dat[, ..SNP_col], 
                   pval = dat[, ..pval_col])
  colnames(df) <- c("rsid", "pval")
  
  # 使用ieugwasr包进行LD clumping
  out <- tryCatch({
    ieugwasr::ld_clump(df, 
                       clump_kb = clump_kb, 
                       clump_r2 = clump_r2, 
                       clump_p = clump_p, 
                       bfile = bfile, 
                       plink_bin = plink_bin, 
                       pop = pop)
  }, silent = TRUE, error = function(x) return(NA))
  
  # 如果clumping失败，返回NA
  if (length(out) == 1) {
    return(NA)
  }
  
  # 返回clumping后的数据
  MRdat <- dat[which(unlist(dat[, ..SNP_col]) %in% out$rsid), ]
  return(MRdat)
}

# -----------------------------------------------------------------------------
# 3. 运行FusioMR-s方法的函数
# 这是一个基于贝叶斯框架的MR方法，使用Gibbs采样估计因果效应
# -----------------------------------------------------------------------------
run_mr_seso <- function(gene, 
                        b_exp,              # exposure的效应量(beta)
                        se_exp,             # exposure的标准误
                        b_out,              # outcome的效应量
                        se_out,             # outcome的标准误
                        a_gamma_prior,      # gamma先验参数a
                        b_gamma_prior,      # gamma先验参数b
                        a_theta_prior,      # theta先验参数a
                        b_theta_prior) {    # theta先验参数b
  
  # 设置MCMC迭代次数和工具变量数量
  niter <- 20000
  K <- length(b_exp)
  
  # 初始化MCMC参数
  start_val <- init_setup(niter, K, 
                          alpha_init = 1, 
                          beta_init = 0, 
                          sigma_gamma_init = 0.1, 
                          sigma_theta_init = 0.1)
  
  # 使用C++编译的Gibbs采样器进行贝叶斯估计
  res <- gibbs_seso_uhp_only_cpp(
    niter, K, 
    start_val$beta_tk, 
    start_val$theta_tk, 
    start_val$gamma_tk, 
    start_val$sigma2_gamma_tk, 
    start_val$sigma2_theta_tk, 
    b_out, b_exp, 
    (se_out)^2, (se_exp)^2, 
    a_gamma_prior, b_gamma_prior, 
    a_theta_prior, b_theta_prior
  )
  
  # 计算后验估计（去除burn-in期，使用后50%的迭代结果）
  if (all(!is.na(res))) {
    ids <- (niter/2 + 1):niter
    bhat <- mean(res$beta_tk[ids], na.rm = TRUE)
    se_bhat <- sd(res$beta_tk[ids], na.rm = TRUE)
    pval <- 2 * (1 - pnorm(abs(bhat) / se_bhat))
  } else { 
    bhat <- se_bhat <- pval <- NA 
  }
  
  # 整理输出结果
  out_fusio <- c(bhat, se_bhat, pval)
  output <- c(gene, K, out_fusio)
  names(output) <- c("gene", "niv", "b_fusio", "se_fusio", "p_fusio")
  
  return(output)
}

# -----------------------------------------------------------------------------
# 4. 处理单个基因的MR分析
# 这是核心函数，整合了数据筛选、clumping和MR分析
# -----------------------------------------------------------------------------
process_gene_bk_only <- function(dat1,                  # 输入数据
                                 gene_list1,            # 基因列表
                                 j,                     # 当前基因索引
                                 pcut = 0.01,           # eQTL p值阈值
                                 clump_kb = 10,         # clumping窗口
                                 clump_r2 = 0.1,        # clumping r²阈值
                                 b_gamma_prior_mean,    # gamma先验均值
                                 b_theta_prior_mean) {  # theta先验均值
  
  message(paste0("Processing gene-", j, ": ", gene_list1[j]))
  
  # 步骤1：筛选当前基因的显著eQTL (p < pcut)
  tmp <- dat1 %>% 
    dplyr::filter(gene_name == gene_list1[j]) %>% 
    dplyr::filter(p <= pcut)
  
  if (nrow(tmp) > 0) {
    # 步骤2：基于数据重新估计超参数
    b_gamma_prior_mean <- max(1e-3, mean(tmp$beta^2) - mean(tmp$se^2))
    b_theta_prior_mean <- max(1e-5, mean(tmp$b_out^2) - mean(tmp$se_out^2))
    
    # 步骤3：LD clumping，选择独立的SNPs
    clumped_df <- clump(tmp, 
                        SNP_col = "rsid", 
                        clump_kb = clump_kb, 
                        clump_r2 = clump_r2, 
                        bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", 
                        pval_col = "p", 
                        plink_bin = genetics.binaRies::get_plink_binary())
    
    # 步骤4：检查clumping结果是否满足MR分析条件
    # 需要至少2个SNPs，且不超过1000个（避免过度拟合）
    if (all(class(clumped_df) == c("data.table", "data.frame")) &&
        nrow(clumped_df) > 1 && nrow(clumped_df) < 1000) {
      
      # 提取效应量和标准误
      b_exp <- clumped_df$beta
      se_exp <- clumped_df$se
      b_out <- clumped_df$b_out
      se_out <- clumped_df$se_out
      
      # 步骤5：设置贝叶斯先验参数
      K <- nrow(clumped_df)
      a_gamma_prior <- a_theta_prior <- max(2, K/4)
      b_gamma_prior <- b_gamma_prior_mean * (a_gamma_prior - 1)
      b_theta_prior <- b_theta_prior_mean * (a_theta_prior - 1)
      
      # 步骤6：运行MR分析
      res <- run_mr_seso(gene_list1[j], 
                         b_exp, se_exp, 
                         b_out, se_out, 
                         a_gamma_prior, b_gamma_prior, 
                         a_theta_prior, b_theta_prior)
      
      # 添加超参数信息到结果中
      res <- c(res, 
               'b_gamma_prior_mean' = b_gamma_prior_mean, 
               'b_theta_prior_mean' = b_theta_prior_mean)
      
      return(res)
      
    } else {
      # clumping后SNP数量不满足条件
      return(c(gene_list1[j], -1, rep(NA, 5)))
    }
    
  } else {
    # 没有显著的eQTL
    return(c(gene_list1[j], -2, rep(NA, 5)))
  }
}

# -----------------------------------------------------------------------------
# 5. 主分析流程
# 对每个染色体、每个脑组织进行MR分析
# -----------------------------------------------------------------------------

# 设置分析参数
pcut <- 0.001          # eQTL显著性阈值
clump_kb <- 50         # LD窗口50kb
clump_r2 <- 0.1        # LD r²阈值0.1

# 定义脑组织列表
tissue_short <- c("Amygdala", "BA24", "Caudate_basal_ganglia", 
                  "Cerebellar_Hemi", "Cerebellum", "Cortex", 
                  "BA9", "Hippocampus", "Hypothalamus", 
                  "Nucleus_accumbens_basal_ganglia", 
                  "Putamen_basal_ganglia", "Spinal_cord", 
                  "Substantia_nigra")

# 选择要分析的组织（这里选择第6个：Cortex）
j <- 6

# 循环分析每条染色体
for (chr in 1:22) {
  
  print(paste0("Processing chromosome ", chr))
  
  # 步骤1：读取harmonized数据（已经匹配好的eQTL和GWAS数据）
  dat <- fread(paste0('harmo_bk/eqtl_ad_', tissue_short[j], '_', chr, '.txt.gz'))
  
  # 步骤2：数据清洗
  dat1 <- dat %>% 
    dplyr::rename(beta = b) %>%                    # 重命名列
    dplyr::filter(se != Inf, se_out != Inf,        # 去除无穷值
                  se > 0, se_out > 0) %>%          # 去除非正值
    drop_na()                                       # 去除缺失值
  
  # 步骤3：获取基因列表
  gene_list1 <- unique(dat1$gene_name)
  message(paste0("Number of genes: ", length(gene_list1)))
  
  # 步骤4：估计全局超参数（基于所有数据）
  b_gamma_prior_mean <- max(1e-3, mean(dat1$beta^2) - mean(dat1$se^2))
  b_theta_prior_mean <- max(1e-5, mean(dat1$b_out^2) - mean(dat1$se_out^2))
  
  # 步骤5：对每个基因并行运行MR分析
  sumtbl1 <- do.call(rbind, pblapply(
    seq_along(gene_list1), 
    function(t) process_gene_bk_only(
      dat1, gene_list1, t, 
      pcut, clump_kb, clump_r2, 
      b_gamma_prior_mean, b_theta_prior_mean
    )
  ))
  
  # 步骤6：保存结果
  output_file <- paste0("sumtbl_bkonly/bkonly_ad_", 
                        pcut, "_", clump_kb, "_", clump_r2, "_", 
                        tissue_short[j], "_", chr, '.txt')
  fwrite(sumtbl1, output_file)
  
  message(paste0("Results saved to: ", output_file))
}

message("Analysis completed!")

# -----------------------------------------------------------------------------
# 6. 结果解读说明
# -----------------------------------------------------------------------------
# 输出文件包含以下列：
# - gene: 基因名称
# - niv: 工具变量数量（独立SNPs数）
#   * -2: 该基因没有显著的eQTL (p < pcut)
#   * -1: clumping后SNP数量 ≤1 或 >1000
#   * >0: 成功进行MR分析的SNP数量
# - b_fusio: FusioMR估计的因果效应
# - se_fusio: 因果效应的标准误
# - p_fusio: 因果效应的p值
# - b_gamma_prior_mean: gamma先验的均值参数
# - b_theta_prior_mean: theta先验的均值参数
# -----------------------------------------------------------------------------