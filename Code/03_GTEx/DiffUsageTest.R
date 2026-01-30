#!/usr/bin/env Rscript

# libraries
library(reshape2)
library(stringr)
library(dplyr)
library(EMT)
library(XNomial)
# data
counts<-read.table("output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv", h = T)
countsLong<-melt(counts, variable.name = "sample", value.name = "count")
names(countsLong)<-c("junction", "sample", "count")
# split junction original name to get new columns
countsLong[c("ID", "gene", "exon_group", "frame", "end")]<-str_split_fixed(countsLong$junction, pattern =  "_",n = 5) 
# create a pair_ID
countsLong$pairID<-paste(countsLong$ID, "_", countsLong$end, sep = "")
# add sample metadata
gtexMetadata<-read.table("data/GTExMetadata.tsv")
names(gtexMetadata)=c("sample", "tissue", "lib_size", "color", "ShortName")
done<-merge(countsLong, gtexMetadata, by = "sample")
# export table for DCExplorer
write.table(done, "output/03_GTEx/AllFrames.pairedCounts_long.tsv", quote = F, row.names = F)
# functions
# Main test loop function
DCE.nomial.tests <- function(input, raw_count_threshold = 2) {
  out_pair <- character()
  # Initialize outputs
  out_p <- numeric()
  out_frames <- numeric()
  out_maxdif <- numeric()
  out_pattern <- character()
  out_n <- numeric()
  out_N <- numeric()
  best_comparison <- character()
  
  # Loop through unique pair IDs
  for (pair_name in unique(input$pairID)) {
    test <- input %>% filter(pairID == pair_name)
    
    # Calculate mean counts for each tissue
    sReads <- test %>%
      group_by(tissue) %>%
      summarise(meanCount = mean(count), .groups = "drop")
    
    # Filter tissues with mean count above the threshold
    t_to_keep <- sReads %>% filter(meanCount > raw_count_threshold) %>% pull(tissue)
    
    if (length(t_to_keep) > 1) {
      test_f <- test %>% filter(tissue %in% t_to_keep)
      
      # Summarize counts and calculate ratios
      TF_Scount <- test_f %>% group_by(tissue, frame) %>% summarise(TF_Scount = sum(count), .groups = "drop")
      T_Scount <- test_f %>% group_by(tissue) %>% summarise(T_Scount = sum(count), .groups = "drop")
      Scount <- merge(TF_Scount, T_Scount, by = "tissue")
      Scount$ratio <- Scount$TF_Scount / Scount$T_Scount
      
      # Calculate max difference in ratios for each frame
      fr_max_diff <- Scount %>%
        group_by(frame) %>%
        summarise(max_diff = max(ratio) - min(ratio), .groups = "drop")
      
      if (max(fr_max_diff$max_diff) >= 0.4) {
        # Perform pairwise tests
        test_f_mat <- dcast(test_f, frame ~ tissue, value.var = "count", fun.aggregate = sum)[, -1]
        test_f_mat = test_f_mat + 1 # avoid p = 1 for expected ratio, which breaks the binomial test
        n <- ncol(test_f_mat)
        N <- choose(n, 2)  # Number of pairwise comparisons
        Z <- data.frame(Comparison = character(N), p = numeric(N), stringsAsFactors = FALSE)
        k <- 1
        
        for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
            Namea <- colnames(test_f_mat)[i]
            Nameb <- colnames(test_f_mat)[j]
            Dataz <- matrix(c(test_f_mat[, i], test_f_mat[, j]), ncol = 2, byrow = FALSE)
            if (nrow(Dataz) == 2) {
              # Binomial test for 2 frames
              test_p <- Dataz[1, 1] / sum(Dataz[, 1])  # Target probability
              test_x <- Dataz[1, 2]  # Observed success
              test_n <- sum(Dataz[, 2])  # Total trials
              test_pVal <- binom.test(test_x, test_n, test_p, alternative = "two.sided")$p.value
              Z$p[k] <- test_pVal * N  # Bonferroni correction
            } else {
              # Multinomial test for >2 frames
              DatazM <- Dataz + 1
              test_pVal <- xmulti2(obs = DatazM[, 1], expr = DatazM[, 2])$pLLR
              Z$p[k] <- test_pVal * N  # Bonferroni correction
            }
            Z$Comparison[k] <- paste(Namea, ":", Nameb)
            k <- k + 1
          }
        }
        
        # Store results
        out_pair <- c(out_pair, pair_name)
        out_n <- c(out_n, n)
        out_N <- c(out_N, N)
        out_p <- c(out_p, min(Z$p))
        mP<-min(Z$p)
        best_comparison <- c(best_comparison, unique(Z[Z$p == mP,]$Comparison)[1])
        #print(Z[Z$p == mP,]$Comparison)
        out_maxdif <- c(out_maxdif, max(fr_max_diff$max_diff))
        out_frames <- c(out_frames, length(unique(test_f$frame)))
        out_pattern <- c(out_pattern, pattern_test(data = test_f, metric = "counts", threshold = 0.6))
        
        # find the best pair
        bestPs<-Z[Z$p == mP,]
        if(length(bestPs$Comparison) > 1){
          bestPtiss<-unique(unlist(str_split(bestPs$Comparison, " : ")))
          summary_data <- test_f %>%
            group_by(tissue, frame) %>%
            summarise(meanVal = median(count), .groups = "drop") %>%
            group_by(tissue) %>%
            mutate(total = sum(meanVal), ratio = meanVal / total)
          # find largest ratio per frame among tissues with the strongest P-values
          filtered_ratios<-summary_data[summary_data$tissue %in% bestPtiss,] %>% group_by(frame) # %>%  slice_max(ratio, n = 1, with_ties = FALSE)
          filtered_ratios<-filtered_ratios[filtered_ratios$ratio > 0.5,]
          max_each<-filtered_ratios[filtered_ratios$tissue %in% bestPtiss,] %>% group_by(frame) %>%  slice_max(meanVal, n = 1, with_ties = FALSE)
          best_t_pair <- max_each[order(max_each$meanVal, decreasing = T),]$tissue[1:2]
          best_comp_final <- paste(best_t_pair[1], " : ", best_t_pair[2], sep = "")
        } else {
          best_comp_final <- best_comparison
        }
        
      }
    }
  }
  
  # Combine results into a data frame
  result <- data.frame(
    pairID = out_pair,
    p_value = out_p,
    num_tissues = out_n,
    num_tests = out_N,
    num_frames = out_frames,
    max_diff = out_maxdif,
    pattern = out_pattern,
    best_comp = best_comp_final,
    stringsAsFactors = FALSE
  )
  
  output = list("Z" = Z, "Scount" = Scount, "result" = result, "mat" = test_f_mat)
  return(output$result)
}
# Helper function for pattern tests
pattern_test <- function(data, metric, threshold = 0.6) {
  #print(metric)
  if(metric == "CPM"){
    summary_data <- data %>%
      group_by(tissue, frame) %>%
      summarise(meanVal = median(count), .groups = "drop") %>%
      group_by(tissue) %>%
      mutate(total = sum(meanVal), ratio = meanVal / total)
  } else {
    summary_data <- data %>%
      group_by(tissue, frame) %>%
      #summarise(meanVal = mean(Mcount), .groups = "drop") %>%
      summarise(meanVal = median(count), .groups = "drop") %>%
      group_by(tissue) %>%
      mutate(total = sum(meanVal), ratio = meanVal / total)
  }
  ifelse(length(unique(summary_data %>% filter(ratio > threshold) %>% pull(frame))) > 1, "***", "NS")
}
# Modified XMulti functions:
xmulti2<-function (obs, expr, statName = "LLR", histobins = F, histobounds = c(0,0),
                   showCurve = T, detail = 1, safety = 1e+09){
  library(XNomial)
  if (length(obs) != length(expr)) 
    stop("\nThere must be the same number of observed and expected values\n")
  if (length(obs) < 2) 
    stop("\nThere must be at least two categories\n")
  ntotal <- sum(obs)
  extotal <- sum(expr)
  if (ntotal <= 0) 
    stop("\nThe total observations must be positive")
  if (extotal <= 0) 
    stop("\nThe sum of the expected numbers must be positive")
  maxprob <- max(expr/extotal)
  statType = which(c("LLR", "Prob", "Chisq") == statName)
  if (ntotal * log(maxprob) <= log(.Machine$double.xmin)){ 
    # replace with xmonte here
    #stop("\nThe numbers of observations are too large causing underflow error. The monte carlo version, \"xmonte\" is recommended for this case.")
    cat("using xmonte...")
    
    return(xmonte2(obs, expr, ntrials = 1e+07, detail = 1))
    stop()
  } else {
    if (histobins == 1) {
      histobins <- 500
    }
    if (histobounds[[2]] == 0) {
      histobounds[[2]] <- qchisq(0.999, length(obs) - 1)
    }
    ntrials <- XNomial:::ntables(obs)
    #print(ntrials)
    if (ntrials > safety) {
      cat("Full enumeration requires examination of", ntrials, 
          "tables.\n")
      stop("This operation could take more than several minutes.\n    • The monte carlo version, \"xmonte\" is recommended for this case. \n    • To override this cutoff, change the parameter 'safety' to something greater than the required number of trials.")
    }
    result = .C("exactMultinomialTest", obs = as.integer(obs), 
                expr = as.double(expr), nn = as.integer(length(obs)), 
                statType = as.integer(statType), pLLR = as.double(0), 
                pProb = as.double(0), pChi = as.double(0), observedLLR = as.double(XNomial:::logLRmultinomial(obs, 
                                                                                                              expr)), observedProb = as.double(dmultinom(obs, prob = expr)), 
                observedChi = as.double(XNomial:::chiStat(obs, expr)), histobins = as.integer(histobins), 
                histobounds = as.double(histobounds), histoData = as.double(double(if (histobins) {
                  histobins
                } else {
                  1
                })), PACKAGE = "XNomial")
    pval <- c(result$pLLR, result$pProb, result$pChi)[statType]
    if (detail == 1) {
      cat("\nP value (", statName, ") = ", formatC(pval), "\n", 
          sep = "")
    }
    if (detail >= 2) {
      cat("\nP value  (LLR)  = ", formatC(result$pLLR))
      cat("\nP value (Prob)  = ", formatC(result$pProb))
      cat("\nP value (Chisq) = ", formatC(result$pChi))
    }
    if (detail >= 3) {
      cat("\n\nObserved: ", obs, "\nExpected ratio: ", expr, 
          "\nTotal number of tables: ", ntrials, "\n")
    }
    if (histobins) {
      statHistogram(result, showCurve)
    }
    value <- result
    value[["nn"]] <- NULL
    if (histobins == 0) {
      value[["histobins"]] <- NULL
      value[["histobounds"]] <- NULL
      value[["histoData"]] <- NULL
    }
    value$statType <- statName
    asymptotic.p.value <- pchisq(result$observedChi, result$nn - 
                                   1, lower.tail = F)
    value <- c(value, asymptotic.p.value = asymptotic.p.value, 
               cases.examined = ntrials)
    invisible(value)
  } # end of if not xmonte
}
xmonte2<-function (obs, expr, ntrials = 1e+05, statName = "LLR", histobins = F, 
                   histobounds = c(0, 0), showCurve = T, detail = 1, safety = 1e+08){
  if (length(obs) != length(expr)) 
    stop("\nThere must be the same number of observed and expected values\n")
  if (length(obs) < 2) 
    stop("\nThere must be at least two categories\n")
  if (ntrials > safety) 
    stop("\nTo run this many trials, you must increase the 'safety' parameter\n")
  statType = which(c("LLR", "Prob", "Chisq") == statName)
  if (length(statType) == 0) {
    statType <- 4
  }
  if (histobins == 1) {
    histobins <- 500
  }
  if (histobounds[[2]] == 0) {
    histobounds[[2]] <- qchisq(0.999, length(obs) - 1)
  }
  invisible(rmultinom(4, obs, expr))
  result = .C("montenomialTest", obs = as.integer(obs), expr = as.double(expr), 
              ntrials = as.integer(ntrials), nn = as.integer(length(obs)), 
              statType = as.integer(statType), pLLR = as.double(0), 
              pProb = as.double(0), pChi = as.double(0), observedLLR = as.double(XNomial:::logLRmultinomial(obs, 
                                                                                                            expr)), observedProb = as.double(dmultinom(obs, prob = expr)), 
              observedChi = as.double(XNomial:::chiStat(obs, expr)), histobins = as.integer(histobins), 
              histobounds = as.double(histobounds), histoData = as.integer(integer(if (histobins) {
                histobins
              } else {
                1
              })), PACKAGE = "XNomial")
  pval <- switch(statType, result$pLLR, result$pProb, result$pChi)
  standardError <- sqrt(pval * (1 - pval)/ntrials)
  if (detail >= 1) {
    cat("\nP value (", statName, ") = ", pval, " ± ", formatC(standardError), 
        sep = "")
  }
  if (detail >= 2) {
    cat("\n", formatC(ntrials), " random trials\n Observed: ", 
        obs, "\n Expected Ratio: ", expr, sep = " ")
  }
  if (histobins) {
    statHistogram(result, showCurve)
  }
  value <- result
  value[["nn"]] <- NULL
  if (histobins == 0) {
    value[["histobins"]] <- NULL
    value[["histobounds"]] <- NULL
    value[["histoData"]] <- NULL
  }
  if (statType != 1) {
    value[["pLLR"]] <- NULL
  }
  if (statType != 2) {
    value[["pProb"]] <- NULL
  }
  if (statType != 3) {
    value[["pChi"]] <- NULL
  }
  value$statType <- statName
  asymptotic.p.value <- pchisq(result$observedChi, result$nn - 
                                 1, lower.tail = F)
  value <- c(value, asymptotic.p.value = asymptotic.p.value)
  value <- append(value, c(standard.error = standardError), 
                  after = which(names(value) == "statType") + 1)
  invisible(value)
}  

# analysis
outAllframes = DCE.nomial.tests(done)
write.table(outAllframes, "output/03_GTEx/output_DCE_All_frames.tsv", quote = F, row.names = F)

# add-on: Normalization for plotting in CPM
library(edgeR)
TMMfactors<-as.data.frame(calcNormFactors(counts[,-which(names(counts) %in% c("needles"))]))
TMMfactors$sample<-rownames(TMMfactors)
names(TMMfactors)<-c("NormFactors", "sample")
write.table(merge(gtexMetadata, TMMfactors, by = "sample"), 
            "output/03_GTEx/GTExMetadataUpdated.txt", quote = F, row.names = F)
