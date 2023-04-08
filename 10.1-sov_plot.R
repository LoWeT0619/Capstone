library(shiny)
library(ggplot2)
library(parallel)
library(dplyr)
library(lumi)
library(minfi)
library(RColorBrewer)
library(MASS)
library(limma)
library(ggcorrplot)
library(carData)
library(car)

load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")

#extract confounders
pd = pd%>%dplyr::select("Sample_Group", "Age", "Parity", "BMI", "GA", "Smoker", "Eth2")

sovdat = data.frame(t(Mvalues))
# sovdatlist <- list()
# for(i in 1:ncol(sovdat)){
#   sovdatlist[[i]] <- sovdat[,i]
# }

Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)

for(i in 2:ncol(pd)){
  varname <- names(pd)[i]
  Ftab[varname] <- numeric()
}

calF <- function(probe = probecol){
  newdata <- pd
  pdnames <- names(newdata)
  newdata$beta <- probe
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  fit <- lm(formstr, data = newdata)
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  F <- aovfit$`F value`
  F <- F[2:(length(F)-1)]
  names(F) <- pdnames
  F <- as.data.frame(F, stringsAsFactors = FALSE)
  F <- as.data.frame(t(F))
  row.names(F) <- 1
  Ftab <- rbind(Ftab, F)
  return(Ftab)
}

Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 40)
Ftab <- do.call(rbind, Ftab)
Fmean <- colMeans(Ftab)
Fmean <- Fmean[order(-Fmean)]
Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)
finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))

save(Ftab, Fmean, finalvars, file = "res2.RData")

load("/home/liuwent/08-SOV_Plot/res2.RData")

sovplot <- function(restab = MSSmean, clustername = 'Case', plottype = 'MSS',
                    textsize = 20){
  resmean <- restab
  samplegroupidx <- match('Sample_Group', resmean$Factor)
  resmean$Factor[samplegroupidx] <- paste0(clustername, '_Control')
  if(plottype == 'MSS'){
    ytitle <- 'Mean Square'
    resmean <- resmean[order(-resmean$MSSstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))
    )
  }else if(plottype == 'pval'){
    ytitle <- '-log2(p-val)'
    resmean <- resmean[order(-resmean$logpval),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))
    )
  }else{
    ytitle <- 'F statistic'
    resmean <- resmean[order(-resmean$Fstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))+
        geom_hline(yintercept=1, linetype="dashed", color = "red")
    )
  }
}

pdf("/nfs/home/liuwent/08-SOV_Plot/04-sov_plot_1219.pdf")
sovplot(restab = Fmean, plottype = 'F', textsize = 7)
dev.off()

Fmean
#The variables that are over 1 are the confounding factors you will address.