# ML LOO on Summarized Experiment
# J.H Badger (as part of McCulloch, et al. 2021)
# November 14, 2021


# Load needed libraries

library(tidyverse)
library(SummarizedExperiment)
library(futile.logger)
library(caret)
library(MLeval)

# summarized experiment object to analyze

lkt <- readRDS("exp.rds")

runML <- function(exp, response="Study_Clin_Response", 
                  positive="Responder", method="rf", index=NULL,
                  Group = "", seed=15, boot=100, metric=max) {
  set.seed(seed)
  train_control <- trainControl(method="cv", number=10, 
                                classProbs=T, savePredictions = T)
  if (is.null(index)) {
    index <- createDataPartition(exp[[response]], p=0.70, list=FALSE)
  }
  models <- list()
  accuracies <- c()
  for(i in 1:boot) {
    training <- exp[, as.vector(index)]
    tData <- t(assay(training))
    tResp <- factor(training[[response]])
    models[[i]] <- train(tData, tResp, trControl = train_control, 
                   method = method, preProcess=c("center", "scale"))
    testing <- exp[,-as.vector(index)] %>% assay() %>% t()
    obs <- factor(exp[,row.names(testing)][[response]])
    predictions <- factor(predict(models[[i]], testing))
    matrix <- confusionMatrix(predictions, obs,
                  positive = positive)
    accuracies[i] <- matrix$overall["Accuracy"]
  }
  accuracies <- round(accuracies, 3)
  target_acc <- metric(accuracies)
  model <- models[[which(accuracies==target_acc)[1]]]
  if (is_null(model)) {
    flog.error("null model detected")
    model <- models[[which(!is.null(models))[1]]]
  }
  ppred <- data.frame(predict(model, testing, type="prob"))
  ppred$obs <- obs
  ppred$Group <- Group
  roc <- evalm(ppred, silent=TRUE, plots="r", showplots = FALSE)
  rstats <- roc$optres[[1]]
  stats <- matrix$overall[1:6]
  stats <- append(stats, rstats["AUC-ROC","Score"])
  names(stats)[length(stats)] <- "AUC"
  list(stats=stats,model=model)
}

testModel <- function(model, exp, seed=15, response="Study_Clin_Response",  positive="Responder", Group = "") {
  set.seed(seed)
  testing <- assay(exp) %>% t()
  obs <- factor(exp[[response]])
  predictions <- factor(predict(model, testing))
  matrix <- confusionMatrix(predictions, obs,
                    positive = positive)
  ppred <- data.frame(predict(model, testing, type="prob"))
  ppred$obs <- obs
  roc <- evalm(ppred, silent=TRUE, plots="r", showplots = FALSE)
  rstats <- roc$optres[[1]]
  stats <- matrix$overall[1:6]
  stats <- append(stats, rstats["AUC-ROC","Score"])
  ppred$Group <- str_c(Group, " p = ", signif(stats[["AccuracyPValue"]],3))
  names(stats)[length(stats)] <- "AUC"
  list(stats=stats, ppred=ppred)
}

results <- NULL
pdf("loo_rocs.pdf")
for(study in sort(unique(lkt$Study))) {
  predictions <- data.frame()
  for(method in c("rf","svmPoly","glm")) {
    flog.info(str_c(study," ", method))
    mdl <- runML(lkt[,lkt$Study != study], method = method)
    res <- testModel(mdl$model, lkt[,lkt$Study==study], Group=method)
    predictions <- rbind(predictions, res$ppred)
    line <- append(c(Study=study, Method=method), res$stats)
    if (is.null(results)) {
      results <- t(data.frame(line))
    } else {
      results <- rbind(results, line)
    }
  }
  roc <- evalm(predictions, silent = TRUE, plots="r", showplots = FALSE)$roc + 
    ggtitle(str_c("LKT ", study, " bootstrapped"))
  print(roc)
}
dev.off()
write_tsv(data.frame(results),"results.tsv")

