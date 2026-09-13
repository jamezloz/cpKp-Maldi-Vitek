source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages({
  library(pROC); library(data.table); library(e1071); library(class)
  library(rpart); library(randomForest); library(ada); library(xgboost)
  library(glmnet)
})

RAW_MODULAR_OBJ <- "resultados/RAW_modular/objects"
RAW_MODULAR_OUT <- "resultados/RAW_modular"
OUT_FILE <- file.path(RAW_MODULAR_OUT, "tables", "raw_extended_classifiers_strain_treatment.csv")

sample_group_split <- function(groups, train_frac = 0.8) {
  ug <- unique(groups)
  tr_g <- sample(ug, size = floor(train_frac * length(ug)))
  list(train = groups %in% tr_g, test = !groups %in% tr_g)
}

select_features_wilcox <- function(X, y, fdr_thresh = 0.05) {
  pvals <- apply(X, 2, function(x) tryCatch(wilcox.test(x ~ y)$p.value, error = function(e) 1))
  fdr   <- p.adjust(pvals, method = "fdr")
  keep  <- which(fdr < fdr_thresh)
  if (length(keep) == 0) keep <- order(pvals)[seq_len(min(10, ncol(X)))]
  keep
}

cv_extended <- function(X, y, groups, model_name, n_iter = 100) {
  pos_label <- levels(y)[2]; neg_label <- levels(y)[1]
  rows <- list()
  for (i in seq_len(n_iter)) {
    sp  <- sample_group_split(groups)
    tr  <- sp$train; te <- sp$test
    ytr <- y[tr];    yte <- y[te]
    if (length(unique(ytr)) < 2 || length(unique(yte)) < 2) next
    keep <- select_features_wilcox(X[tr, , drop=FALSE], ytr)
    Xtr  <- X[tr, keep, drop=FALSE]; Xte <- X[te, keep, drop=FALSE]
    Xtr_sc <- scale(Xtr)
    sc_ctr <- attr(Xtr_sc,"scaled:center"); sc_scl <- attr(Xtr_sc,"scaled:scale")
    sc_scl[sc_scl==0] <- 1
    Xtr_sc[is.nan(Xtr_sc)] <- 0
    Xte_sc <- scale(Xte, center=sc_ctr, scale=sc_scl)
    probs <- NULL; pred <- NULL
    tryCatch({
      if (model_name=="svm") {
        fit<-svm(Xtr_sc,ytr,kernel="radial",probability=TRUE)
        pp<-predict(fit,Xte_sc,probability=TRUE); pred<-pp
        probs<-attr(pp,"probabilities")[,pos_label]
      } else if (model_name=="knn") {
        k<-max(1,round(sqrt(nrow(Xtr_sc))))
        pred<-knn(Xtr_sc,Xte_sc,cl=ytr,k=k,prob=TRUE)
        raw_p<-attr(pred,"prob"); probs<-ifelse(pred==pos_label,raw_p,1-raw_p)
      } else if (model_name=="dt") {
        fit<-rpart(y~.,data=data.frame(y=ytr,Xtr),method="class",
                   control=rpart.control(minsplit=5,cp=0.001))
        pp<-predict(fit,data.frame(Xte),type="prob"); probs<-pp[,pos_label]
        pred<-factor(ifelse(probs>=0.5,pos_label,neg_label),levels=levels(y))
      } else if (model_name=="extra_trees") {
        fit<-randomForest(x=Xtr,y=ytr,ntree=500,mtry=ncol(Xtr),replace=FALSE)
        probs<-predict(fit,Xte,type="prob")[,pos_label]
        pred<-predict(fit,Xte,type="response")
      } else if (model_name=="adaboost") {
        Xtr_ada<-as.data.frame(Xtr); colnames(Xtr_ada)<-paste0("f",colnames(Xtr_ada))
        Xte_ada<-as.data.frame(Xte); colnames(Xte_ada)<-paste0("f",colnames(Xte_ada))
        fit<-ada(x=Xtr_ada,y=ytr,iter=100,loss="exponential")
        probs_mat<-predict(fit,Xte_ada,type="probs")
        probs<-probs_mat[,which(levels(ytr)==pos_label)]
        pred<-factor(ifelse(probs>=0.5,pos_label,neg_label),levels=levels(y))
      } else if (model_name=="logistic") {

        cvfit<-cv.glmnet(Xtr_sc,ytr,family="binomial",alpha=0,nfolds=5)
        probs<-as.numeric(predict(cvfit,Xte_sc,s="lambda.min",type="response"))
        pred<-factor(ifelse(probs>=0.5,pos_label,neg_label),levels=levels(y))
      } else if (model_name=="naive_bayes") {
        fit<-naiveBayes(Xtr_sc,ytr)
        pp<-predict(fit,Xte_sc,type="raw"); probs<-pp[,pos_label]
        pred<-predict(fit,Xte_sc,type="class")
      } else if (model_name=="rf") {
        fit<-randomForest(x=Xtr,y=ytr,ntree=500)
        probs<-predict(fit,Xte,type="prob")[,pos_label]
        pred<-predict(fit,Xte,type="response")
      } else if (model_name=="xgb") {
        dtr<-xgb.DMatrix(data=Xtr,label=as.numeric(ytr==pos_label))
        dte<-xgb.DMatrix(data=Xte,label=as.numeric(yte==pos_label))
        fit<-xgb.train(params=list(objective="binary:logistic",eval_metric="auc",nthread=1),
                       data=dtr,nrounds=100,verbose=0)
        probs<-predict(fit,dte)
        pred<-factor(ifelse(probs>=0.5,pos_label,neg_label),levels=levels(y))
      }
    }, error=function(e) message(sprintf("  [iter %d %s]: %s",i,model_name,e$message)))
    if (is.null(pred)||is.null(probs)) next
    TP<-sum(pred==pos_label&yte==pos_label); FN<-sum(pred==neg_label&yte==pos_label)
    TN<-sum(pred==neg_label&yte==neg_label); FP<-sum(pred==pos_label&yte==neg_label)
    sens<-if((TP+FN)>0) TP/(TP+FN) else NA_real_
    spec<-if((TN+FP)>0) TN/(TN+FP) else NA_real_
    ba<-mean(c(sens,spec),na.rm=TRUE)
    auc<-tryCatch(as.numeric(pROC::auc(pROC::roc(as.numeric(yte==pos_label),probs,quiet=TRUE))),
                  error=function(e) NA_real_)
    rows[[length(rows)+1]] <- data.table(
      iter=i, model=model_name,
      balanced_accuracy=ba, auc=auc,
      sensitivity=sens, specificity=spec, n_features=length(keep)
    )
  }
  rbindlist(rows)
}

main <- function() {

  ds      <- readRDS(file.path(RAW_MODULAR_OBJ, "raw_peak_strain_treatment_dataset.rds"))
  X_str   <- ds$matrix
  meta_str <- as.data.frame(ds$metadata)

  strain_ids_block <- meta_str$strain_id

  models <- c("rf","xgb","svm","knn","dt","extra_trees","adaboost","logistic","naive_bayes")
  dir.create(dirname(OUT_FILE), recursive=TRUE, showWarnings=FALSE)
  if (file.exists(OUT_FILE)) file.remove(OUT_FILE)

  save_res <- function(res, task, level) {
    if (nrow(res)==0) return(invisible(NULL))
    res[, task:=task][, level:=level]
    fwrite(res, OUT_FILE, append=file.exists(OUT_FILE))
    cat(sprintf("  saved %d rows\n", nrow(res)))
  }

  cat("=== CRKp/CSKp — strain x treatment (539 obs) ===\n")
  y_rs <- factor(meta_str$label_res_sen, levels=c("SEN","RES"))
  for (mod in models) {
    cat(sprintf("  %s\n", mod))
    res <- cv_extended(X_str, y_rs, strain_ids_block, model_name=mod, n_iter=100)
    save_res(res, "res_sen", "strain_treatment")
  }

  cat("=== KPC/NDM — strain x treatment ===\n")
  keep_kn <- meta_str$mechanism %in% c("KPC","NDM")
  y_kn    <- factor(meta_str$mechanism[keep_kn], levels=c("NDM","KPC"))
  Xm      <- X_str[keep_kn, , drop=FALSE]
  grp_kn  <- strain_ids_block[keep_kn]
  for (mod in models) {
    cat(sprintf("  %s\n", mod))
    res <- cv_extended(Xm, y_kn, grp_kn, model_name=mod, n_iter=100)
    save_res(res, "kpc_ndm", "strain_treatment")
  }

  cat("\nDone. Results in:", OUT_FILE, "\n")

  results <- fread(OUT_FILE)
  cat("\n=== SUMMARY ===\n")
  summary_tbl <- results[, .(BA=round(mean(balanced_accuracy,na.rm=TRUE),3),
                              AUC=round(mean(auc,na.rm=TRUE),3)),
                          by=.(task,model)][order(task,-BA)]
  print(summary_tbl)
}

set.seed(42)
main()
