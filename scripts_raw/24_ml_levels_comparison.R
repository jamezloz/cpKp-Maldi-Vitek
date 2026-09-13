suppressPackageStartupMessages({
  library(pROC); library(data.table); library(e1071); library(class)
  library(rpart); library(randomForest); library(ada); library(xgboost)
  library(glmnet)
})

BASE     <- "resultados/RAW_modular"
OBJ_DIR  <- file.path(BASE, "objects")
OUT_DIR  <- file.path(BASE, "tables/levels_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_ITER <- as.integer(Sys.getenv("ML_N_ITER", "30"))
TAG    <- Sys.getenv("ML_TAG", "")
OUT_FILE <- file.path(OUT_DIR, sprintf("ml_levels_iterations%s.csv", TAG))
SUM_FILE <- file.path(OUT_DIR, sprintf("ml_levels_summary%s.csv", TAG))
MODELS <- c("rf","xgb","svm","knn","dt","extra_trees","adaboost","logistic","naive_bayes")

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

cv_extended <- function(X, y, groups, model_name, n_iter = N_ITER) {
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
      iter=i, model=model_name, balanced_accuracy=ba, auc=auc,
      sensitivity=sens, specificity=spec, n_features=length(keep)
    )
  }
  rbindlist(rows)
}

mat_rep <- readRDS(file.path(OBJ_DIR, "raw_peak_matrix.rds"))
det_rep <- readRDS(file.path(OBJ_DIR, "raw_peak_detected.rds"))
meta    <- as.data.table(readRDS(file.path(OBJ_DIR, "raw_metadata.rds")))
meta[, sample := basename(sample)]
stopifnot(nrow(mat_rep) == nrow(meta))
cat(sprintf("Replicados: %d espectros, %d cepas, %d features\n\n",
            nrow(mat_rep), uniqueN(meta$strain_id), ncol(mat_rep)))

aggregate_level <- function(rows, keys) {
  idx_by <- split(rows, keys[rows])
  X <- do.call(rbind, lapply(idx_by, function(ix) colMeans(mat_rep[ix, , drop = FALSE])))
  M <- rbindlist(lapply(idx_by, function(ix) meta[ix[1], .(strain_id, treatment,
                                                           mechanism, label_res_sen)]))
  list(X = X, meta = M)
}

build_level <- function(name) {
  all_rows <- seq_len(nrow(meta))
  switch(name,
    replicate        = list(X = mat_rep,
                            meta = meta[, .(strain_id, treatment, mechanism, label_res_sen)]),
    strain_treatment = aggregate_level(all_rows,
                         paste(meta$strain_id, meta$treatment, sep = "__")),
    strain_full      = aggregate_level(all_rows, meta$strain_id),
    col_only         = aggregate_level(which(meta$treatment == "COL"), meta$strain_id),
    ext_only         = aggregate_level(which(meta$treatment == "EXT"), meta$strain_id)
  )
}

LEVELS <- strsplit(Sys.getenv("ML_LEVELS",
  "replicate,strain_treatment,strain_full,col_only,ext_only"), ",")[[1]]
cat(sprintf("Niveles: %s | iteraciones: %d\n\n", paste(LEVELS, collapse=", "), N_ITER))
if (file.exists(OUT_FILE)) file.remove(OUT_FILE)
set.seed(42)

for (nm in LEVELS) {
  lv <- build_level(nm)
  M  <- lv$meta
  cat(sprintf("===== %s: %d obs, %d cepas =====\n", nm, nrow(lv$X), uniqueN(M$strain_id)))

  y_rs <- factor(M$label_res_sen, levels = c("SEN","RES"))
  for (mod in MODELS) {
    res <- cv_extended(lv$X, y_rs, M$strain_id, model_name = mod)
    if (nrow(res)) {
      res[, `:=`(task = "res_sen", nivel = nm,
                 n_obs = nrow(lv$X), n_cepas = uniqueN(M$strain_id))]
      fwrite(res, OUT_FILE, append = file.exists(OUT_FILE))
      cat(sprintf("  res_sen  %-12s BA=%.3f AUC=%.3f\n", mod,
                  mean(res$balanced_accuracy, na.rm=TRUE), mean(res$auc, na.rm=TRUE)))
    }
  }

  keep_kn <- M$mechanism %in% c("KPC","NDM")
  y_kn <- factor(M$mechanism[keep_kn], levels = c("NDM","KPC"))
  Xk   <- lv$X[keep_kn, , drop = FALSE]
  gk   <- M$strain_id[keep_kn]
  for (mod in MODELS) {
    res <- cv_extended(Xk, y_kn, gk, model_name = mod)
    if (nrow(res)) {
      res[, `:=`(task = "kpc_ndm", nivel = nm,
                 n_obs = nrow(Xk), n_cepas = uniqueN(gk))]
      fwrite(res, OUT_FILE, append = file.exists(OUT_FILE))
      cat(sprintf("  kpc_ndm  %-12s BA=%.3f AUC=%.3f\n", mod,
                  mean(res$balanced_accuracy, na.rm=TRUE), mean(res$auc, na.rm=TRUE)))
    }
  }
}

r <- fread(OUT_FILE)
s <- r[, .(n_obs = n_obs[1], n_cepas = n_cepas[1],
           BA = round(mean(balanced_accuracy, na.rm=TRUE), 3),
           AUC = round(mean(auc, na.rm=TRUE), 3),
           sens = round(mean(sensitivity, na.rm=TRUE), 3),
           spec = round(mean(specificity, na.rm=TRUE), 3),
           n_feats = round(mean(n_features, na.rm=TRUE), 1)),
       by = .(task, nivel, model)][order(task, nivel, -BA)]
fwrite(s, SUM_FILE)

cat("\n===== MEJOR MODELO POR NIVEL =====\n")
print(s[, .SD[which.max(BA)], by = .(task, nivel)][order(task, -BA)])
cat("\nGuardado en:", OUT_DIR, "\n")
