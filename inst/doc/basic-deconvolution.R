## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache=TRUE,cache.path='c/',fig.path='fig/')

## ----cache=TRUE---------------------------------------------------------------
library('hspe')
names(shen_orr_ex)

## ----cache=TRUE---------------------------------------------------------------
?hspe

## ----cache=TRUE---------------------------------------------------------------
data = shen_orr_ex$data$log[,c(1:10,201:210,401:410)]
mixture_proportions = shen_orr_ex$annotation$mixture

## ----cache=TRUE---------------------------------------------------------------
pure_samples = list(Liver=c(1,2,3),Brain=c(4,5,6),Lung=c(7,8,9))
out = hspe(Y=data, pure_samples = pure_samples,optim_opts=list(pkg="nloptr"))

## ----cache=TRUE---------------------------------------------------------------
true_proportions = mixture_proportions[-(1:9),]
matplot(true_proportions,out$estimates, xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

## ----cache=TRUE---------------------------------------------------------------
mixture_samples = data[-(1:9),]
reference_samples = data[1:9,]

out = hspe(Y=mixture_samples, reference=reference_samples,pure_samples = pure_samples)

matplot(true_proportions,out$estimates, xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

## ----cache=TRUE---------------------------------------------------------------
ref_reduced = t(sapply(pure_samples,function(x)colMeans(reference_samples[x,,drop=FALSE])))

out = hspe(Y=mixture_samples, reference=ref_reduced)

matplot(true_proportions,out$estimates, xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=mixture_samples, references = ref_reduced)

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=mixture_samples, references = ref_reduced,marker_method = "diff")

## ----cache=TRUE---------------------------------------------------------------
out$n_markers

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=3)

out$n_markers

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=c(1,2,3))

out$n_markers

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=.075)

out$n_markers

out = hspe(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=c(.1,.15,.05))

out$n_markers

## ----cache=TRUE---------------------------------------------------------------
marker_genes = list(c(1,2,3),
                    c(4,5,6),
                    c(7,8,9))

out = hspe(Y=mixture_samples, references = ref_reduced,markers=marker_genes)
out$n_markers

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y = mixture_samples,references = ref_reduced,markers=mrkrs,n_markers=.1)

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=lin_scale_mix,references = lin_scale_ref,inv_scale = base::identity,
                  seed=1234,markers = mrkrs$L)
head(out$estimates)

## ----cache=TRUE---------------------------------------------------------------
ahs_scale_mix = asinh(lin_scale_mix)
ahs_scale_ref = asinh(lin_scale_ref)

out = hspe(Y=ahs_scale_mix,references = ahs_scale_ref,inv_scale = base::sinh,
                  seed=1234,markers = mrkrs$L)
head(out$estimates)

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=lin_scale_mix,references = lin_scale_ref,inv_scale = base::identity,fit_scale=base::sqrt,seed=1234,markers = mrkrs$L)
head(out$estimates)

## ----cache=TRUE---------------------------------------------------------------
out = hspe(Y=lin_scale_mix,references = lin_scale_ref,inv_scale = base::identity,fit_scale=base::log,loss_fn="L2",seed=1234,markers = mrkrs$L)
head(out$estimates)

