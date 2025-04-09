## ----setup, echo = FALSE, results = "hide", message = FALSE--------------
set.seed(290875)

dir <- system.file("applications", package = "tbm")
datadir <- system.file("rda", package = "TH.data")


sapply(c("tram", "survival", "tbm", "mboost", "lattice", 
         "latticeExtra", "partykit"), library, char = TRUE)

trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

knitr::opts_chunk$set(echo = TRUE, results = 'markup', error = FALSE,
                      warning = FALSE, message = FALSE,
                      tidy = FALSE, cache = FALSE, size = "small",
                      fig.width = 6, fig.height = 4, fig.align = "center",
                      out.width = NULL, ###'.6\\linewidth', 
                      out.height = NULL,
                      fig.scap = NA)
knitr::render_sweave()  # use Sweave environments
knitr::set_header(highlight = '')  # do not \usepackage{Sweave}
## R settings
options(prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)  # JSS style
options(width = 75)
library("colorspace")
col <- diverge_hcl(2, h = c(246, 40), c = 96, l = c(65, 90))
fill <- diverge_hcl(2, h = c(246, 40), c = 96, l = c(65, 90), alpha = .3)

## ----citation, echo = FALSE----------------------------------------------
year <- substr(packageDescription("tbm")$Date, 1, 4)
version <- packageDescription("tbm")$Version

## ----src-appl, results = "hide"------------------------------------------
system.file("applications", package = "tbm")

## ----src-sim, results = "hide"-------------------------------------------
system.file("simulations", package = "tbm")

## ----results, echo = FALSE-----------------------------------------------
ds <- c("Beetles Exctinction Risk",
        "Birth Weight Prediction",
        "Body Fat Mass",
        "CAO/ARO/AIO-04 DFS",
        "Childhood Malnutrition",
        "Head Circumference",
        "PRO-ACT ALSFRS",
        "PRO-ACT OS")
names(ds) <- c("ex_beetles.rda", "ex_fetus.rda", "ex_bodyfat.rda", "ex_CAO.rda", 
        "ex_india.rda", "ex_heads.rda", "ex_ALSFRS.rda", "ex_ALSsurv.rda")
models <- c(expression(paste("L ", bold(vartheta)(bold(x)))),
       expression(paste("L ", beta(bold(x)))),
       expression(paste("N ", bold(vartheta)(bold(x)))),
       expression(paste("N ", beta(bold(x)))),
       expression(paste("T ", bold(vartheta)(bold(x)))),
       expression(paste("T ", beta(bold(x)))),
       expression(paste("Tree ", bold(vartheta)(bold(x)))),
       expression(paste("F ", bold(vartheta)(bold(x)))))
names(models) <- c("glm_ctm", "glm_tram", "gam_ctm", "gam_tram", 
              "tree_ctm", "tree_tram", "trtf_tree", "trtf_forest")
mor <- c("gam_ctm", "glm_ctm", "tree_ctm", "gam_tram", "glm_tram", 
        "tree_tram", "trtf_tree", "trtf_forest")

pit <- function(object, data) {
    cf <- predict(object, newdata = data, coef = TRUE)
    tmp <- object$model
    ret <- numeric(NROW(data))
    for (i in 1:NROW(data)) {
        coef(tmp) <- cf[i,]
        ret[i] <- predict(tmp, newdata = data[i,,drop = FALSE], type = "distribution")
    }
    ret
}

pitplot <- function(object, data, main = "") {
    u <- pit(object, data)
    plot(1:length(u) / length(u), sort(u), xlab = "Theoretical Uniform Quantiles", 
         ylab = "PIT Quantiles", main = main, ylim = c(0, 1), xlim = c(0, 1))
    abline(a = 0, b = 1, col = "red")
} 

plot.ctm <- function(x, which = sort(unique(selected(object))), 
                     obs = TRUE, ...) {

    object <- x
    q <- mkgrid(x$model, n = 100)
    X <- model.matrix(x$model$model, data = q)
    pfun <- function(which) {
        stopifnot(length(which) == 1)
        df <- model.frame(object, which = which)[[1]]
        df <- lapply(df, function(x) seq(from = min(x), 
                                         to = max(x), length = 50))
        df2 <- do.call("expand.grid", c(df, q))
        class(object) <- class(object)[-1]
        tmp <- predict(object, newdata = df, which = which)
        CS <- diag(ncol(tmp[[1]]))
        CS[upper.tri(CS)] <- 1
        cf <- tmp[[1]] %*% CS
        df2$pr <- as.vector(t(X %*% t(cf)))
        ret <- cbind(df2, which = names(df)[1])
        names(ret)[1:2] <- c("x", "y")
        ret
    }

    ret <- do.call("rbind", lapply(which, pfun))
    at <- quantile(ret$pr, prob = 0:10/10)

    pf <- function(x, y, z, subscripts, ...) {
        xname <- as.character(unique(ret[subscripts, "which"]))
        xo <- model.frame(object, which = xname)[[1]][[xname]]
        yo <- object$response
        panel.levelplot.raster(x, y, z, subscripts, ...)
        panel.contourplot(x, y, z, subscripts, contour = TRUE, ...)
        if (obs)
            panel.xyplot(x = xo, y = yo, pch = 20)
    }
    levelplot(pr ~ x + y | which, data = ret, panel = pf,
              scales = list(x = list(relation = "free")), region = TRUE, 
              at = at, col.regions = grey.colors(length(at)), ...)
}

### source("movie.R")

## ----riskplot, eval = FALSE, echo = FALSE--------------------------------
# x <- risk[,mor] - ll0
# med <- apply(x, 2, median, na.rm = TRUE)
# boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
#         ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
#         ylim = ylim)
# abline(v = c(3.5, 6.5), col = "lightgrey")
# abline(h = 0, col = "lightgrey")
# abline(h = max(med), col = col[2])
# axis(1, at = 1:ncol(x), labels = models[colnames(x)])
# axis(2)
# box()

## ----beetles-results, echo = FALSE, fig.width = 7, fig.height = 6--------
load(file.path(dir, file <- "ex_beetles.rda"))
ylim <- c(-.5, .5)
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----beetles-setup, echo = FALSE-----------------------------------------
load(file.path(datadir, "beetles.rda"))

ldata <- bmodel
ldata$TS <- NULL

X <- model.matrix(~ species - 1, data = ldata)
ldata$species <- NULL

nvars <- c("mean_body_size", "mean_elev", "flowers",
          "niche_decay", "niche_diam", "niche_canopy", "distribution")
fvars <- colnames(ldata)[!colnames(ldata) %in% c("RL", nvars)]
xvars <- c(nvars, fvars)

ldata[nvars] <- lapply(nvars, function(x) as.numeric(scale(ldata[[x]])))

fm_gam <- c(
"ctm" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 2) + ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 2)", collapse = "+"))),
"stm" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 1) + ",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 1)", collapse = "+"))))

fm_glm <- c(
"ctm" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 2) +",
    paste("bols(", xvars, ", df = 2)", collapse = "+"))),
"stm" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 1) +",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"))))

fm_tree <- as.formula(paste("RL ~ ", paste(xvars, collapse = "+")))

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]

## ----beetles-model0, echo = TRUE, cache = TRUE---------------------------
m_mlt <- Polr(RL ~ 1, data = ldata)
logLik(m_mlt)

## ----beetles-null, echo = FALSE------------------------------------------
q <- mkgrid(m_mlt, 200)[[1]]
p <- c(predict(m_mlt, newdata = data.frame(RL = q), type = "density"))
names(p) <- as.character(q)
barplot(p, xlab = "Red List Status", ylab = "Distribution")

## ----beetles-model, cache = TRUE-----------------------------------------
fm_tree
bm <- stmboost(m_mlt, formula = fm_tree, data = ldata,
               method = quote(mboost::blackboost))[626]

## ----insampleriskplot, echo = FALSE--------------------------------------
plot(-risk(bm), xlab = "Boosting Iteration", ylab = "Log-likelihood", type = "l")

## ----beetles-prop, echo = FALSE------------------------------------------
x <- tabulate(do.call("c", sapply(get("ens", environment(bm$predict)), 
    function(x) nodeapply(x$model, ids = nodeids(x$model, terminal = FALSE), FUN = function(x) split_node(x)$varid - 1))), nbins = length(vars))
names(x) <- vars
x

## ----beetles-plot, echo = FALSE, fig.width = 6, fig.height = 8-----------
idx <- sapply(levels(ldata$RL), function(x) 
    which(ldata$RL == x)[1])
q <- mkgrid(m_mlt, n = 200)[[1]]
d <- predict(bm, newdata = ldata[idx,], type = "density", q = q)
layout(matrix(1:6, nrow = 3, byrow = TRUE))
out <- sapply(1:6, function(x) barplot(d[,x], xlab = "Red List",
              ylab = "Density", main = paste(gsub("_", " ", rownames(ldata)[idx[x]]), " (RL:", ldata$RL[idx[x]], ")", sep = "")))

## ----beetles-mlt, cache = TRUE-------------------------------------------
po <- blackboost(fm_tree, data = ldata, family = PropOdds())[626]
max(abs(risk(bm) - risk(po)))

## ----fetus-results, echo = FALSE, fig.width = 7, fig.height = 6----------
load(file.path(dir, file <- "ex_fetus.rda"))
ylim <- c(0, 1.7)
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----fetus-setup, echo = FALSE-------------------------------------------
load(file.path(datadir, "fetus.rda"))

fetus$birthweight <- as.double(fetus$birthweight)
ldata <- fetus

xvars <- colnames(ldata)
xvars <- xvars[xvars != "birthweight"]

ldata[xvars] <- lapply(xvars, function(x) as.numeric(scale(ldata[[x]])))


fm_gam <- c("ctm" = as.formula(paste("birthweight ~ ",
                paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
                paste("bbs(", xvars, ", center = TRUE, df = 2)", collapse = "+"))),
            "stm" = as.formula(paste("birthweight ~ ", 
                paste("bols(", xvars, ", intercept = FALSE)", collapse = "+"), "+",
                paste("bbs(", xvars, ", center = TRUE, df = 1)", collapse = "+"))))
fm_glm <- c("ctm" = as.formula(paste("birthweight ~ ",
                paste("bols(", xvars, ", df = 2)", collapse = "+"))),
            "stm" = as.formula(paste("birthweight ~ ", 
                paste("bols(", xvars, ", intercept = FALSE)", collapse = "+"))))

fm_tree <- birthweight ~ .

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]

## ----fetus-model0, cache = TRUE------------------------------------------
m_mlt <- BoxCox(birthweight ~ 1, data = ldata, extrapolate = TRUE)
logLik(m_mlt)

## ----fetus-null, echo = FALSE--------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "density", K = 200,
     xlab = "Birth Weight (in gram)", ylab = "Density", col = "black")

## ----fetus-model, cache = TRUE-------------------------------------------
fm_gam[["stm"]]
bm <- stmboost(m_mlt, formula = fm_gam[["stm"]], data = ldata,
               method = quote(mboost::mboost))[253]

## ----fetus-risk, echo = FALSE--------------------------------------------
plot(-risk(bm), xlab = "Boosting Iteration", ylab = "Log-likelihood", type = "l")

## ----fetus-pit, echo = FALSE, fig.width = 5, fig.height = 5--------------
pitplot(bm, ldata, main = ds[file])

## ----fetus-prop, echo = FALSE--------------------------------------------
tab <- matrix(tabulate(selected(bm), nbins = length(bm$baselearner)), ncol = 1)
rownames(tab) <- names(bm$baselearner)
colnames(tab) <- ""
tab

## ----fetus-plot, echo = FALSE, fig.width = 6, fig.height = 6-------------
w <- sort(unique(selected(bm)))
mbm <- bm
class(mbm) <- class(mbm)[-1]
nc <- ceiling(length(w) / 2) * 2
layout(matrix(1:nc, ncol = 2))
plot(mbm, which = w, ylim = c(-3, 3))

## ----fetus-base, echo = FALSE--------------------------------------------
tmp <- as.mlt(m_mlt)
coef(tmp) <- nuisance(bm)
plot(tmp, newdata = data.frame(1), K = 200, type = "trafo", xlab = "Birth weight (in gram)", 
     ylab = "Transformation h", col = "black")

## ----fetus-dens, echo = FALSE--------------------------------------------
idx <- order(ldata$birthweight)[floor(c(.1, .25, .5, .75, .9) * nrow(ldata))]
q <- mkgrid(m_mlt, n = 200)[[1]]
matplot(x = q, y = d <- predict(bm, newdata = ldata[idx,], type = "density", q = q), type =
        "l", col = rgb(.1, .1, .1, .6), lty = 1, xlab = "Birth Weight (in gram)", ylab = "Density")
j <- sapply(ldata[idx, "birthweight"], function(y) which.min((q - y)^2))
points(ldata[idx, "birthweight"], sapply(1:5, function(i) d[j[i], i]), pch = 19, col = col[1])

## ----bodyfat-results, echo = FALSE, fig.width = 7, fig.height = 6--------
load(file.path(dir, file <- "ex_bodyfat.rda"))
ylim <- c(0, 1.7)
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----bodyfat-setup, echo = FALSE-----------------------------------------

data("bodyfat", package = "TH.data")
ldata <- bodyfat

xvars <- colnames(ldata)
xvars <- xvars[xvars != "DEXfat"]

ldata[xvars] <- lapply(xvars, function(x) as.numeric(scale(ldata[[x]])))


fm_gam <- c("ctm" = as.formula(paste("DEXfat ~ ", 
                paste("bbs(", xvars, ")", collapse = "+"))),
            "stm" = as.formula(paste("DEXfat ~ ", 
                paste("bols(", xvars, ", intercept = FALSE)", collapse = "+"), "+",
                paste("bbs(", xvars, ", center = TRUE, df = 1)", collapse = "+"))))
fm_glm <- c("ctm" = as.formula(paste("DEXfat ~ ",
                paste("bols(", xvars, ")", collapse = "+"))),
            "stm" = as.formula(paste("DEXfat ~ ", 
                paste("bols(", xvars, ", intercept = FALSE)", collapse = "+"))))

fm_tree <- DEXfat ~ .

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]

## ----bodyfat-model0, cache = TRUE----------------------------------------
m_mlt <- BoxCox(DEXfat ~ 1, data = ldata, prob = c(.1, .99))
logLik(m_mlt)

## ----bodyfat-null, echo = FALSE------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "density", K = 200,
     xlab = "Body Fat Mass (in kilogram)", ylab = "Density", col = "black")

## ----bodyfat-model, cache = TRUE-----------------------------------------
fm_gam[["ctm"]]
bm <- ctmboost(m_mlt, formula = fm_gam[["ctm"]], data = ldata,
               method = quote(mboost::mboost))[1000]

## ----bodyfat-risk, echo = FALSE------------------------------------------
plot(-risk(bm), xlab = "Boosting Iteration", ylab = "Log-likelihood", type = "l")

## ----bodyfat-pit, echo = FALSE, fig.width = 5, fig.height = 5------------
pitplot(bm, ldata, main = ds[file])

## ----bodyfat-ctmplot, echo = FALSE, fig.width = 6, fig.height = 6, dev = "png"----
plot.ctm(bm, ylab = "Body Fat Mass")

## ----bodyfat-prop, echo = FALSE------------------------------------------
tab <- matrix(tabulate(selected(bm), nbins = length(bm$baselearner)), ncol = 1)
rownames(tab) <- names(bm$baselearner)
colnames(tab) <- ""
tab

## ----bodyfat-dens, echo = FALSE------------------------------------------
idx <- order(ldata$DEXfat)[floor(c(.1, .25, .5, .75, .9) * nrow(ldata))]
q <- mkgrid(m_mlt, n = 200)[[1]]
matplot(x = q, y = d <- predict(bm, newdata = ldata[idx,], type = "density", q = q), type =
        "l", col = rgb(.1, .1, .1, .6), lty = 1, xlab = "Body Fat Mass (in kilogram)", ylab = "Density")
j <- sapply(ldata[idx, "DEXfat"], function(y) which.min((q - y)^2))
points(ldata[idx, "DEXfat"], sapply(1:5, function(i) d[j[i], i]), pch = 19, col = col[1])

## ----CAO-results, echo = FALSE, fig.width = 7, fig.height = 6------------
load(file.path(dir, file <- "ex_CAO.rda"))
ylim <- NULL
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----CAO-setup, echo = FALSE---------------------------------------------
load(file.path(datadir, "Primary_endpoint_data.rda"))

ldata <- CAOsurv

xvars <- c("strat_t", "strat_n", "randarm",
           "geschlecht", "age", "groesse", "gewicht",
           "ecog_b", "mason",
           "gesamt_t", "gesamt_n", 
           "histo", "grading_b", "ap_anlage",
           "stenose",
           "gesamt_n_col", "UICC", "bentf")

ldata <- ldata[, c("DFS", xvars)]
ldata <- ldata[complete.cases(ldata),]
ldata$UICC <- ldata$UICC[,drop = TRUE]
ldata$ecog_b <- ldata$ecog_b[,drop = TRUE]


nvars <- c("age", "groesse", "gewicht")
fvars <- xvars[!xvars %in% nvars]
xvars <- c(nvars, fvars)

ldata[nvars] <- lapply(nvars, function(x) as.numeric(scale(ldata[[x]])))

fm_gam <- c(
"ctm" = as.formula(paste("DFS ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 2)", collapse = "+"))),
"stm" = as.formula(paste("DFS ~ ",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 1)", collapse = "+"))))

fm_glm <- c(
"ctm" = as.formula(paste("DFS ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"))),
"stm" = as.formula(paste("DFS ~",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"))))

fm_tree <- as.formula(paste("DFS ~ ", paste(xvars, collapse = "+")))

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]

## ----CAO-model0, cache = TRUE--------------------------------------------
m_mlt <- Coxph(DFS ~ 1, data = ldata, prob = c(0, .9))
logLik(m_mlt)

## ----CAO-null, echo = FALSE----------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "survivor", K = 200,
     xlab = "Time (in days)", ylab = "Time", col = "black")

## ----CAO-model, cache = TRUE, eval = FALSE-------------------------------
# fm_glm[["stm"]]
# bm <- stmboost(m_mlt, formula = fm_glm[["stm"]], data = ldata,
#                method = quote(mboost::mboost))[963]

## ----CAO-mlt, cache = TRUE, eval = FALSE---------------------------------
# bmCoxPH <- mboost(fm_glm[["stm"]], data = ldata, family = CoxPH())[963]

## ----india-results, echo = FALSE, fig.width = 6, fig.height = 6----------
load(file.path(dir, file <- "ex_india.rda"))
ylim <- NULL
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----india-setup, echo = FALSE-------------------------------------------
load(file.path(datadir, "india.rda"))

xvars <- c("cage", "breastfeeding", "mbmi", "mage", "medu", 
           "edupartner", "csex", "ctwin", "cbirthorder", 
           "munemployed", "mreligion", "mresidence", "deadchildren",
           "electricity", "radio", "television", "refrigerator", 
           "bicycle", "motorcycle", "car")

fvars <- xvars[sapply(xvars, function(x) length(unique(kids[[x]]))) < 6]
nvars <- xvars[!xvars %in% fvars]
kids[fvars] <- lapply(fvars, function(f) factor(kids[[f]]))
kids[nvars] <- lapply(nvars, function(n) scale(kids[[n]]))

fm_gam <- c(
"ctm" = as.formula(paste("stunting ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 2)", collapse = "+"))),
"stm" = as.formula(paste("stunting ~",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 1)", collapse = "+"))))

fm_glm <- c(
"ctm" = as.formula(paste("stunting ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"))),
"stm" = as.formula(paste("stunting ~",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"))))

fm_tree <- as.formula(paste("stunting ~ ", paste(xvars, collapse = "+")))

kids$stunting <- as.double(kids$stunting)
ldata <- kids

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]


## ----india-model0, cache = TRUE------------------------------------------
m_mlt <- BoxCox(stunting ~ 1, data = ldata, prob = c(.05, .975), extrapolate = TRUE)
logLik(m_mlt)

## ----india-null, echo = FALSE--------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "density", K = 200,
     xlab = "Stunting", ylab = "Density", col = "black")

## ----india-model, cache = TRUE, eval = FALSE-----------------------------
# fm_tree
# bm <- ctmboost(m_mlt, formula = fm_tree, data = ldata,
#                method = quote(mboost::blackboost))[515]

## ----heads-results, echo = FALSE, fig.width = 7, fig.height = 6----------
load(file.path(dir, file <- "ex_heads.rda"))
ylim <- NULL
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----heads-setup, echo = FALSE-------------------------------------------
### Paul Eilers: Using 1/2 is better
data("db", package = "gamlss.data")
db$lage <- with(db, age^(1/3))

ldata <- db

fm_gam <- c("ctm" = head ~ bbs(lage),
            "stm" = head ~ bols(lage, intercept = FALSE) + bbs(lage, center = TRUE, df = 1))
fm_glm <- c("ctm" = head ~ bols(lage),
            "stm" = head ~ bols(lage, intercept = FALSE))

fm_tree <- head ~ .
N <- nrow(ldata)

## ----heads-model0, cache = TRUE------------------------------------------
m_mlt <- BoxCox(head ~ 1, data = ldata)
logLik(m_mlt)

## ----heads-model, cache = TRUE-------------------------------------------
fm_gam[["ctm"]]
bm <- ctmboost(m_mlt, formula = fm_gam[["ctm"]], data = ldata, 
               method = quote(mboost::mboost))[1000]

## ----heads-risk, echo = FALSE, cache = TRUE------------------------------
plot(-risk(bm), xlab = "Boosting Iteration", ylab = "Log-likelihood", type = "l")

## ----heads-pit, echo = FALSE, fig.width = 5, fig.height = 5, cache = TRUE, dev = "png"----
pitplot(bm, ldata, main = ds[file])

## ----heads-plot, fig.width = 6, fig.height = 4, echo = FALSE, dev = "png"----
l <- with(db, seq(from = min(lage), to = max(lage), length = 100))
q <- mkgrid(m_mlt, n = 200)[[1]]
pr <- expand.grid(head = q, lage = l^3)
pr$p <- c(predict(bm, newdata = data.frame(lage = l), q = q, type = "distribution"))
pr$cut <- factor(pr$lage > 2.5)
levels(pr$cut) <- c("Age < 2.5 yrs", "Age > 2.5 yrs")

pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, 
        at = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6)/ 100, ...)
    panel.xyplot(x = db$age, y = db$head, pch = 20,
                 col = rgb(.1, .1, .1, .1), ...)
}
print(contourplot(p ~ lage + head | cut, data = pr, panel = pfun, region = FALSE,
            xlab = "Age (years)", ylab = "Head circumference (cm)", 
            scales = list(x = list(relation = "free"))))

## ----ALSFRS-results, echo = FALSE, fig.width = 7, fig.height = 6---------
load(file.path(dir, file <- "ex_ALSFRS.rda"))
ylim <- NULL
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----ALSFRS-setup, echo = FALSE------------------------------------------
load(file.path(datadir, "ALSFRSdata.rda"))
ldata <- ALSFRSdata[complete.cases(ALSFRSdata),]

xvars <- names(ldata)
xvars <- xvars[xvars != "ALSFRS.halfYearAfter"]
fvars <- xvars[sapply(ldata[xvars], is.factor)]
nvars <- xvars[!xvars %in% fvars]

ldata[nvars] <- lapply(nvars, function(x) as.numeric(scale(ldata[[x]])))

fm_gam <- c(
"ctm" = as.formula(paste("ALSFRS.halfYearAfter ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 2)", collapse = "+"))),
"stm" = as.formula(paste("ALSFRS.halfYearAfter ~",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 1)", collapse = "+"))))

fm_glm <- c(
"ctm" = as.formula(paste("ALSFRS.halfYearAfter ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"))),
"stm" = as.formula(paste("ALSFRS.halfYearAfter ~",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"))))

fm_tree <- as.formula(paste("ALSFRS.halfYearAfter ~ ", paste(xvars, collapse = "+")))

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]


## ----ALSFRS-model-0, cache = TRUE----------------------------------------
m_mlt <- Colr(ALSFRS.halfYearAfter ~ 1, data = ldata, prob = c(.05, .9), extrapolate = TRUE)
logLik(m_mlt)

## ----ALSFRS-null, echo = FALSE-------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "density", K = 200,
     xlab = "ALSFRS at 6 months", ylab = "Density", col = "black")

## ----ALSFRS-model, cache = TRUE------------------------------------------
bm <- ctmboost(m_mlt, formula = fm_glm[["ctm"]], data = ldata, 
               method = quote(mboost::mboost))[613]

## ----ALSFRS-risk, echo = FALSE-------------------------------------------
plot(-risk(bm), xlab = "Boosting Iteration", ylab = "Log-likelihood", type = "l")

## ----ALSFRS-prop, echo = FALSE-------------------------------------------
tab <- matrix(tabulate(selected(bm), nbins = length(bm$baselearner)), ncol = 1)
rownames(tab) <- names(bm$baselearner)
colnames(tab) <- ""
tab[tab > 0,]

## ----ALSFRS-plot, results = "hide", echo = FALSE-------------------------
cf <- mboost:::coef.mboost(bm)
q <- mkgrid(m_mlt, n = 100)[[1]]
X <- model.matrix(m_mlt$model, data = data.frame(ALSFRS.halfYearAfter = q))
vars <- c("time_onset_treatment", "ALSFRS.Start") 
layout(matrix(1:length(vars), ncol = 2))
out <- sapply(vars, function(v) {
    cf <- matrix(cf[[grep(v, names(cf))]], ncol = 2, byrow = TRUE)
    cf <- cumsum(cf[,2])
    plot(q, X %*% cf, main = v, xlab = "ALSFRS at 6 months",
         ylab = "Varying coefficient", type = "l")
})

## ----ALSsurv-results, echo = FALSE, fig.width = 7, fig.height = 6--------
load(file.path(dir, file <- "ex_ALSsurv.rda"))
ylim <- NULL
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = models[colnames(x)])
axis(2)
box()

## ----ALSsurv-setup, echo = FALSE-----------------------------------------
load(file.path(datadir, "ALSsurvdata.rda"))

ldata <- ALSsurvdata[complete.cases(ALSsurvdata),]
ldata$y <- with(ldata, Surv(survival.time, cens))
ldata$survival.time <- NULL
ldata$cens <- NULL

nvars <- c("time_onset_treatment", "age", "height")
fvars <- names(ldata)[!names(ldata) %in% c("y", nvars)]
xvars <- c(nvars, fvars)

ldata[nvars] <- lapply(nvars, function(x) as.numeric(scale(ldata[[x]])))


fm_gam <- c(
"ctm" = as.formula(paste("y ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 2)", collapse = "+"))),
"stm" = as.formula(paste("y ~ ",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 1)", collapse = "+"))))

fm_glm <- c(
"ctm" = as.formula(paste("y ~ ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"))),
"stm" = as.formula(paste("y ~",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"))))

fm_tree <- as.formula(paste("y ~ ", paste(xvars, collapse = "+")))

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]

## ----ALSsurv-model0, cache = TRUE----------------------------------------
m_mlt <- Coxph(y ~ 1, data = ldata)
logLik(m_mlt)

## ----ALSsurv-null, echo = FALSE------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "survivor", K = 200,
     xlab = "ALS Overall Survival (in days)", ylab = "Probability", col = "black")

## ----ALSsurv-model, cache = TRUE, eval = FALSE---------------------------
# fm_tree
# bm <- stmboost(m_mlt, formula = fm_tree, data = ldata,
#                control = boost_control(nu = 0.01),
#                method = quote(mboost::blackboost))[451]

## ----ALSsurv-mlt, cache = TRUE, eval = FALSE, eval = FALSE---------------
# blackboost(fm_tree, data = ldata, family = CoxPH())

## ----dgp-----------------------------------------------------------------
library("tram")
### set-up RNG
set.seed(27031105)
### order of Bernstein polynomials
order <- 6
### support of distributions
sup <- c(-4, 6)
### bounds (essentially for plots)
bds <- c(-8, 8)
### shift effects: main effect of grp, x, and interaction effect
beta <- c(-2, 2, -1)

## ----dgp-data------------------------------------------------------------
### two groups
grp <- gl(2, 1)
### uniform predictor
x <- seq(from = 0, to = 1, length.out = 1000)
d <- expand.grid(grp = grp, x = x)
### transformation of x: sin(2 pi x) (1 + x)
d$g <- with(d, sin(x * 2 * pi) * (1 + x))
### generate some response (this is NOT the 
### response used for the simulations)
X <- model.matrix(~ grp * x, data = d)[,-1]
d$y <- rnorm(nrow(d), mean = X %*% beta)

## ----dgp-Linear-beta-----------------------------------------------------
### h_Y
(cf0 <- seq(from = sup[1], to = sup[2], length = order + 1) + 
        sin(seq(from = 0, to = 2*pi, length = order + 1)) * 
        (seq(from = 0, to = 2*pi, length = order + 1) <= pi) * 2)
m1 <- BoxCox(y ~ grp * x, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m1)
(cf[] <- c(cf0, beta))
coef(m1) <- cf

## ----dgp-Nonlinear-beta--------------------------------------------------
m2 <- BoxCox(y ~ grp * g, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m2)
cf[] <- c(cf0, beta)
coef(m2) <- cf

## ----dgp-Linear-parm-----------------------------------------------------
m3 <- BoxCox(y | grp * x ~ 1, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m3)

### beta_1
(cf1 <- sin(seq(from = 0, to = pi / 2, length.out = order + 1)) * beta[1])
### beta_2
(cf2 <- sqrt(seq(from = 0, to = 2, length.out = order + 1)) / sqrt(2) * beta[2])
### beta_3
(cf21 <- sin(seq(from = 0, to = pi / 2, length.out = order + 1)) * beta[3])

cf[] <- c(cf0, cf1, cf2, cf21)
coef(m3) <- cf

## ----dgp-Nonlinear-parm--------------------------------------------------
m4 <- BoxCox(y | grp * g ~ 1, data = d, order = order, 
             model_only = TRUE, support = sup, bounds = bds)
cf <- coef(m4)

cf[] <- c(cf0, cf1, cf2, cf21)
coef(m4) <- cf

## ----plot-dgp, echo = FALSE,  fig.width = 6, fig.height = 8, dev = "png"----

xlim <- c(-5, 5)
q <- mkgrid(m1, n = 500)[["y"]]
x <- 0:10 / 10
nd <- expand.grid(grp = grp, x = x)
ndq <- expand.grid(y = q, grp = grp, x = x)
ndq$xconfig <- with(ndq, interaction(grp, factor(x)))
nd$g <- with(nd, sin(x * 2 * pi) * (1 + x))

ndq$d1 <- c(predict(m1, newdata = nd, type = "density", q = q))
ndq$d2 <- c(predict(m2, newdata = nd, type = "density", q = q))
ndq$d3 <- c(predict(m3, newdata = nd, type = "density", q = q))
ndq$d4 <- c(predict(m4, newdata = nd, type = "density", q = q))

ndq$t1 <- c(predict(m1, newdata = nd, type = "trafo", q = q))
ndq$t2 <- c(predict(m2, newdata = nd, type = "trafo", q = q))
ndq$t3 <- c(predict(m3, newdata = nd, type = "trafo", q = q))
ndq$t4 <- c(predict(m4, newdata = nd, type = "trafo", q = q))

colr <- rep(grey.colors(length(x), alpha = .9), each = 2)
xlim <- bds

ret <- do.call("rbind", list(ndq)[rep(1, 4)])
ret$d <- with(ndq, c(d1, d2, d3, d4))
ret$t <- with(ndq, c(t1, t2, t3, t4))
ret$model <- gl(4, nrow(ndq))
levels(ret$model) <- mlev <- c("GLM Tram", "GAM Tram", "DR", "CTM")
ret$model <- factor(as.character(ret$model), levels = rev(mlev), labels = rev(mlev))
levels(ret$grp) <- c("Group 1", "Group 2")

lev <- c(expression(paste("Nonlinear ", bold(vartheta)(bold(x)))),
         expression(paste("Linear ", bold(vartheta)(bold(x)))),
         expression(paste("Nonlinear ", beta(bold(x)))),
         expression(paste("Linear ", beta(bold(x)))))

sc <- function(which.given, ..., factor.levels) {
    if (which.given == 1) {
        strip.default(which.given = which.given, ..., factor.levels)
    } else {
        if (which.given == 2)
            strip.default(which.given = 2, ..., factor.levels = lev)
   }
}


obj <- (xyplot(t ~ y | grp + model, data = ret, type = "l", groups = xconfig, col = colr, strip = sc,
     lwd = 2, lty = 1, xlim = xlim, layout = c(4, 2), ylab = "Transformation h", xlab = "Response Y", 
     scales = list(y = list(relation = "free"))))

obj <- useOuterStrips(obj, strip.left = function(which.given, ..., factor.levels) 
    strip.custom(horizontal = FALSE)(which.given = which.given, ..., factor.levels = lev))

plot(obj <- update(obj, legend = list(right = list(fun = "draw.colorkey", 
                   args = list(list(at = 1:22 / 22, col = colr, title = "x"))))))

trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression(x), 0.5, 1.07, hjust=0.5, vjust=1)
trellis.unfocus()

## ----simulate, eval = FALSE----------------------------------------------
# y1 <- simulate(m1, newdata = d, n = 100)
# y2 <- simulate(m2, newdata = d, n = 100)
# y3 <- simulate(m3, newdata = d, n = 100)
# y4 <- simulate(m4, newdata = d, n = 100)

## ----dgp-setup, echo = FALSE---------------------------------------------
dir <- system.file("simulations", package = "tbm")
ctm <- list.files(path = dir, pattern = "ctm.*rda", full = TRUE)
tram <- list.files(path = dir, pattern = "tram.*rda", full = TRUE)

out <- c()

for (f in c(ctm, tram)) {

    load(f)

    out <- rbind(out, ret)

}

out$model <- factor(out$model, levels = c("ctm_gam", "ctm_glm", "ctm_tree", "tram_gam", "tram_glm", "tram_tree"),
                    labels = c("gam_ctm", "glm_ctm", "tree_ctm", "gam_tram", "glm_tram", "tree_tram"))

out$PNON <- factor(out$PNON)
nlab <- paste("N =", nlev <- rev(sort(unique(out$NOBS))))
out$NOBS <- factor(out$NOBS, levels = nlev, labels = nlab)
out$m <- factor(out$m)
out$todistr <- factor(out$todistr)
out$order <- factor(out$order)

out$boost_ll[grep("rror", out$boost_ll)] <- NA
out$boost_ll <- as.double(out$boost_ll)
out$true_ll <- as.double(out$true_ll)
out$mlt_ll <- as.double(out$mlt_ll)
out$mstop <- as.double(out$mstop)

out$todistr <- factor(as.character(out$todistr), levels = c("normal", "logistic", "minextrval"), 
                      labels = c("Normal", "Logistic", "MEV"))

lev <- rev(c(expression(paste("Nonlinear ", bold(vartheta)(bold(x)))),
         expression(paste("Linear ", bold(vartheta)(bold(x)))),
         expression(paste("Nonlinear ", beta(bold(x)))),
         expression(paste("Linear ", beta(bold(x))))))

sc <- function(which.given, ..., factor.levels) {
    if (which.given == 2) {
        strip.default(which.given = which.given, ..., factor.levels)
    } else {
        if (which.given == 1)
            strip.default(which.given = 1, ..., factor.levels = lev)
   }
}

  mypanel <- function(x, y, groups, subscripts, ...) {
    med <- tapply(y, x, median, na.rm = TRUE)
    panel.bwplot(x = x, y = y, col = col[(med == max(med)) + 1], 
                 fill = col[(med == max(med)) + 1], pch = "|", ...)
    panel.abline(v = c(3.5), col = "lightgrey")
    panel.abline(h = 0, col = "lightgrey")
    panel.abline(h = max(med), col = col[2])
    panel.bwplot(x = x, y = y, col = col[(med == max(med)) + 1], 
                 fill = col[(med == max(med)) + 1], pch = "|", ...)
#    tapply(1:length(y), groups[subscripts], function(i) {
#        llines(x = 1:nlevels(x), y = y[i][order(x[i])], 
#               col = rgb(0.1, 0.1, 0.1, 0.1))
#    })
  }




## ----dgp-1, results = "asis", echo = FALSE-------------------------------
lto <- levels(out$todistr)
lp <- levels(out$PNON)
lo <- levels(out$order)

ret <- c()

i <- 1
for (l1 in lto) {
    for(l2 in lp) {
       for (l3 in lo) {
           os <- subset(out, todistr == l1 & PNON == l2 & order == l3)
           os <- subset(os, boost_ll - true_ll > -10000)

           txt <- paste("J = ", l2, "M = ", l3)
           if (l1 == "Normal") {
               main <- expression(F[Z] == Phi)
           } else if (l1 == "Logistic") {
               main <- expression(F[Z] == expit)
           } else {
               main <- expression(F[Z] == MEV)
           }
           pdf(fig <- paste("sfig", i, ".pdf", sep = ""), width = 12, height = 8)
           pl <- bwplot(I(boost_ll - true_ll) ~ model |  m + NOBS, data = os, main = main, strip = sc, 
                        ylab = "Out-of-sample log-likelihood (centered)",
                        xlab = txt, panel = mypanel,
                        scales = list(y = list(relation = "free"),
                                      x = list(at = 1:6, label = models[levels(os$model)])))
           plot(pl)
           dev.off()
           ret <- c(ret, "\\begin{sidewaysfigure}",
               "\\begin{center}",
               paste("\\includegraphics{", fig, "}", sep = ""),
               "\\end{center}",
               "\\end{sidewaysfigure}", " ")
           i <- i + 1
        }
    }
}
writeLines(ret)

## ----sessionInfo, echo = FALSE, results = "hide"-------------------------
sessionInfo()

