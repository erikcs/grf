# National Study of Learning Mindsets
data.all = read.csv("../acic18/synthetic_data.csv")
data.all$schoolid = factor(data.all$schoolid)
DF = data.all[,-1]
school.id = as.numeric(data.all$schoolid)
school.mat = model.matrix(~ schoolid + 0, data = data.all)
school.size = colSums(school.mat)
Y = DF$Y
X.raw = DF[,-(1:2)]
C1.exp = model.matrix(~ factor(X.raw$C1) + 0)
XC.exp = model.matrix(~ factor(X.raw$XC) + 0)
X = cbind(X.raw[,-which(names(X.raw) %in% c("C1", "XC"))], C1.exp, XC.exp)
data.NSLM = list(X = X, Y = Y)

# Angrist & Evans
data = read.csv("../aos19/angrist_evans/angev80_recode_run1_line525.csv.xz")
FEATURES=data.frame(data$twoa.agem,
                    data$twoa.agefstm,
                    data$twoa.educm,
                    as.numeric(data$twoa.blackm),
                    as.numeric(data$twoa.hispm),
                    as.numeric(data$twoa.othracem),
                    twoa.incomed=round(data$twoa.incomed))
#labor income: twoa.incomem, worked for pay: twoa.workedm
DF.all=data.frame(
  X=FEATURES,
  Y=1 - as.numeric(data$twoa.workedm), # The outcome is whether the mother did not work
  W=as.numeric(data$twoa.kidcount > 2),
  I=as.numeric(data$twoa.samesex))

# remove ~15% of data with missing father's income
# roughly 4% of fathers have zero income after removing missing
#
# only consider married women
is.ok = !is.na(data$twoa.incomed) & (data$twoa.marital==0)
DF=DF.all[is.ok,1:8]
data.angrist = list(X = DF[, -8], Y = DF[,8])

#' Simulate regression data
#'
#' Regression benchmarks based on simulations DGPs in grf::generate_causal_data as well as
#' public datasets from
#' 1) National Study of Learning Mindsets from the 2018 Atlantic Causal Inference Conference
#' train: 8 000, test: 2 391
#' 2) Labour decisions dataset from Angrist & Evans (1998)
#' train: 280 000, test: 54 535
#'
generate_regression_data <- function(n, p,
                                     dgp = c("kunzel", "ai1", "ai2", "NSLM", "AngristEvans"),
                                     ntrain = NULL) {
  dgp <- match.arg(dgp)
  if (dgp %in% c("kunzel", "ai1", "ai2")) {
    out = grf::generate_causal_data(n, p, dgp = dgp)
  } else if (dgp == "NSLM") {
    if (is.null(ntrain)) {
      ntrain = sample(1:nrow(data.NSLM$X), 8000)
      out = list(X=data.NSLM$X[ntrain, ], Y=data.NSLM$Y[ntrain], ntrain=ntrain)
    } else {
      out = list(X=data.NSLM$X[-ntrain, ], Y=data.NSLM$Y[-ntrain], ntrain=NULL)
    }
  } else if (dgp == "AngristEvans") {
    if (is.null(ntrain)) {
      ntrain = sample(1:nrow(data.angrist$X), 280000)
      out = list(X=data.angrist$X[ntrain, ], Y=data.angrist$Y[ntrain], ntrain=ntrain)
    } else {
      out = list(X=data.angrist$X[-ntrain, ], Y=data.angrist$Y[-ntrain], ntrain=NULL)
    }
  }

  out
}
