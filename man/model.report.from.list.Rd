% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_linear_models.R
\name{model.report.from.list}
\alias{model.report.from.list}
\title{Generate a model report and a plot from a model list (after model.subset)}
\usage{
model.report.from.list(
  dataset,
  model.list,
  out.col = "output",
  leave.out = "",
  predict = F,
  what.model = NULL
)
}
\arguments{
\item{dataset}{a dataframe with outcome column (must be named 'output')}

\item{model.list}{minimum # of features (default = 2)}

\item{out.col}{column number of y value}

\item{leave.out}{name of observations to leave out (e.g. 'p_Br') and retrain the model}

\item{predict}{if leave.out is not empty, should a prediction
of it be computed?}

\item{what.model}{if not in an interactive session, what model should be used?}
}
\value{
models list, CV results for k=3/5/LOO and a plot of choice. interactive.
}
\description{
Screen for models, cross validate them and plot. Designed for interactive work.
}
