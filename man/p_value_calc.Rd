% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{p_value_calc}
\alias{p_value_calc}
\title{P value calculation function for models used in comparison_plot() function}
\usage{
p_value_calc(
  data_path,
  window_size = 3,
  model,
  parameter,
  p_val_method,
  repetition
)
}
\arguments{
\item{data_path}{The path to the data file.}

\item{window_size}{The size of the window for analysis.}

\item{model}{The model to be used for p-value calculation.}

\item{parameter}{The parameters for the model.}

\item{p_val_method}{The method to be used for p-value calculation ("Parametric" or "Non-Parametric").}

\item{repetition}{The number of repetitions for the analysis.}
}
\value{
A list containing the calculated p-values and other related information.
}
\description{
P value calculation function for models used in comparison_plot() function
}
