\name{classify.states}
\alias{classify.states}
\title{Classification of Signal States}
\description{
Function to classify single-molecule biophysics signal vs time data into distinct states.
}
\usage{
classify.states(data,signal.step=1000)
}
\arguments{
\item{data}{A numerical vector (time-series) of data to be classified.}
\item{signal.step}{A single numerical value, representing the difference in signal values between states. DEFAULT = 1000.}
}
\value{
Returns a data frame with the following columns: id = The state classification id of each data point in the original data vector, & avg = The average signal value for the state id of each data point in the original data vector.
}
