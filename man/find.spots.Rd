\name{find.spots}
\alias{find.spots}
\title{Microscopy Image Spot Finder}
\description{
Function to identify bivariate local maxima and approximate their volume.
}
\usage{
find.spots(data,box.size,low.lim,high.lim,fill.radius)
}
\arguments{
\item{data}{A numerical matrix, where row and column axes are the two variables of interest.}
\item{box.size}{A single numerical value, representing the radius (number of data points) of a box around each focal data point while scanning for local maxima.}
\item{low.lim}{A single numerical matrix, representing the minimum approximated integration volume around a local maximum necessary to pass identification.}
\item{high.lim}{A single numerical matrix, representing the maximum approximated integration volume around a local maximum necessary to pass identification.}
\item{fill.radius}{A single numerical value, representing the radius (number of data points) of a box around each local maximum that numerical integration should be performed within.}
}
\value{
Returns a data frame with the following columns: x = x-coordinates of identified maxima for on an 'image' plot of the parent matrix, y = y-coordinates of identified maxima for on an 'image' plot of the parent matrix, row = row indices of identified maxima, col = column indices of identified maxima, & vol = approximated integration volume of identified maxima.
}
