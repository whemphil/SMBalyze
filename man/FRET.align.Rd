\name{FRET.align}
\alias{FRET.align}
\title{Image Alignment}
\description{
A function to align reference images for dual camera single-molecule microscopy experiments.
}
\usage{
FRET.align(path.to.file='./',file.name=NULL,alignment=list('r'=0.97,'theta'=7),search.radius=6,integration.radius=6,spot.min=NULL,spot.max=Inf)
}
\arguments{
\item{path.to.file}{Path to the directory where the import file(s) are located. Write-out files will also be exported to this directory. DEFAULT = (working directory)}
\item{file.name}{A character string vector of length two, indicating the names of the tiff files containing the stacked reference images to be aligned -- should be in the order Cy3-camera then Cy5-camera images. DEFAULT = (all tiff files in directory, sorted alphabetically)}
\item{alignement}{Either a list containing numeric values for the variables 'r' and 'theta', which will be used as the starting parameter estimates for regression-based auto alignment, or the character string 'manual', which will default to no alignment adjustment and begin an iterative manual alignment. DEFAULT = (auto-alignment with empirically determined starting estimates)}
\item{search.radius}{A single numerical value, representing the box radius used to scan for local maxima during spot identification.  DEFAULT = 6.}
\item{integration.radius}{A single numerical value, representing the radius (number of pixels) around spots' maxima to be numerically integrated. DEFAULT = 6.}
\item{spot.min}{A single numerical value, representing the minimum integration volume for identification as a spot. DEFAULT = (5-fold the scaled image background).}
\item{spot.max}{A single numerical value, representing the maximum integration volume for identification as a spot. DEFAULT = Inf.}
}
\value{
Returns nothing.

Files:

\item{Alignment_Parameters.RData}{Contains the alignment parameters 'r' and 'theta' -- used to align images in downstream analyses.}
}
