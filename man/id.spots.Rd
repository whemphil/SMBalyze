\name{id.spots}
\alias{id.spots}
\title{Spot Identifier}
\description{
A function to identify particles of interest from the stacked tiff microscopy image file of a single-molecule biophysics experiment, and to export signal vs time data for those particles.
}
\usage{
id.spots(path.to.file,file.name,time.step,spot.box=6,spot.radius=6,spot.min=300,spot.max=Inf,spot.picking='composite',r.pick=2.9,min.pick=5)
}
\arguments{
\item{path.to.file}{Path to the directory where the import file(s) are located. Write-out files will also be exported to this directory.}
\item{file.name}{The name of the tiff file that contains the stacked images to be analyzed.}
\item{time.step}{A single numerical value, representing the exposure time-step between images in the input file.}
\item{spot.box}{A single numerical value, representing the box radius used to scan for local maxima during particle identification.  DEFAULT = 6.}
\item{spot.radius}{A single numerical value, representing the radius (number of pixels) around particles' maxima to be numerically integrated. DEFAULT = 6.}
\item{spot.min}{A single numerical value, representing the minimum integration volume for identification as a particle. DEFAULT = 300.}
\item{spot.max}{A single numerical value, representing the maximum integration volume for identification as a particle. DEFAULT = Inf.}
\item{spot.picking}{A single character string, designating the spot-picking approach to be used; can be 'composite' (an average of all image frames is used for spot identification) or 'all.frames' (spots are identified on each image frame then clustered across time to identify conserved spots). DEFAULT = 'composite'.}
\item{r.pick}{A single numerical value, representing the radius (in pixels) for spot-merging if 'spot.picking'='all.frames'. DEFAULT = 2.9.}
\item{min.pick}{A single numerical value, representing the minimum number of frames a spot must be identified in to pass as conserved. DEFAULT = 5.}
}
\value{
Returns a list with the following components: composite image data = a numerical matrix containing data for a composite image of parent microscopy images, initial particle summary = a data frame containing indices and volume information for the identified particles, particle traces = a matrix providing signal vs time traces for all identified particles, & particle images data = an array providing image data for identified particles across time frames.

Files:

\item{initial_particle_summary.txt}{Contains indices and volume information for the identified particles.}
\item{initial_particle_traces.txt}{Contains signal vs time traces for all identified particles.}
\item{Initial-Particle_Data.RData}{Contains information from the analysis, including composite image data, initial particle summary data, particle image data, exposure time between images, number of frames, pixel dimensions of images, and radius of integration box.}
}
