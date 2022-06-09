\name{id.spots}
\alias{id.spots}
\title{Spot Identifier}
\description{
A function to identify particles of interest from the stacked tiff microscopy image file of a single-molecule biophysics experiment, and to export signal vs time data for those particles.
}
\usage{
id.spots(time.step,path.to.file='./',file.name='*.tif',spot.box=6,spot.radius=6,spot.min=NULL,spot.max=Inf,spot.picking='composite')
}
\arguments{
\item{time.step}{A single numerical value, representing the exposure time-step between images in the input file.}
\item{path.to.file}{Path to the directory where the import file(s) are located. Write-out files will also be exported to this directory. DEFAULT = (working directory)}
\item{file.name}{The name of the tiff file that contains the stacked images to be analyzed. DEFAULT = (any tiff file)}
\item{spot.box}{A single numerical value, representing the box radius used to scan for local maxima during particle identification.  DEFAULT = 6.}
\item{spot.radius}{A single numerical value, representing the radius (number of pixels) around particles' maxima to be numerically integrated. DEFAULT = 6.}
\item{spot.min}{A single numerical value, representing the minimum integration volume for identification as a particle. DEFAULT = 200 (spot.picking='comopsite') or 0 (spot.picking='all.frames').}
\item{spot.max}{A single numerical value, representing the maximum integration volume for identification as a particle. DEFAULT = Inf.}
\item{spot.picking}{A single character string, designating the spot-picking approach to be used; can be 'composite' (spots are identified on a composite image and filtered by integration volume on the composite image) or 'all.frames' (spots are identified on a composite image then filtered by integration volume on the composite image and dynamic integration volume across all image slices). DEFAULT = 'composite'.}
}
\value{
Returns a list with the following components: composite image data = a numerical matrix containing data for a composite image of parent microscopy images, initial particle summary = a data frame containing indices and volume information for the identified particles, particle traces = a matrix providing signal vs time traces for all identified particles, & particle images data = an array providing image data for identified particles across time frames.

Files:

\item{initial_particle_summary.txt}{Contains indices and volume information for the identified particles.}
\item{initial_particle_traces.txt}{Contains signal vs time traces for all identified particles.}
\item{Initial-Particle_Data.RData}{Contains information from the analysis, including composite image data, initial particle summary data, particle image data, exposure time between images, number of frames, pixel dimensions of images, and radius of integration box.}
}
