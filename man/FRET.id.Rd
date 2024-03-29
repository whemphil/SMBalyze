\name{FRET.id}
\alias{FRET.id}
\title{Dual-Camera Spot Identifier}
\description{
A function to independently identify and pair particles of interest from the stacked tiff image files of a single-molecule microscopy experiment using two cameras/filters, and to export signal vs time data for those particles.
}
\usage{
FRET.id(time.step,path.to.file='./',Cy3.file.name=NULL,Cy5.file.name=NULL,align.file.name='Alignment_Parameters.RData',Cy3.search.box=6,Cy5.search.box=6,Cy3.integration.radius=6,Cy5.integration.radius=6,Cy3.spot.min=3500,Cy5.spot.min=1e4,Cy3.spot.max=Inf,Cy5.spot.max=Inf,border.spots=F)}
\arguments{
\item{time.step}{A single numerical value, representing the exposure time-step between images in the input file.}
\item{path.to.file}{A single character string indicating the path to the directory where the import file(s) are located. Write-out files will also be exported to this directory. DEFAULT = (working directory)}
\item{Cy3.file.name}{A single character string indicating the name of the tiff file that contains the stacked images from camera 1 to be analyzed. NOTE -- these designations should match the file name designations used for the FRET.align function. DEFAULT = (first alphabetical tiff file of all tiff files in directory)}
\item{Cy5.file.name}{A single character string indicating the name of the tiff file that contains the stacked images from camera 2 to be analyzed. NOTE -- these designations should match the file name designations used for the FRET.align function. DEFAULT = (second alphabetical tiff file of all tiff files in directory)}
\item{align.file.name}{A single character string indicating the name of the alignment file generated by the FRET.align function. DEFAULT = 'Alignment_Parameters.RData'}
\item{Cy3.search.box}{A single numerical value, representing the box radius used to scan for local maxima during spot identification on the Cy3-camera images.  DEFAULT = 6.}
\item{Cy5.search.box}{A single numerical value, representing the box radius used to scan for local maxima during spot identification on the Cy5-camera images.  DEFAULT = 6.}
\item{Cy3.integration.radius}{A single numerical value, representing the radius (number of pixels) around spots' maxima to be numerically integrated on the Cy3-camera images. DEFAULT = 6.}
\item{Cy5.integration.radius}{A single numerical value, representing the radius (number of pixels) around spots' maxima to be numerically integrated on the Cy5-camera images. DEFAULT = 6.}
\item{Cy3.spot.min}{A single numerical value, representing the minimum integration volume for identification as a spot on the Cy3-camera images. DEFAULT = 3500.}
\item{Cy5.spot.min}{A single numerical value, representing the minimum integration volume for identification as a spot on the Cy5-camera images. DEFAULT = 10000.}
\item{Cy3.spot.max}{A single numerical value, representing the maximum integration volume for identification as a spot on the Cy3-camera images. DEFAULT = Inf.}
\item{Cy5.spot.max}{A single numerical value, representing the maximum integration volume for identification as a spot on the Cy5-camera images. DEFAULT = Inf.}
\item{border.spots}{A logical, indicating if spots close to the image border should be included. DEFAULT = FALSE.}
}
\value{
Returns a list with the following components each for the Cy3-camera and Cy5-camera: composite image data = a numerical matrix containing data for a composite image of parent microscopy images, initial particle summary = a data frame containing indices and volume information for the identified particles, particle traces = a matrix providing signal vs time traces for all identified particles, & particle images data = an array providing image data for identified particles across time frames.

Files:

\item{initial_particle_Cy3-summary.txt}{Contains indices on the Cy3-images for the identified particles.}
\item{initial_particle_Cy3-traces.txt}{Contains signal vs time traces on the Cy3-images for all identified particles.}
\item{initial_particle_Cy5-summary.txt}{Contains indices on the Cy5-images for the identified particles.}
\item{initial_particle_Cy5-traces.txt}{Contains signal vs time traces on the Cy5-images for all identified particles.}
\item{Initial-Particle_Data.RData}{Contains information from the analysis, with Cy3-images and Cy5-images versions where appropriate, including composite image data, initial particle summary data, particle image data, exposure time between images, number of frames, pixel dimensions of images, and radius of integration box.}
\item{Identified_Spots.pdf}{An image depicting the identified particles' locations on Cy3-camera and Cy5-camera images. NOTE -- green markers were only identified on the Cy3-camera images, red markers on the Cy5-camera images, and blue markers on both.}
}
