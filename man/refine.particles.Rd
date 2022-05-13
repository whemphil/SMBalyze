\name{refine.particles}
\alias{refine.particles}
\title{Particle Refinement and Analysis}
\description{
Function to analyze bining states of particles from single molecule biophysics microscopy experiments.
}
\usage{
refine.particles(path.to.file,file.name='Initial-Particle_Data.RData',skip.manual='n',signal.step=1000)
}
\arguments{
\item{path.to.file}{Path to the directory where the import file(s) are located. Write-out files will also be exported to this directory.}
\item{file.name}{The name of the tiff file that contains the stacked images to be analyzed.}
\item{time.step}{A single numerical value, representing the exposure time-step between images in the input file.}
\item{spot.box}{A single numerical value, representing the 'box.size' argument in the 'find.spots' function.  DEFAULT = 6.}
\item{spot.radius}{A single numerical value, representing the radius (number of pixels) around particles' maxima to be numerically integrated. Feeds into the 'fill.radius' argument in the 'find.spots' function. DEFAULT = 6.}
\item{spot.min}{A single numerical value, representing the minimum integration volume for identification as a particle. Feeds into the 'low.lim' argument of the 'find.spots' function. DEFAULT = 300.}
\item{spot.max}{A single numerical value, representing the maximum integration volume for identification as a particle. Feeds into the 'high.lim' argument of the 'find.spots' function. DEFAULT = Inf.}
}
\value{
Returns a list with the following components:  image.avg = a numerical matrix containing data for a composite image of parent microscopy images, residence.data = a data frame describing manually identified binding states for selected particles, including columns for particle id, start time, end time, residence time, and signal average for each state, refined.particle.traces = a matrix containing signal vs time data for each selected particle, refined.particle.trace.rolls = a matrix containing smoothed signal vs time data for each selected particle, refined.spots = a data frame containing indices and volume information for all selected particles, refined.particle.snaps = an array providing image data for selected particles across time frames, state.calls = a numerical matrix containing automatically called state vs time information for all initially identified particles, residence.calls = a data frame containing the automatically called sequential binary binding states and residence times for each initially identified particle, & dwell.calls = a data frame containing the automatically called sequential multi-binding states and residence times for each initially identified particle.

Files:

\item{Refined-Particle_Data.RData}{An RData file containing all necessary results and input for subsequent functions.}
\item{residence_data.txt}{A table describing manually identified binding states for selected particles, including columns for particle id, start time, end time, residence time, and signal average for each state.}
\item{selected_particle_traces.txt}{A table containing signal vs time data for each selected particle.}
\item{selected_particle_smoothed-traces.txt}{A table containing smoothed signal vs time data for each selected particle.}
\item{selected_particle_summary.txt}{A table containing indices and volume information for all selected particle composites.}
\item{all-particle_state-calls.txt}{A table containing automatically called state vs time information for all initially identified particles.}
\item{all-particle_dwell-calls.txt}{A table containing the automatically called sequential multi-binding states and residence times for each initially identified particle.}
\item{all-particle_residence-calls.txt}{A table containing the automatically called sequential binary binding states and residence times for each initially identified particle.}
}

