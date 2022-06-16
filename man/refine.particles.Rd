\name{refine.particles}
\alias{refine.particles}
\title{Particle Refinement and Analysis}
\description{
Function to analyze binding states of particles from single molecule biophysics microscopy experiments.
}
\usage{
refine.particles(path.to.file='./',file.name='Initial-Particle_Data.RData',skip.manual='n',signal.step=NULL,auto.filter='none',classification.strategy='classic',background.subtraction='lower.quartile')
}
\arguments{
\item{path.to.file}{Path to the directory where the import file(s) are located. Write-out files will also be exported to this directory. DEFAULT = (working directory)}
\item{file.name}{The name of the file that contains the selected particle traces to be analyzed. DEFAULT = 'Initial-Particle_Data.RData'}
\item{skip.manual}{A single character string, either 'y' (yes) or 'n' (no), defining if manual particle selection should be skipped.  DEFAULT = 'n'.}
\item{signal.step}{A single numerical value, representing the change in signal associated with a single binding event when the 'classification.strategy' argument is set to 'supervised'. DEFAULT = NULL, which approximates the parameter from signal variance in the data.}
\item{auto.filter}{A single character string, defining which particle traces should be automatically filtered out if skip.manual='n'. Can be 'unbound', to filter out particle traces that remain in an unbound state the entire time, or 'all.stable', to filter out particle traces that remain in ANY single state the entire time. DEFAULT = 'none'.}
\item{classification.strategy}{A single character string, designating the approach used for binding state classification from signal vs time data; can be 'classic', 'uni.dbscan', or 'var.shift'. DEFAULT = 'classic'.}
\item{background.subtraction}{A single character string, designating the approach used for subtracting background from the signal vs time plots to set unbound state averages to ~0 signal; can be 'k.mode' (uses the kernel density mode from signal vs time data as the baseline), 'basal.states' (uses the median of particles' minimum state signals as the baseline), or 'lower.quartile' (uses the 25th percentile of state signals across all data as the baseline). DEFAULT = 'lower.quartile'.}
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
\item{selected_particle_state-calls.txt}{A table containing automatically called state vs time information for each selected particle.}
\item{selected_particle_residence-calls.txt}{A table containing the automatically called sequential binary binding states and residence times for each selected particle.}
\item{selected_particle_dwell-calls.txt}{A table containing the automatically called sequential multi-binding states and residence times for each selected particle.}
}
