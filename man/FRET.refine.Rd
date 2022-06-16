\name{FRET.refine}
\alias{FRET.refine}
\title{Dual-Camera Particle Refinement and Analysis}
\description{
Function to analyze binding states of particles from single molecule biophysics microscopy experiments with a two-camera setup.
}
\usage{
FRET.refine(path.to.file='./',file.name='Initial-Particle_Data.RData',auto.filter='none',spot.types='all',classification.strategy='classic',background.subtraction='lower.quartile')
}
\arguments{
\item{path.to.file}{A single character string indicating the path to the directory where the import file(s) are located. Write-out files will also be exported to this directory. DEFAULT = (working directory)}
\item{file.name}{A single character string indicating the name of the file that contains the selected particle traces to be analyzed. DEFAULT = 'Initial-Particle_Data.RData'}
\item{auto.filter}{A single character string, defining which particle traces should be automatically filtered out. Can be 'unbound', to filter out particle traces that remain in an unbound state the entire time, or 'all.stable', to filter out particle traces that remain in ANY single state the entire time. NOTE -- this functionality is applied to the Cy3-camera and Cy5-camera traces independently, and only particles that meet this criteria in both are filtered out. DEFAULT = 'none'.}
\item{spot.types}{A single character string, designating which types of identified particles to analyze; can be 'Cy3' to selectively analyze particles identified only on the Cy3-camera images, 'Cy5' to selectively analyze particles identified only on the Cy5-camera images, or 'DUAL' to selectively analyze particles identified on both camera images. DEFAULT = 'all'.}
\item{classification.strategy}{A single character string, designating the approach used for binding state classification from signal vs time data; can be 'classic', 'uni.dbscan', or 'var.shift'. DEFAULT = 'classic'.}
\item{background.subtraction}{A single character string, designating the approach used for subtracting background from the signal vs time plots to set unbound state averages to ~0 signal; can be 'k.mode' (uses the kernel density mode from signal vs time data as the baseline), 'basal.states' (uses the median of particles' minimum state signals as the baseline), or 'lower.quartile' (uses the 25th percentile of state signals across all data as the baseline). DEFAULT = 'lower.quartile'.}
}
\value{
Returns a list with the following components, with Cy3-camera and Cy5-camera counterparts where applicable:  image.avg = a numerical matrix containing data for a composite image of parent microscopy images, residence.data = a data frame describing manually identified binding states for selected particles, including columns for particle id, start time, end time, and residence time for each state, refined.particle.traces = a matrix containing signal vs time data for each selected particle, refined.particle.trace.rolls = a matrix containing smoothed signal vs time data for each selected particle, refined.spots = a data frame containing indices and volume information for all selected particles, refined.particle.snaps = an array providing image data for selected particles across time frames, state.calls = a numerical matrix containing automatically called state vs time information for all initially identified particles, residence.calls = a data frame containing the automatically called sequential binary binding states and residence times for each initially identified particle, & dwell.calls = a data frame containing the automatically called sequential multi-binding states and residence times for each initially identified particle.

Files:

\item{Refined-Particle_Data.RData}{An RData file containing all necessary results.}
\item{residence_data.txt}{A table describing manually identified binding states for selected particles, including columns for particle id, start time, end time, residence time, and signal average for each state.}
\item{selected_Cy3_traces.txt}{A table containing signal vs time data from the Cy3-camera images for each selected particle.}
\item{selected_Cy3_smoothed-traces.txt}{A table containing smoothed signal vs time data from the Cy3-camera images for each selected particle.}
\item{selected_Cy3-particle_summary.txt}{A table containing indices and volume information from the Cy3-camera images for all selected particle composites.}
\item{all-particle_Cy3-state-calls.txt}{A table containing automatically called state vs time information from the Cy3-camera images for all initially identified particles.}
\item{all-particle_Cy3-dwell-calls.txt}{A table containing the automatically called sequential multi-binding states and residence times from the Cy3-camera images for each initially identified particle.}
\item{all-particle_Cy3-residence-calls.txt}{A table containing the automatically called sequential binary binding states and residence times from the Cy3-camera images for each initially identified particle.}
\item{selected_Cy5_traces.txt}{A table containing signal vs time data from the Cy5-camera images for each selected particle.}
\item{selected_Cy5_smoothed-traces.txt}{A table containing smoothed signal vs time data from the Cy5-camera images for each selected particle.}
\item{selected_Cy5-particle_summary.txt}{A table containing indices and volume information from the Cy5-camera images for all selected particle composites.}
\item{all-particle_Cy5-state-calls.txt}{A table containing automatically called state vs time information from the Cy5-camera images for all initially identified particles.}
\item{all-particle_Cy5-dwell-calls.txt}{A table containing the automatically called sequential multi-binding states and residence times from the Cy5-camera images for each initially identified particle.}
\item{all-particle_Cy5-residence-calls.txt}{A table containing the automatically called sequential binary binding states and residence times from the Cy5-camera images for each initially identified particle.}
}
