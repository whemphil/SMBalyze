\name{calc.kn1}
\alias{calc.kn1}
\title{Dissociation Rate Calculation}
\description{
A function to calculate dissociation rate from residence time data.
}
\usage{
calc.kn1(path.to.file='./',file.name='Refined-Particle_Data.RData',use.auto.times='n',min.residence=0,max.residence=Inf)
}
\arguments{
\item{path.to.file}{Path to the directory where the input file(s) are located. Write-out files will also be exported to this directory. DEFAULT = (working directory)}
\item{file.name}{The name of the Rdata output file from the 'refine.particles' function. DEFAULT = 'Refined-Particle_Data.RData'.}
\item{use.auto.times}{Should the automatcally called residence data, as opposed to manually selected, be used for analysis (y/n)? DEFAULT = 'n'.}
\item{min.residence}{The minimum residence times to be included in analysis. DEFAULT = 0.}
\item{max.residence}{The maximum residence times to be included in analysis. DEFAULT = Inf.}
}
\value{

Files:

\item{Dissociation-Rate_Analysis.pdf}{Graphs of the residence time distribution data.}
}
