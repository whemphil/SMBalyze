\name{blind.input}
\alias{blind.input}
\title{Input File Blinding}
\description{
A function to randomize input tiff files for subsequent blinded analyses. 
}
\usage{
blind.input(path.to.files='./')
}
\arguments{
\item{path.to.files}{Path to the directory where the tiff file(s) to be imported and blinded are located. Write-out files will also be exported to this directory. DEFAULT = (working directory).}
}
\value{
Files:

\item{*.tif}{Outputs a tiff file for each original tiff file, using an alphabetical index system.}
\item{file_key.txt}{A text file matching the blinded alphabetically indexed tiff files to the original tiff files.}
\item{file_visibles.txt}{A text file giving all }
}
\note{
Input tiff files should be named [visibles]_[invisibles].tif, where [visibles] is information necessary for analysis to be provided for each blinded file and [invisibles] is the additional information to be hidden during the file blinding; the [visibles] and [invisibles] strings should not contain the '_' character.
}


