\name{corner_coords}
\alias{corner_coords}
\title{Get the corner coordinates}
\usage{
  corner_coords(tr, coords_fun = "plot_phylo3_nodecoords",
    tmplocation = "manual")
}
\arguments{
  \item{tr}{A tree object in \code{\link[ape]{phylo}}
  format.}

  \item{coords_fun}{The name of the function to use to get
  node coordinates. Default: "plot_phylo3_nodecoords".}

  \item{tmplocation}{Default is "manual", which looks for
  whatever is hard-coded. Most users should probably use
  the system.file command in the examples, below.}
}
\value{
  \code{corners_list}
}
\description{
  Gets the coordinates of the corners when the tree is
  plotted.
}
\details{
  Because this function needs to use a modified version of
  the APE plot.phylo function, and for complex reasons
  APE's .C functions cannot be used elsewhere without
  causing problems with R CMD check, this function is left
  up to user specification.  Basically, the user puts in
  the name of the function, which is available in the
  extension data (extdata) directory of the package.  The
  defaults work on the developer's machine, other users may
  have to e.g. change "manual" to "auto".
}
\note{
  Go BEARS!
}
\examples{
blah=1

# Set location like this if you don't have plot_phylo3_nodecoords
# hardcoded/sourced elsehwhere
# tmplocation = system.file("extdata/", package="BioGeoBEARS")
#
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
}

