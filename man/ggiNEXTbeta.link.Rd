% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iNEXTlink.R
\name{ggiNEXTbeta.link}
\alias{ggiNEXTbeta.link}
\title{ggplot2 extension for outcome from \code{iNEXT.beta.link}}
\usage{
ggiNEXTbeta.link(
  outcome,
  type = c("B", "D"),
  diversity = "TD",
  scale = "free",
  main = NULL,
  transp = 0.4,
  stript.size = 11,
  text.size = 13
)
}
\arguments{
\item{outcome}{the outcome from \code{"iNEXTbeta.link"}}

\item{type}{selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ;
\code{type = 'D'} for plotting 4 turnover dissimilarities.}

\item{scale}{Are scales shared across all facets (\code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})? Default is \code{"free"}.}

\item{main}{The title of the plot.}

\item{transp}{a value between 0 and 1 controlling transparency. \code{transp = 0} is completely transparent, default is 0.4.}

\item{measurement}{character indicating the label of y-axis.}
}
\value{
a figure for Beta diversity or dissimilarity diversity.

a ggplot2 object
}
\description{
\code{ggiNEXTbeta.link}: ggplot for Interpolation and extrapolation of Beta diversity with order q
}
\examples{
\dontrun{
data(Norfolk)
beta1 = iNEXTbeta.link(networks = puerto.rico$data, level = seq(0.5, 1, 0.5), datatype='abundance',q = c(0, 1, 2),
                       diversity = 'TD', nboot = 10, conf = 0.95)
beta2 = iNEXTbeta.link(networks = puerto.rico$data, level = seq(0.5, 1, 0.5), datatype='abundance',q = c(0, 1, 2),
                       diversity = 'PD', nboot = 10, conf = 0.95,
                       row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)

ggiNEXTbeta.link(beta1,diversity = 'TD', type = 'B')
ggiNEXTbeta.link(beta1,diversity = 'TD', type = 'D')
ggiNEXTbeta.link(beta2,diversity = 'PD', type = 'B')
ggiNEXTbeta.link(beta2,diversity = 'PD', type = 'D')
}
}
