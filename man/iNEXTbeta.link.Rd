% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iNEXTlink.R
\name{iNEXTbeta.link}
\alias{iNEXTbeta.link}
\title{Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity
Function \code{iNEXTbeta.link} Interpolation and extrapolation of Beta diversity with order q}
\usage{
iNEXTbeta.link(
  data,
  diversity = "TD",
  level = seq(0.5, 1, 0.5),
  datatype = c("abundance", "incidence_raw"),
  q = c(0, 1, 2),
  nboot = 20,
  conf = 0.95,
  row.tree = NULL,
  col.tree = NULL
)
}
\arguments{
\item{data}{a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
and thus N matrices, users must first merge the N matrices by species identity to obtain a large
merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
data. The row (species) names of data must match the species names in the phylogenetic tree and
thus cannot be missing.}

\item{diversity}{a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic',
and 'FD' = 'Functional' under certain threshold.}

\item{datatype}{data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).}

\item{q}{a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.}

\item{nboot}{a positive integer specifying the number of bootstrap replications when assessing
sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50}

\item{conf}{a positive number < 1 specifying the level of confidence interval. Default is 0.95.}

\item{type}{desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.}
}
\value{
A list of seven lists with three-diversity and four-dissimilarity.
}
\description{
Interpolation (rarefaction) and extrapolation of Chao et al.’s (2021) network diversity and mean network diversity
Function \code{iNEXTbeta.link} Interpolation and extrapolation of Beta diversity with order q
}
\examples{
\dontrun{
# example
data(puerto.rico)
beta1 = iNEXTbeta.link(data = puerto.rico$data, level = seq(0.5, 0.9, 0.4), datatype='abundance',q = c(0, 1, 2),
                       diversity = 'TD', nboot = 10, conf = 0.95)
beta2 = iNEXTbeta.link(networks = puerto.rico$data, level = seq(0.5, 0.9, 0.4), datatype='abundance',q = c(0, 1, 2),
                       data = 'PD', nboot = 10, conf = 0.95,
                       row.tree = puerto.rico$row.tree, col.tree = puerto.rico$col.tree)
}
}
\references{
Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J.(2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters 8, 148-159. (pdf file) Spanish translation in pp. 85-96 of Halffter, G. Soberon, J., Koleff, P. and Melic, A. (eds) 2005 Sobre Diversidad Biologica: el Sognificado de las Diversidades Alfa, Beta y Gamma. m3m-Monografias 3ercer Milenio, vol. 4, SEA, CONABIO, Grupo DIVERSITAS & CONACYT, Zaragoza. IV +242 pp.
Chiu, C.-H., Jost, L. and Chao*, A. (2014). Phylogenetic beta diversity, similarity, and differentiation measures based on Hill numbers. Ecological Monographs 84, 21-44.
}
