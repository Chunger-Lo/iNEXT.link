% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iNEXTlink.R
\name{ggSpec.link}
\alias{ggSpec.link}
\title{ggplot for Evenness
\code{ggSpec.link} The figure for estimation of Evenness with order q\cr}
\usage{
ggSpec.link(outcome)
}
\arguments{
\item{outcome}{a table generated from \code{Spec.link} function\cr}
}
\value{
a figure of estimated sample completeness with order q\cr
}
\description{
ggplot for Evenness
\code{ggSpec.link} The figure for estimation of Evenness with order q\cr
}
\examples{
\dontrun{
data(Norfolk)
Est <- Spec.link(data = Norfolk, diversity = 'TD', datatype = "abundance", q = c(0,1,2),
                 nboot = 30, method = "Estimated")
Emp <- Spec.link(data = Norfolk, diversity = 'TD', datatype = "abundance", q = c(0,1,2),
                 nboot = 30, method = "Empirical")
ggSpec(output = Est)
ggSpec(output = Emp)
}
}
\references{
Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
}
