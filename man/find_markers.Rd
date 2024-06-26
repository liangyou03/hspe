% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_markers_fn.R
\name{find_markers}
\alias{find_markers}
\title{Find marker genes for each cell type.}
\usage{
find_markers(
  Y,
  references = NULL,
  pure_samples = NULL,
  marker_method = "ratio"
)
}
\arguments{
\item{Y}{Expression matrix.

(Required) Two-dimensional numeric. Must implement \code{as.matrix}.

Each row contains expression measurements for a particular sample. Each columm contains the measurements of the same gene over all individuals. Can either contain just the mixture samples to be deconvolved or both the mixture samples and the reference samples. See \code{pure_samples} and \code{references} for more details.}

\item{references}{Cell-type reference expression matrix.

(Optional) Two-dimensional numeric. Must implement \code{as.matrix}. Must have same number of columns as \code{Y}. Columns must correspond to columns of \code{Y}.

Each row contains expression measurements for a reference profile of a particular cell type. Columns contain measurements of reference profiles of a gene. Optionally may merge this matrix with \code{Y} and use \code{pure_samples} to indicate which rows of \code{Y} are pure samples. If \code{pure_samples} is not specified \code{references} must be specified. In this case each row of \code{references} is assumed to be a distinct cell-type. If both \code{pure_samples} and \code{references} are specified then \code{pure_samples} specifies to which cell-type each row of \code{references} corresponds.}

\item{pure_samples}{The pure sample indicies.

(Optional) List of one-dimensional integer. Must implement \code{as.list}.

The i-th element of the top-level list is a vector of indicies (rows of \code{Y} or \code{references}) that are pure samples of type i. If \code{references} is not specified then this argument identifies which rows of \code{Y} correspond to pure reference samples of which cell-types. If \code{references} is specified then this makes same idenficiation but for the \code{references} matrix instead.}

\item{marker_method}{Method used to rank marker genes.

(Optional) One-dimensional string.

The method used to rank genes as markers. If not supplied defaults to ``ratio''. Only used if markers are not provided to argument ``markers''. Options are
\itemize{
\item{'ratio'}{ selects and ranks markers by the ratio of the mean expression of each gene in each cell type to the mean of that gene in all other cell types.}
\item{'regression '}{ selects and ranks markers by estimated regression coefficients in a series of regressions with single covariate that is indicator of each type.}
\item{'diff'}{ selects and ranks markers based upon the difference, for each cell type, between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
\item{'p.value'}{ selects and ranks markers based upon the p-value of a t-test between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
}}
}
\value{
List with four elements. ``L'' is respective ranked markers for each cell type and ``V'' is the corresponding values of the ranking method (higher are better) used to determine markers and sort them, ``M'' is the matrix used to create the other two arguments after sorting and subsetting, and ``sM'' is a sorted version of M.
}
\description{
Find marker genes for each cell type.
}
\examples{
truth = shen_orr_ex$annotation$mixture
pure_samples <- lapply(1:3, function(i) {
   which(truth[, i] == 1)
})
Y <- shen_orr_ex$data$log
find_markers(Y=Y,pure_samples=pure_samples,
,marker_method='ratio')
}
