\name{panoramic_spatialstats}
\alias{panoramic_spatialstats}
\title{Compute Pairwise Spatial Statistics for PANORAMIC}
\usage{
panoramic_spatialstats(
  prep,
  pairs = "auto",
  radii_um,
  stat = "local_comp_enrichment",
  nsim = 100,
  correction = "translate",
  seed = 123,
  boot = c("approx", "block"),
  tile_size = NULL,
  nx = NULL,
  ny = NULL,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
}
\arguments{
\item{prep}{List of prepared SpatialExperiment objects from \code{panoramic_prepare()}.}

\item{pairs}{Either \code{"auto"} or a data.frame with columns \code{ct1}, \code{ct2}.}

\item{radii_um}{Numeric vector of radii (microns).}

\item{stat}{Summary statistic. The default \code{"local_comp_enrichment"}
uses an edge-corrected CSR null. \code{"local_comp_global_enrichment"} is
an opt-in global-composition (random-label) alternative conditional on the
observed cell locations. Other supported values are \code{"Lcross"},
\code{"Kcross"}, \code{"Lest"}, and \code{"Kest"}.}

\item{nsim}{Integer number of Loh bootstrap replicates per sample/pair.}

\item{correction}{Edge correction method passed to spatstat summaries. It is
ignored for local-composition statistics.}

\item{seed}{Optional seed for reproducible bootstrap sampling.}

\item{boot}{Bootstrap mode: \code{"approx"} or \code{"block"}.}

\item{tile_size}{Optional tile size for block bootstrap in coordinate units.}

\item{nx, ny}{Optional tile counts for block bootstrap when \code{boot = "block"}.
When neither is supplied, a 4-by-4 grid is used.}

\item{BPPARAM}{BiocParallel backend used across samples.}

\item{verbose}{If \code{TRUE}, show bootstrap console output.}
}
\value{
A \code{SummarizedExperiment} with assays:
\itemize{
\item \code{yi}: centered estimates per \code{(ct1, ct2, radius_um)} feature and sample.
\item \code{vi}: variance estimates aligned to \code{yi}.
}
\code{rowData} contains \code{ct1}, \code{ct2}, \code{radius_um}, \code{stat};
\code{colData} contains \code{sample} and \code{group}.
}
\description{
Compute pairwise spatial summary curves and bootstrap variances across samples
for all requested cell-type pairs and radii.
}
\details{
For L-function statistics, PANORAMIC computes variances via the delta method
from corresponding K-function bootstrap estimates and centers by \code{r}.

For \code{stat = "local_comp_enrichment"}, PANORAMIC reports the
edge-corrected observed-minus-expected local target composition in percentage
points. The expected proportion is the expected target count in the overlap
area of the radius disc and the tissue window, divided by the observed local
neighbor count. For self-pairs, the anchor is excluded from the expected target
density.

\code{stat = "local_comp_global_enrichment"} instead subtracts the
sample-wide target proportion (excluding the anchor for self-pairs). It is a
random-label null conditional on the observed cell locations, and can be useful
when the question concerns label mixing within an observed tissue architecture.
It does not by itself adjust for spatially varying abundance or establish a
direct cell-cell interaction. Both local-composition statistics use a spatial
bootstrap for uncertainty.
}
\examples{
toy <- panoramic_simulate_dataset(seed = 1)
prep <- panoramic_prepare(
  spe_list = toy$spe_list,
  design = toy$design,
  cell_type = "cell_type",
  min_cells = 2,
  window = "rect",
  BPPARAM = BiocParallel::SerialParam()
)
se_stats <- panoramic_spatialstats(
  prep = prep,
  radii_um = c(10, 20),
  nsim = 5,
  correction = "translate",
  seed = 1,
  BPPARAM = BiocParallel::SerialParam()
)
se_stats
}
