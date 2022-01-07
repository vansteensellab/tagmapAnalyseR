#' Run sliding window over chromosome for transposon readouts
#'
#' Given genomic coordinates (chromosome, start, end) a sliding window runs over
#' the genome with a specific window size and step size. Readouts of this window
#' are calculated with specific "readout" functions, e.g.
#' \code{\link{slidingWindowReadoutCounts}} returns the absolute number of
#' insertions within a window, however possibilities are countless.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param data Insertion site data (as formatted by TagMap pipeline).
#' @param chromosome Chromosome, formatted as "chr(number)".
#' @param windowStart Chromosome position (integer) to calculate window
#' coordinates from.
#' @param windowEnd Chromosome position (integer) to calculate window
#' coordinates to.
#' @param windowSize Size of each window, e.g. 1e3 or 1e4 bp long.
#' @param stepSize Size that each step increases, from minimal steps of 1 bp to
#' any other custom value.
#'
#' @return Returns data table for putative insertion sites outputted by the
#' TagMap pipeline.
#' \describe{
#' \item{windowStart, windowEnd}{Chromosome start/end positions of a specific
#' window}
#' \item{readout}{Readout of the window, such as number of integrations}
#' }
#'
#' @importFrom tibble tibble
#' @importFrom dplyr %>%
#' @importFrom utils tail
#' @export
chromosomeSlidingWindow <- function(data, chromosome, windowStart, windowEnd, windowSize, stepSize) {
  # Calculate sliding window start positions.
  windows <- seq(
    from = windowStart,
    to = (windowEnd - windowSize),
    by = stepSize
  )

  # Check if modulus of last window end - window size over the step size is
  # non-zero, add a final window start point to avoid the windows ending
  # prematurely.
  if (((windowEnd - windowSize) %% stepSize) != 0) {
    windows <- c(windows, (tail(windows, 1) + windowSize))
  }

  windowReadout <- vector(length = length(windows))
  # Iterate windows to get readout specified by custom function.
  for (i in 1:length(windows)) {
    windowReadout[i] <- slidingWindowReadoutCounts(
      data = data,
      chromosome = chromosome,
      windowStart = windows[i],
      windowEnd = (windows[i] + windowSize)
    )
  }

  return(tibble::tibble(
    windowStart = windows,
    windowEnd = windows + windowSize,
    readout = windowReadout
  ))
}




#' Calculate the number of insertions within a sliding window
#'
#' Subsets the complete data object to only contain insertion site data within
#' the genomic coordinates of the sliding window.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param data Insertion site data (as formatted by TagMap pipeline).
#' @param chromosome Chromosome, formatted as "chr(number)".
#' @param windowStart Window start position.
#' @param windowEnd Window end position.
#'
#' @return Returns number of insertions in the subsetted data object.
#'
#' @importFrom dplyr filter %>%
#' @export
slidingWindowReadoutCounts <- function(data, chromosome, windowStart, windowEnd) {
  return(nrow(data %>% dplyr::filter(chr == chromosome & start >= windowStart & end <= windowEnd)))
}
