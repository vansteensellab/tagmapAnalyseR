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
#' @param mode Readout mode, eiter 'counts' for number of insertions per sliding
#' window or 'fraction' for number of insertions per sliding across experimental
#' conditions.
#' @param numerator Experimental condition(s) for numerator of fraction of
#' insertions (only required for fraction readout).
#' @param denominator Experimental condition(s) for denominator of fraction of
#' insertions (only required for fraction readout).
#'
#' @return Returns data table for putative insertion sites outputted by the
#' TagMap pipeline.
#' \describe{
#' \item{windowStart, windowEnd}{Chromosome start/end positions of a specific
#' window}
#' \item{readout}{Readout of the sliding window, e.g. number of integrations}
#' }
#'
#' @importFrom tibble tibble
#' @importFrom dplyr %>%
#' @importFrom utils tail
#' @export
chromosomeSlidingWindow <- function(data, chromosome, windowStart, windowEnd, windowSize = 10000, stepSize = 1000, mode, numerator = NULL, denominator = NULL) {
  # Test readout mode formatting
  if (!tolower(mode) %in% c('counts', 'fraction')){
    stop("Readout mode incorrect, use either 'counts' or 'fraction'.")
  }

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
    # Readout: Number of insertions
    if(tolower(mode) == 'counts'){
      windowReadout[i] <- slidingWindowReadoutCounts(
        data = data,
        chromosome = chromosome,
        windowStart = windows[i],
        windowEnd = (windows[i] + windowSize)
      )
    }
    # Readout: Fraction of insertions between experimental conditions
    if(tolower(mode) == 'fraction'){
      windowReadout[i] <- slidingWindowReadoutCounts(
        data = data,
        chromosome = chromosome,
        windowStart = windows[i],
        windowEnd = (windows[i] + windowSize),
        numerator = numerator,
        denominator = denominator
      )
    }
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
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param data Insertion site data (as formatted by TagMap pipeline), columns
#' should minimally include- and be named 'chr', 'start' and 'end.
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



#' Calculate the fraction of insertions between experimental conditions
#' within a sliding window
#'
#' Subsets the complete data object to only contain insertion site data within
#' the genomic coordinates of the sliding window. Next, insertions are divided
#' by experimental conditions ('class' column) into the numerator insertions and
#' denominator insertions (e.g. effect vs. control, respectively). The fraction
#' of insertions is calculated by dividing the numerator by the denominator, a
#' pseudocount of 1 is added for zero-values of the denominator.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param data Insertion site data (as formatted by TagMap pipeline), columns
#' should minimally include and be formatted as 'chr', 'start', 'end' and
#' 'class'.
#' @param chromosome Chromosome, formatted as "chr(number)".
#' @param windowStart Window start position.
#' @param windowEnd Window end position.
#' @param numerator Names of experimental condition(s) in column 'class' of data
#' to select numerator insertions for the fraction.
#' @param denominator of experimental condition(s) in column 'class' of
#' data to select denominator insertions for the fraction.
#'
#' @return Returns the fraction of insertions between experimental conditions.
#'
#' @importFrom dplyr filter %>%
#' @export
slidingWindowReadoutFraction <- function(data, chromosome, windowStart, windowEnd, numerator, denominator, pseudocount = 1) {
  # Select subsets of numerator and denominator insertions
  insertions <- data %>% dplyr::filter(chr == chromosome & start >= windowStart & end <= windowEnd)
  insertionsNumerator <- nrow(insertions %>% dplyr::filter(class %in% numerator))
  insertionsDenominator <- nrow(insertions %>% dplyr::filter(class %in% denominator))

  # Add pseudocount of 1 to zero-values of the denominator
  if(insertionsDenominator == 0){
    insertionsDenominator <- insertionsDenominator + pseudocount
  }

  # Return fraction
  return(insertionsNumerator / insertionsDenominator)
}
