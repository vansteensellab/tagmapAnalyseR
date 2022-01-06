#' Read putative insertion sites called by pipeline.
#'
#' Reads putative insertions called by the pipeline into a data.table, according to the
#' specific column naming convention, to be used for further analyses.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param inputFile Path to input file.
#'
#' @return Returns data table for putative insertion sites outputted by the
#' TagMap pipeline.
#' \describe{
#' \item{read_name}{Read name}
#' \item{chr}{Read chromosome}
#' \item{read_start}{Read start position}
#' \item{read_end}{Read end position}
#' \item{mapq}{Read mapping quality (mapq)}
#' \item{strand}{Read strandedness}
#' \item{TIS_seq}{Insertion site overhang sequence}
#' \item{TIS_start}{Insertion site start position}
#' \item{TIS_end}{Insertion site end position}
#' \item{sample}{Experiment sample name, e.g. 'WT78_SB_Fw'}
#' }
#'
#' @importFrom data.table data.table fread
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr replace_na
#' @export
readPutativeInsertions <- function(inputFile) {
  return(data.table::fread(inputFile,
    header = F,
    sep = "\t",
    col.names = c(
      "read_name", "chr", "read_start",
      "read_end", "strand", "mapq",
      "TIS_seq", "TIS_start", "TIS_end",
      "reaction"
    )
  ) %>%
    dplyr::mutate(TIS_seq = tidyr::replace_na(TIS_seq, "NA")))
}



#' Find transposon insertion sites with ambiguous positions
#'
#' Identify the reads for insertion sites where the exact insert position is
#' ambiguous due to locally shifted read positions.
#'
#' Reads for which the previous or next read on the genome has shifted by
#' \emph{up to N bp} are detected by this function. Such reads should indicate
#' that despite the local shift, reads map to the identical insertion site.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param dt A data.table as returned by readPutativeInsertions.
#' @param padding Padding added around the TIS position that is used to evaluate
#' whether the next upstream or downstream read overlaps within a range of N bp.
#'
#' @return Returns read values for shifted insertions:
#' \describe{
#' \item{read_name}{Read name}
#' \item{chr}{Chromosome}
#' \item{TIS_start}{Insertion site start position}
#' \item{TIS_end}{Insertion site end position}
#' \item{strand}{Read strandedness}
#' \item{TIS_seq}{Insertion site overhang sequence}
#' }
#'
#' @importFrom dplyr arrange filter select %>%
#' @importFrom data.table data.table
#' @export
findAmbiguousInsertionSites <- function(dt, padding = 2) {
  dt <- dplyr::arrange(dt, dt$chr, dt$TIS_start)

  ambiguousInsertions <- data.table::as.data.table(dplyr::filter(
    dt,
    (dplyr::lead(dt$TIS_start) - dt$TIS_start >= 1 &
      dplyr::lead(dt$TIS_start) - dt$TIS_start <= padding) |
      (dt$TIS_start - dplyr::lag(dt$TIS_start) >= 1 &
        dt$TIS_start - dplyr::lag(dt$TIS_start) <= padding)
  ) %>%
    dplyr::select(read_name, chr, TIS_start, TIS_end, strand, TIS_seq))

  return(ambiguousInsertions)
}



#' Update insertion site data
#'
#' Update data for insertion sites, including the exact inserted sequence
#' ("overhang") when read data is ambiguous and the exact position when read
#' data displays neighboring yet shifted reads (returned by
#' \code{\link{findAmbiguousInsertionSites}}).
#'
#' Insertion site data are updated using a decision tree.
#'
#' If argument \code{ambiguousInsertions} is passed, the ambiguous insertion
#' site positions are investigated: the position with highest number of reads is
#' returned if possible, otherwise the overhang sequences are compared to the
#' one expected from argument \code{overhang}. If all ambiguous positions match
#' the expected overhang, a random position will be returned, otherwise it is
#' checked if one position uniquely matches the overhang which will be returned
#' or that any positions matches the overhang at all, when none are present then
#' a random position will be returned and if more positions uniquely match the
#' overhang then a random position will be returned.
#'
#' If no argument \code{ambiguousInsertions} is passed, insertion site overhang
#' sequences are investigated: if only one unique sequence is present, then this
#' will be returned, otherwise the sequence with the highest number of reads is
#' returned if possible and otherwise a random sequence is returned.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param dt A data table with either all or subsets of insertion site data, as
#' processed in the function mapInsertionSites.
#' @param gr A GRanges object with detailed paired-end read data underlying the
#' insertion sites in dt.
#' @param overhang Insertion site overhang sequence, e.g. 'TTAA' for PiggyBac.
#' @param readsPerTIS A data table of read counts and genomic coordinates per
#' unique insertion site, typically generated in the function mapInsertionSites.
#' @param ambiguousInsertions An optional data table with shifted insertion site
#' locations outputted by the function \code{\link{findAmbiguousInsertionSites}}
#' . When left out, only overhang sequence values are updated.
#'
#' @return Returns the input data table \code{dt} with updated values.
#'
#' @importFrom dplyr filter select mutate %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table data.table as.data.table
#' @export
updateInsertionSiteData <- function(dt, gr, overhang, readsPerTIS, ambiguousInsertions = NULL) {
  colnames(dt)[which(colnames(dt) == "seqnames")] <- "chr"

  # Iterate over insertion sites to update data for.
  for (i in 1:nrow(dt)) {
    # Extract for the insertion site at index i: overhang sequence(s), genomic
    # location and counts of supporting reads.

    if (is.null(ambiguousInsertions)) {
      # Update non-unique insertion site sequences (no shifted sites given).

      indexedData <- data.table::data.table(
        chr = dt[i, "chr"],
        start = dt[i, "start"],
        end = dt[i, "end"],
        table(gr[unlist(dt[i, "revmap"]), ]$TIS_seq)
      )
      colnames(indexedData) <- c("chr", "start", "end", "TIS_seq", "count")
      indexedData <- indexedData[, c("TIS_seq", "chr", "start", "end", "count")]
    } else if (!is.null(ambiguousInsertions)) {
      # Update non-unique insertion site sequences and shifted insertion sites.

      indexedData <- unique(
        ambiguousInsertions %>%
          dplyr::filter(read_name %in% gr[unlist(dt[i, "revmap"]), ]$read_name) %>%
          dplyr::select(TIS_seq, chr, TIS_start, TIS_end)
      )
      colnames(indexedData) <- c("TIS_seq", "chr", "start", "end")

      count <- c(1:nrow(indexedData))
      for (c in 1:nrow(indexedData)) {
        count[c] <- readsPerTIS[TIS_start == indexedData[c,]$start & TIS_end == indexedData[c,]$end]$n
      }
      indexedData$count <- count
    }
    print(indexedData)


    # Determine insertion site updates with the following decision tree:
    # IF shifted insertion site data is absent
    #   YES:
    #     IF only one insertion site sequence is present
    #       YES: return this single insertion site with sequence
    #       NO:
    #         IF all have the equal read counts:
    #           YES: select insertion site sequence at random
    #           NO:  select insertion site sequence with maximum read count
    #   NO:
    #     IF all shifted positions have equal read counts:
    #       YES:
    #         IF all shifted positions match the overhang sequence:
    #           YES: randomly select one insertion site position
    #           NO:
    #             IF only one position matches the overhang:
    #               YES: select the insertion site position matching overhang
    #               NO:
    #                  IF any of the reads matches the overhang:
    #                     YES: select insertion site position at random
    #                     NO: select insertion site position matching overhang at random
    #       NO: select insertion site position with maximum read count
    updated <- data.table::as.data.table(ifelse(
      is.null(ambiguousInsertions),
      yes = ifelse(nrow(indexedData) == 1,
                   yes = list(indexedData),
                   no = ifelse(length(unique(indexedData$count)) == 1,
                               yes = list(indexedData[round(stats::runif(1, 1, nrow(indexedData)))]),
                               no = list(indexedData[which.max(count), ])
                   )
      ),
      no = ifelse(length(unique(indexedData$count)) == 1,
                  yes = ifelse(all(unique(indexedData$TIS_seq) == overhang),
                               yes = list(indexedData[round(stats::runif(1, 1, nrow(indexedData))), ]),
                               no = ifelse(nrow(indexedData[TIS_seq == overhang]) == 1,
                                           yes = list(indexedData[TIS_seq == overhang]),
                                           no = ifelse(nrow(indexedData[TIS_seq == overhang]) == 0,
                                                       yes = list(indexedData[round(stats::runif(1, 1, nrow(indexedData)))]),
                                                       no = list(indexedData[TIS_seq == overhang][round(stats::runif(1, 1, nrow(indexedData[TIS_seq == overhang])))])
                                                       )
                               )
                  ),
                  no = list(indexedData[which.max(count), ])
      )
    ))


    # Apply updates to insertion site sequence and/or shifted position.
    if (is.null(ambiguousInsertions)) {
      dt[i, "TIS_seq"] <- updated$TIS_seq
    } else if (!is.null(ambiguousInsertions)) {
      dt[i, c("start", "end", "TIS_seq")] <- updated[, c("start", "end", "TIS_seq")]
      dt[i, "width"] <- length(updated$start:updated$end)
    }
  }
  return(dt)
}




#' Read BAM file(s)
#'
#' Function to read one or two BAM files with a Samtools command line call. The
#' data of two BAM files are combined into a single output.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param bamFilePath Full path(s) to BAM file(s).
#' @param samtoolsPath Full path to SAMTOOLS.
#' @param flags BAM flag(s) to filter reads by.
#'
#' @return Returns data frame with SAM-formatted data (see specification at
#' at \url{https://samtools.github.io/hts-specs/SAMv1.pdf}, however the default
#' QUAL column is removed here.
#'
#' @importFrom utils read.table
#' @export
readBam <- function(bamFilePath, samtoolsPath, flags = 0) {
  # Evaluate if only 1 or 2 BAM files are supplied
  if (!length(bamFilePath) %in% c(1, 2)) {
    stop("Only 1 or 2 BAM files allowed!")
  }
  samFormat <- c(
    "QNAME", "FLAG", "RNAME", "POS", "MAPQ",
    "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ"
  )

  # Read either 1 or 2 BAM files
  if (length(bamFilePath) == 1) {
    cmd <- system(paste(samtoolsPath, " view -f ", flags, " ", bamFilePath[1],
                        " | cut -f-10", sep = ""),
                  intern = TRUE)
  } else if (length(bamFilePath) == 2) {
    cmd <- c(system(paste(samtoolsPath, " view -f ", flags, " ", bamFilePath[1],
                          " | cut -f-10", sep = ""),
                    intern = TRUE),
             system(paste(samtoolsPath, " view -f ", flags, " ", bamFilePath[2],
                          " | cut -f-10", sep = ""),
                    intern = TRUE)
             )
  }
  bam <- utils::read.table(text = cmd, sep = "\t", col.names = samFormat)

  # Return formatted data
  return(bam)
}




#' Switch read strand
#'
#' Function to switch read strand when calling the Tagmentation reaction strand.
#' Reverse strands are switched to forward and vice versa, the output is
#' directly returned.
#'
#' This function is used for \code{\link{mapInsertionSites}} where the direction
#' of one Tagmentation reaction can be determined from the reads of the opposite
#' reaciton
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param strand Character of read strand, either '+' or '-'.
#'
#' @return Returns character '+' or '-'.
#' @export
switchReadStrand <- function(strand){
  switched <- ifelse(strand == '+', yes = '-', no = '+')
  return(switched)
}
