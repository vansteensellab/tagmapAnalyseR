#' Map transposon insertion sites.
#'
#' Maps transposon insertion sites (TIS) in high detail from locations called by
#' the pipeline.
#'
#' Input data table \code{dt} of filtered R2 reads (MAPQ >= 10, proper pairs, no
#' PCR duplicates) where the transposon overhang sequence and position is
#' reported, are processed in this function. TISes can be filtered based on read
#' depth, ambiguous insertion positions or overhang sequences are corrected
#' for. More detailed metadata per TIS are reported, including
#' Tagmentation-reaction specific counts, median MAPQ and strand orientation.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param dt A data.table with TIS data taken directly from
#' \code{\link{readPutativeInsertions}}.
#' @param bam Full path to BAM file(s).
#' @param overhang TIS overhang sequence, e.g. 'TTAA' for PiggyBac.
#' @param depth Minimal read depth to filter TISes by.
#' @param gapWidth Minimal permitted gap between reads for merging reads
#' into a single insertion site location.
#' @param ignoreStrand Boolean indicating whether to ignore strand direction
#' and merge reads nonetheless (TRUE) or to separate merging by strand
#' direction (FALSE)
#' @param samtoolsPath Full path to SAMTOOLS.
#' @param ambiguousInsertions A data table with ambiguous insertion site locations
#' taken directly from \code{\link{findAmbiguousInsertionSites}}.
#'
#' @return A data table with the columns:
#' \describe{
#' \item{chr, start, end, seq, strand}{Genomic coordinates, overhang sequence
#' and strand of the insertion site.}
#' \item{region_start, region_end}{Exact genomic position of reads around TIS.}
#' \item{read_names}{Names of individual reads supporting a particular TIS.}
#' \item{read_count}{Number of reads supporting a TIS.}
#' \item{mapq}{Median mapq of reads supporting a TIS.}
#' }
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame reduce as.data.frame
#' @importFrom dplyr rowwise mutate filter count select left_join coalesce %>%
#' @importFrom stats runif
#' @export
mapInsertionSites <- function(dt, bam, overhang, depth = 1, gapWidth = 1, ignoreStrand = TRUE, samtoolsPath = NULL, ambiguousInsertions = NULL) {
  # VARIABLE DECLARATION AND TESTING ----------------------------------------


    # Test for function path to samtools
  if (is.null(samtoolsPath)) {
    SAMTOOLS <- system("which samtools", intern = TRUE)
  } else if (!is.null(samtoolsPath)) {
    SAMTOOLS <- samtoolsPath
  } else {
    stop("SAMTOOLS NOT FOUND")
  }


  # Test overhang sequence
  if (is.null(overhang)) {
    stop("NO INSERTION SITE OVERHANG GIVEN")
  }


  # Test shifted insertion site position data
  if (!exists('ambiguousInsertions')) {
    ambiguousInsertions <- findAmbiguousInsertionSites(dt)
  } else if (is.null(ambiguousInsertions)){
    ambiguousInsertions <- findAmbiguousInsertionSites(dt)
  }


  # Calculate length of overhang sequence
  overhangLength <- nchar(overhang)


  # Calculate reference table of read counts per unique insertion site
  readsPerTIS <- dt %>% dplyr::count(chr, TIS_start, TIS_end)


  # Evaluate which Tagmentation reactions are in the data
  hasForwardReaction <- 'Fw' %in% dt$reaction
  hasReverseReaction <- 'Rv' %in% dt$reaction


  # Column names expected for particular single or double Tagmentation reactions
  COLS_FWRV <- c('read_count_1.x','mapq_1.x','read_count_2.x','mapq_2.x')
  COLS_FWRV_AMBIGUOUS <- c('seqnames','start','end','width',
    'TIS_seq.y','strand.x','region_start','region_end','read_names.x','revmap.x',
    'read_count.x','mapq.x','read_count_1.x','mapq_1.x','read_count_2.x',
    'mapq_2.x')
  COLS_FWRV_UNAMBIGUOUS <- c('seqnames','start.x','end.x','width.x',
    'TIS_seq','strand.x','region_start','region_end','read_names.x','revmap.x',
    'read_count.x','mapq.x','read_count_1.x','mapq_1.x','read_count_2.x',
    'mapq_2.x')
  RENAME_FWRV <- c('chr','start','end','width','TIS_seq','strand',
    'region_start', 'region_end', 'read_names','revmap','read_count','mapq',
    'read_count_1', 'mapq_1','read_count_2', 'mapq_2')
  COLS_SINGLE_AMBIGUOUS <- c('seqnames','start','end','width','TIS_seq.y',
    'strand.x','region_start','region_end','read_names.x','revmap.x',
    'read_count.x','mapq.x')
  COLS_SINGLE_UMAMBIGUOUS <- c('chr','start.x','end.x','width.x','TIS_seq',
    'strand.x','region_start','region_end','read_names.x','revmap.x','read_count.x','mapq.x')
  RENAME_SINGLE <- c('chr','start','end','width','TIS_seq','strand',
    'region_start','region_end','read_names','revmap','read_count',
    'mapq')
  AMBIGUOUS_STRAND <- '*'




  # MAP BASIC INSERTION SITE DATA -------------------------------------------


  # Create GenomicRanges object
  dtGenomicRanges <- GenomicRanges::makeGRangesFromDataFrame(
    dt,
    keep.extra.columns = TRUE,
    start.field = "TIS_start",
    end.field = "TIS_end",
    strand.field = "strand"
  )


  # Reduce to overlapping insertion sites
  dtReducedRanges <- GenomicRanges::reduce(
    dtGenomicRanges,
    min.gapwidth = gapWidth,
    ignore.strand = ignoreStrand,
    with.revmap = TRUE
  )


  # Get read counts, median MAPQ and read names, remove sites below minimum read depth
  readsR1 <- readBam(bamFilePath = bam,
                     samtoolsPath = SAMTOOLS,
                     flags = '64')
  # readsR1 = Sequence data for R1 reads in SAM-format to derive exact start
  #           and end of paired-end reads R1/R2. Flag = 64 selects R1 reads
  dtOutput <- GenomicRanges::as.data.frame(dtReducedRanges) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(read_count = length(unlist(revmap))) %>%
    dplyr::filter(read_count >= depth) %>%
    dplyr::mutate(mapq = median(dtGenomicRanges[unlist(revmap),]$mapq)) %>%
    dplyr::mutate(read_names = list(dtGenomicRanges[unlist(revmap),]$read_name))
  # dtOutput = Main object with filtered data on insertion sites


  # CRUCIAL: if dtOutput is already empty at this stage, further processing will
  # be skipped to return an empty object instead.
  if ( nrow(dtOutput) != 0 ) {
    # Get start/end of reads in region around insertion site
    dtOutput <- dtOutput %>%
      dplyr::rowwise() %>%
      dplyr::mutate(region_start = min(readsR1[readsR1$QNAME %in% read_names,c('POS','PNEXT')],
                                       dtGenomicRanges[unlist(revmap),]$read_start)) %>%
      dplyr::mutate(region_end = max(readsR1[readsR1$QNAME %in% read_names,c('POS','PNEXT','TLEN')] %>%
                                       mutate(end = min(POS,PNEXT) + abs(TLEN) - 1) %>% select(end),
                                     dtGenomicRanges[unlist(revmap), ]$read_end))


    # Add column to store insertion site overhang sequence
    dtOutput['TIS_seq'] <- NA




    # PROCESS FORWARD/REVERSE TAGMENTATION REACTION SEPARATELY ----------------


    # Get Tagmentation reaction-specific read count and median MAPQ
    if (hasForwardReaction & hasReverseReaction) {
      # Read count for forward (read_count_1) and reverse (read_count_2) reactions
      forwardReaction <- dt %>% dplyr::filter(reaction == 'Fw') %>% dplyr::pull('read_name')
      reverseReaction <- dt %>% dplyr::filter(reaction == 'Rv') %>% dplyr::pull('read_name')
      dtOutput <- dtOutput %>% dplyr::rowwise() %>%
        dplyr::mutate(read_count_1 = sum(read_names %in% forwardReaction),
                      read_count_2 = sum(read_names %in% reverseReaction))

      # Median MAPQ for forward (mapq_1) and reverse (mapq_2) reactions
      mapq_1 <- c(1:nrow(dtOutput))
      mapq_2 <- c(1:nrow(dtOutput))
      for (indexTIS in 1:nrow(dtOutput)) {
        mapq_1[indexTIS] <- median((GenomicRanges::as.data.frame(dtGenomicRanges[unlist(dtOutput[indexTIS,'revmap']),]) %>% dplyr::filter(reaction == 'Fw'))$mapq)
        mapq_2[indexTIS] <- median((GenomicRanges::as.data.frame(dtGenomicRanges[unlist(dtOutput[indexTIS,'revmap']),]) %>% dplyr::filter(reaction == 'Rv'))$mapq)
      }
      dtOutput$mapq_1 <- mapq_1
      dtOutput$mapq_2 <- mapq_2
      dtOutput <- dtOutput %>% dplyr::mutate(mapq_1 = tidyr::replace_na(mapq_1, '-'),
                                             mapq_2 = tidyr::replace_na(mapq_2, '-'))
    }




    # CORRECT FOR SHIFTED INSERTION SITES -------------------------------------


    # Find sites with ambiguous genomic positions due to conflicting reads
    dtAmbiguousPositions <- dtOutput[dtOutput$width > overhangLength, ]
    # dtAmbiguousPositions = Subset of insertion sites with a wider than expected
    #                        insertion site length


    # Correct data if shifted insertion sites are reported
    if (nrow(dtAmbiguousPositions) != 0) {
      dtAmbiguousPositions <- updateInsertionSiteData(dt = dtAmbiguousPositions,
                                                      gr = dtGenomicRanges,
                                                      overhang = overhang,
                                                      readsPerTIS = readsPerTIS,
                                                      ambiguousInsertions = ambiguousInsertions)

      # Re-join data
      dtOutput <- dplyr::left_join(dtOutput, dtAmbiguousPositions,
                                   by = c("region_start", "region_end")) %>%
        dplyr::mutate(start = dplyr::coalesce(start.y, start.x),
                      end = dplyr::coalesce(end.y, end.x),
                      width = dplyr::coalesce(width.y, width.x))

      # Check for columns indicating a double Tagmentation reaction (Fw/Rv samples)
      if ( !all(COLS_FWRV %in% colnames(dtOutput)) ) {
        dtOutput <- dtOutput %>% dplyr::select(COLS_SINGLE_AMBIGUOUS)
        colnames(dtOutput) <- RENAME_SINGLE
      } else if ( all(COLS_FWRV %in% colnames(dtOutput)) ) {
        dtOutput <- dtOutput %>% dplyr::select(COLS_FWRV_AMBIGUOUS)
        colnames(dtOutput) <- RENAME_FWRV
      }
    }




    # ADD INSERTION SITE OVERHANG SEQUENCE TO DATA ----------------------------


    # Find insertion sites without overhang sequence value
    dtUnambiguousPositions <- dtOutput[is.na(dtOutput$TIS_seq),]
    # dtUnambiguousPositions = (Subset of) insertion sites with NA values for
    #                          the overhang sequence


    # Correct data
    if (nrow(dtUnambiguousPositions) != 0) {
      dtUnambiguousPositions <- updateInsertionSiteData(dt = dtUnambiguousPositions,
                                                        gr = dtGenomicRanges,
                                                        overhang = overhang,
                                                        readsPerTIS = readsPerTIS)


      # Re-join data
      dtOutput <- dplyr::left_join(dtOutput, dtUnambiguousPositions,
                                   by = c("region_start", "region_end")) %>%
        dplyr::mutate(TIS_seq = dplyr::coalesce(TIS_seq.y, TIS_seq.x))


      # Check for columns indicating a double Tagmentation reaction (Fw/Rv samples)
      if ( !all(COLS_FWRV %in% colnames(dtOutput)) ) {
        dtOutput <- dtOutput %>% dplyr::select(COLS_SINGLE_UMAMBIGUOUS)
        colnames(dtOutput) <- RENAME_SINGLE
      } else if ( all(COLS_FWRV %in% colnames(dtOutput)) ) {
        dtOutput <- dtOutput %>% dplyr::select(COLS_FWRV_UNAMBIGUOUS)
        colnames(dtOutput) <- RENAME_FWRV
      }
    }




    # ADD TAGMENTATION REACTION STRAND ----------------------------------------


    # Create vector to store strands per insertion site
    strand <- c(1:nrow(dtOutput))


    # Iterate over insertion sites, determine dominant strands for forward and/or
    # reverse reaction reads, and return the strand per insertion site as '+', '-'
    # or '*' in ambiguous circumstances
    for (indexTIS in 1:nrow(dtOutput)) {
      # Forward Tagmentation reaction reads
      if (hasForwardReaction) {
        fwStrands <- dplyr::pull(dt[read_name %in% unlist(dtOutput[indexTIS,'read_names']) & reaction == 'Fw', 'strand'])
        fwDominantStrand <- names(which.max(table(fwStrands)))
        fwDominantStrandFraction <- sum(fwStrands == fwDominantStrand, na.rm = TRUE) / length(fwStrands)
      }

      # Reverse Tagmentation reaction reads
      if (hasReverseReaction) {
        rvStrands <- dplyr::pull(dt[read_name %in% unlist(dtOutput[indexTIS,'read_names']) & reaction == 'Rv', 'strand'])
        rvDominantStrand <- names(which.max(table(rvStrands)))
        rvDominantStrandFraction <- sum(rvStrands == rvDominantStrand, na.rm = TRUE) / length(rvStrands)
      }

      # Combined forward and reverse reactions
      if (hasForwardReaction & hasReverseReaction) {
        if (is.null(fwDominantStrand) | is.null(rvDominantStrand)) {
          if (is.null(fwDominantStrand)) {
            strand[indexTIS] <- switchReadStrand(rvDominantStrand)
          } else if (is.null(rvDominantStrand)) {
            strand[indexTIS] <- fwDominantStrand
          }
        } else {
          if ((fwDominantStrand == "+" & rvDominantStrand == "-") | (fwDominantStrand == "-" & rvDominantStrand == "+")) {
            strand[indexTIS] <- fwDominantStrand
          } else {
            strand[indexTIS] <- AMBIGUOUS_STRAND
          }
        }

        # Reverse reaction
      } else if (!hasForwardReaction & hasReverseReaction) {
        strand[indexTIS] <- switchReadStrand(rvDominantStrand)

        # Forward reaction
      } else if (hasForwardReaction & !hasReverseReaction) {
        strand[indexTIS] <- fwDominantStrand
      }
    }
    # Finalise column names.
    dtOutput$strand <- strand
    dtOutput <- dtOutput %>% select(-revmap)
    colnames(dtOutput)[which(colnames(dtOutput) == "TIS_seq")] <- "seq"
  } else {
    dtOutput <- NULL
  }
  return(dtOutput)
}
