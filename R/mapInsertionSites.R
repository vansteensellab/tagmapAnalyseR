#' Map transposon insertion sites.
#'
#' Maps transposon insertion sites (TIS) in high detail from locations called by
#' the pipeline.
#'
#' Input data table \code{dt} of filtered R2 reads (MAPQ >= 10, proper pairs, no
#' PCR duplicates) where the transposon overhang sequence and position is
#' reported, are processed in this function. TISes can be filtered based on read
#' depth, ambiguous insertion positions or overhang sequences are corrected
#' for, and more detailed metadata per TIS are reported.
#'
#' @author Koen Rademaker, \email{k.rademaker@nki.nl}
#' @param dt A data.table with TIS data taken directly from
#' \code{\link{readCalledInsertions}}.
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
#' \item{chr, start, end}{Exact genomic position of the insertion site.}
#' \item{TIS_seq}{Overhang sequence of the TIS.}
#' \item{strand}{Strand on which the insertion has taken place.}
#' \item{region_start, region_end}{Exact genomic position of reads around TIS.}
#' \item{read_names}{Names of individual reads supporting a particular TIS.}
#' \item{revmap}{TO-DO}
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


  # Column names expected for particular single or double Tagmentation reactions
  COL_MARKER_TWO_TAGMAP <- c('read_count_1.x','mapq_1.x','read_count_2.x','mapq_2.x')
  COL_SELECT_TWO_TAGMAP_SHIFTED <- c('seqnames','start','end','width',
    'TIS_seq.y','strand.x','region_start','region_end','read_names.x','revmap.x',
    'read_count.x','mapq.x','read_count_1.x','mapq_1.x','read_count_2.x',
    'mapq_2.x')
  COL_SELECT_TWO_TAGMAP <- c('chr.x','start.x','end.x','width.x',
    'TIS_seq','strand.x','region_start','region_end','read_names.x','revmap.x',
    'read_count.x','mapq.x','read_count_1.x','mapq_1.x','read_count_2.x',
    'mapq_2.x')
  COL_RENAME_TWO_TAGMAP <- c('chr','start','end','width','TIS_seq','strand',
    'region_start', 'region_end', 'read_names','revmap','read_count','mapq',
    'read_count_1', 'mapq_1','read_count_2', 'mapq_2')
  COL_SELECT_ONE_TAGMAP_SHIFTED <- c('seqnames','start','end','width','TIS_seq.y',
    'strand.x','region_start','region_end','read_names.x','revmap.x',
    'read_count.x','mapq.x')
  COL_SELECT_ONE_TAGMAP <- c('chr.x','start.x','end.x','width.x','TIS_seq',
    'strand.x','region_start','region_end','read_names.x','revmap.x','read_count.x','mapq.x')
  COL_RENAME_ONE_TAGMAP <- c('chr','start','end','width','TIS_seq','strand',
    'region_start','region_end','read_names','revmap','read_count',
    'mapq')



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

  # Get start/end of reads in region around insertion site
  dtOutput <- dtOutput %>%
    dplyr::rowwise() %>%
    dplyr::mutate(region_start = min(readsR1[readsR1$QNAME %in%
                                             read_names,c('POS','PNEXT')])) %>%
    dplyr::mutate(region_end = max(readsR1[readsR1$QNAME %in%
                                             read_names,c('POS','PNEXT','TLEN')] %>%
                                     mutate(end = min(POS,PNEXT) + abs(TLEN) - 1) %>%
                                     select(end),
                                   dtGenomicRanges[unlist(revmap), ]$read_end))

  # Add column to store insertion site overhang sequence
  dtOutput['TIS_seq'] <- NA




  # PROCESS FORWARD/REVERSE TAGMENTATION REACTION SEPARATELY ----------------



  # Get specific read count, median MAPQ and concordance
  if (length(unique(dt$sample) == 2)) {
    readsOfFirstSample <- dt %>% dplyr::filter(sample == unique(dt$sample)[1]) %>% dplyr::pull('read_name')
    readsOfSecondSample <- dt %>% dplyr::filter(sample == unique(dt$sample)[2]) %>% dplyr::pull('read_name')

    dtOutput <- dtOutput %>% dplyr::rowwise() %>%
      dplyr::mutate(read_count_1 = sum(read_names %in% readsOfFirstSample),
                    read_count_2 = sum(read_names %in% readsOfSecondSample))

    dtOutput <- dtOutput %>% dplyr::rowwise() %>%
      dplyr::mutate(mapq_1 = median(dtGenomicRanges[read_names %in% readsOfFirstSample,]$mapq),
                    mapq_2 = median(dtGenomicRanges[read_names %in% readsOfSecondSample,]$mapq))

    # Concordance (1/2)
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
    if ( !all(COL_MARKER_TWO_TAGMAP %in% colnames(dtOutput)) ) {
      dtOutput <- dtOutput %>% dplyr::select(COL_SELECT_ONE_TAGMAP_SHIFTED)
      colnames(dtOutput) <- COL_RENAME_ONE_TAGMAP
    } else if ( all(COL_MARKER_TWO_TAGMAP %in% colnames(dtOutput)) ) {
      dtOutput <- dtOutput %>% dplyr::select(COL_SELECT_TWO_TAGMAP_SHIFTED)
      colnames(dtOutput) <- COL_RENAME_TWO_TAGMAP
    }
  }




  # ADD INSERTION SITE OVERHANG SEQUENCE TO DATA ----------------------------



  # Find insertion sites without overhang sequence value
  dtNonShiftedPositions <- dtOutput[is.na(dtOutput$TIS_seq),]
  # dtNonShiftedPositions = (Subset of) insertion sites with NA values for
  #                         the overhang sequence

  # Correct data
  dtNonShiftedPositions <- updateInsertionSiteData(dt = dtNonShiftedPositions,
                                                    gr = dtGenomicRanges,
                                                    overhang = overhang,
                                                    readsPerTIS = readsPerTIS)

  # Re-join data
  dtOutput <- dplyr::left_join(dtOutput, dtNonShiftedPositions,
                               by = c("region_start", "region_end")) %>%
              dplyr::mutate(TIS_seq = dplyr::coalesce(TIS_seq.y, TIS_seq.x))

  # Check for columns indicating a double Tagmentation reaction (Fw/Rv samples)
  if ( !all(COL_MARKER_TWO_TAGMAP %in% colnames(dtOutput)) ) {
    dtOutput <- dtOutput %>% dplyr::select(COL_SELECT_ONE_TAGMAP)
    colnames(dtOutput) <- COL_RENAME_ONE_TAGMAP
  } else if ( all(COL_MARKER_TWO_TAGMAP %in% colnames(dtOutput)) ) {
    dtOutput <- dtOutput %>% dplyr::select(COL_SELECT_TWO_TAGMAP)
    colnames(dtOutput) <- COL_RENAME_TWO_TAGMAP
  }

  return(dtOutput)
}
