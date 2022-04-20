# tagmapAnalyseR

Stand-alone package integrated into transposon mapping pipeline [tagmap_hopping](https://github.com/krademaker/tagmap_hopping/tree/snakefile_overhaul) for Tagmentation-Based Mapping (Tagmap) data, as well as stand-alone functionality.

## Background

TagmapAnalyseR is a module inside the tagmap_hopping pipeline, and is particularly responsible to map (processed) TagMap reads to genomic locations and figure out where transposons (e.g. PiggyBac or Sleeping Beauty) integrated into the genome. TagmapAnalyseR combines several tricks to derive a consensus mapping location for integrations are read coverage may be ambiguous towards the precise location. More details can be obtained from the author upon request (methods section in thesis).

## Usage

Dependencies for pipeline are described in the respective repository, see its [map_insertions.R](https://github.com/krademaker/tagmap_hopping/blob/snakefile_overhaul/src/scripts/map_insertions.R) script for practical application of tagmapAnalyseR.

**Installation**

```devtools::install_github("https://github.com/krademaker/tagmapAnalyseR")```


**Load processed read coverage of (putative) integrations**

```dt <- tagmapAnalyseR::readPutativeInsertions(input)```

**Obtain integrations with ambiguous read coverage**

```ambiguous <- tagmapAnalyseR::findAmbiguousInsertionSites(dt, padding = (nchar(overhang_sequence)*2)+2)```

**Map integrations**

```
mapped <- tagmapAnalyseR::mapInsertionSites(dt = dt,
                                            bam = bam_path,
                                            overhang = overhang_sequence,
                                            samtoolsPath = samtools_path,
                                            depth = as.integer(minimal_read_depth),
                                            ambiguousInsertions = ambiguous)
```

## Contact / Support 

Reach out through the issues section or reach out directly to the [author](https://github.com/krademaker).
