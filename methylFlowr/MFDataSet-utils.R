regions <- function(obj) obj@regions
patterns <- function(obj) obj@patterns
components <- function(obj) obj@components

nregions <- function(obj, by.component=TRUE) {
  if (by.component)
    return(table(regions(obj)$cid))
  nrow(regions(obj))
}

npatterns <- function(obj, by.component=TRUE) {
  if (by.component)
    return(table(patterns(obj)$cid))
  nrow(patterns(obj))
}

counts <- function(obj, level=c("region","component"), kind=c("raw","normalized"))
  {
    if (level == "region") {
      switch(kind,
             raw=regions(obj)$coverage,
             normalized=regions(obj)$norm_coverage)
    } else {
      switch(kind,
             raw=components(obj)$coverage,
             normalized=components(obj)$coverage)
    }
  }


