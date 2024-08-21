#' StandardizeGeneNames
#'
#' Standardize gene names for a vector of gene names
#'
#' @param genes A character vector containing gene names
#' @param organism Currently, only mus musculus is supported
#' @return A character vector containing symbols converted from synonyms
#' @importFrom RExtra mapvalues %!in% + duplicated
#' @export
#' @examples
#' StandardizeGeneNames(genelist)

StandardizeGeneNames <- function(
  genes,
  organism = c("mmu"),
  force_unique = FALSE) {

  # Input validation
  if (force_unique & any(duplicated(genes))) {
    stop("Input genes may not contain duplicates when force_unique == TRUE")
  }

  # Select conversion table
  if (organism == "mmu") {
    data = GeneNames::mmu
  }

  # Create df dataframe
  df <- data.frame(input = genes)

  # Keep all input genes that are already symbols
  already_symbol = df$input %in% data$all_symbols

  # First pass only unique conversions
  convtab = data$unique
  df$converted <- df$input
  df$converted[!already_symbol] <- mapvalues(
    df$input[!already_symbol],
    convtab$synonym,
    convtab$symbol)

  # Remove duplicates after first pass conversion
  if (force_unique) {
    duplicates = !already_symbol &
      duplicated(df$converted, allItems = T)
    df$converted[duplicates] = df$input[duplicates]
  }

  # Again keep all input genes that are already symbols
  already_symbol = df$converted %in% data$all_symbols

  # For the second pass
  # First remove symbols from non-unique table if already present
  # Then remove remaining duplicated synonyms (mapvalues cannot handle)
  # Only synonym where duplication was removed due to step 1 remain
  convtab = data$non_unique
  convtab = convtab[convtab$symbol %!in% df$converted,]
  convtab = convtab[!duplicated(convtab$synonym, allItems = T),]
  df$converted[!already_symbol] <- mapvalues(
    df$input[!already_symbol],
    convtab$synonym,
    convtab$symbol)

  # Remove duplicates after second pass conversion
  if (force_unique) {
    duplicates = !already_symbol &
      duplicated(df$converted, allItems = T)
    df$converted[duplicates] = df$input[duplicates]

  }

  pct_cnv <- round(sum(df$input != df$converted)/nrow(df) * 100, 2)
  print("Percentages of genes changed: " + pct_cnv + "%")

  pct_nonsymbol <- round(sum(df$converted %!in% data$all_symbols)/nrow(df) * 100, 2)
  print("Percentage of remaining genes not recognized as symbols: " + pct_nonsymbol + "%")

  return(df$converted)
}

