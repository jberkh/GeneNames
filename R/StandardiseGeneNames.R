#' StandardiseGeneNames
#'
#' Standardise gene names for a vector of gene names
#'
#' @param genes A character vector containing gene names
#' @param organism Currently, only mus musculus is supported
#' @return A character vector containing symbols converted from synonyms
#' @importFrom RExtra mapvalues %!in% + duplicated
#' @export
#' @examples
#' StandardiseGeneNames(genelist)

StandardiseGeneNames <- function(genes, organism = c("mmu"), force_unique = F) {
  mmu = GeneNames::mmu

  # Create df dataframe
  df <- data.frame(input = genes)

  # Keep all input genes that are already symbols
  already_symbol = df$input %in% mmu$all_symbols

  # First pass only unique conversions
  df$pass1 <- df$input
  df$pass1[!already_symbol] <- mapvalues(
    df$input[!already_symbol],
    mmu$unique$Synonym[!is.na(mmu$unique$Synonym)],
    mmu$unique$Symbol[!is.na(mmu$unique$Synonym)])

  # Force unique output (no duplicates)
  if (force_unique) {
    duplicate = duplicated(df$pass1, allItems = T)
    df$pass1[duplicate] <- df$input[duplicate]
  }

  pct_p1 <- round(sum(df$input != df$pass1)/nrow(df) * 100, 2)
  print("Percentages of genes changed: " + pct_p1 + "%")

  pct_nonsymbol <- round(sum(df$pass1 %!in% mmu$all_symbols)/nrow(df) * 100, 2)
  print("Percentage of input genes not recognized as symbols: " + pct_nonsymbol + "%")

  return(df$pass1)
}

