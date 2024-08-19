#' Generate conversion table: MGI
#'
#' @keywords internal

if (interactive()) {

  library(tidyverse)
  library(RExtra)
  library(here)

  # Get data from Jax
  URL = "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
  df <- vroom(URL, row_names = F)

  # Create alternate symbol, dash replaced w/ dot
  df$alt <- str_replace_all(df$`Marker Symbol`, "-", "\\.")

  # Subset relevant rows & columns
  INCL = c("Pseudogene", "Gene")
  df <- df[df$`Marker Type` %in% INCL, c(1,7,12,13)]
  names(df) <- c("MGI", "symbol", "synonym", "alt")

  # Loop over rows, split synonyms into multiple rows
  data<- NULL
  tmp <- NULL
  for (i in 1:nrow(df)) {
    # Progress
    if (i %% 10000 == 0) {
      print(as.character(i) + " / " + nrow(df))
    }

    # If the alt-symbol occurs in the symbol column, add to the symonym column
    if (df$alt[i] %!in% df$symbol) {
      df$synonym[i] <- df$synonym[i] + "|" + df$alt[i]
    }

    # If no synonyms, record row in tmp
    # If there are synonyms, record a row for each synonym
    if (length(grep("\\|", df$synonym[i])) == 0) {
      tmp <- rbind(tmp, df[i, 1:3])
    } else {
      synonyms <- unique(unlist(str_split(df$synonym[i], "\\|")))
      for (synonym in synonyms[synonyms != "NA"]) {
        tmp <- rbind(tmp, unlist(c(df[i, 1:2], "synonym" = synonym)))
      }
    }
    # Record tmp to output dataframe after N entries
    # No direct recording to data for better performance
    if (nrow(tmp) > 25000) {
      data <- rbind(data, tmp)
      tmp <- NULL
    }
    # Last iteration, finalise output dataframe
    if (i == nrow(df)) {
      data <- rbind(data, tmp)
      rm(tmp)
    }
  }

  # Remove symbols without synonyms
  data = data[!is.na(data$synonym),]

  # Remove rows where synonym and symbol are equal
  data = data[data$symbol != data$synonym,]

  # Get non-unique synonyms
  # These are synonyms that occur for multiple different symbols
  non_unique = data[duplicated(data$synonym),]$synonym

  # Create and save list with symbol->synonyms
  mmu <- list(
    unique = data[data$synonym %!in% non_unique,],
    non_unique = data[data$synonym %in% non_unique, ],
    all_symbols = unique(data$symbol)
  )
  save(mmu, file = here("./data/mmu.rda"), compress = T, compression_level = 9)
}
