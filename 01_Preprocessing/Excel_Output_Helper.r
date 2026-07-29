# ============================================================
# STANDARDIZED EXCEL OUTPUT FORMATTING HELPER
# ============================================================
# This script contains helper functions to ensure consistent
# Excel output formatting across all analysis scripts.
#
# USER REQUIREMENTS:
# - One decimal place for all numeric values
# - Font size 14, Helvetica Neue, black text (user to apply)
# - No actual Excel formatting (just data in rows/columns)
# - Keep % and total n in SEPARATE columns
# - Display percentages as "50.0%" not "50.0" or "0.5"
# ============================================================

library(dplyr)

# Format a number to 1 decimal place
format_number <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.1f", x))
}

# Format a percentage to 1 decimal place with % sign
format_percent <- function(x) {
  ifelse(is.na(x), NA_character_, paste0(sprintf("%.1f", x), "%"))
}

# Format a count (integer, no decimals)
format_count <- function(x) {
  ifelse(is.na(x), NA_character_, as.character(as.integer(x)))
}

# Calculate count and percentage, return as separate columns
# Usage: Used in summary tables where you want "N" and "Percent" columns
calc_n_pct <- function(condition, total = NULL) {
  n <- sum(condition, na.rm = TRUE)
  if (is.null(total)) {
    total <- length(condition[!is.na(condition)])
  }
  pct <- if (total > 0) (n / total) * 100 else NA_real_
  
  list(
    N = n,
    Percent = format_percent(pct)
  )
}

# Create a formatted summary table with separate N and Percent columns
# Input: data frame with N_xxx and Pct_xxx columns (numeric)
# Output: data frame with N_xxx (character) and Pct_xxx (character with %)
format_summary_table <- function(df) {
  df %>%
    mutate(across(starts_with("N_") & where(is.numeric), format_count)) %>%
    mutate(across(starts_with("Pct_") & where(is.numeric), format_percent)) %>%
    mutate(across(ends_with("_Pct") & where(is.numeric), format_percent)) %>%
    mutate(across(matches("Percent$") & where(is.numeric), format_percent))
}

# Apply standard formatting to all numeric columns in a data frame
# Excludes count columns (those starting with "N_" or "Total")
# Percentages (containing "Pct" or "Percent") get % sign
format_excel_output <- function(df) {
  df %>%
    mutate(across(
      where(is.numeric) & (starts_with("N_") | starts_with("Total_") | ends_with("_N")),
      format_count
    )) %>%
    mutate(across(
      where(is.numeric) & (contains("Pct") | contains("Percent") | contains("percent")),
      format_percent
    )) %>%
    mutate(across(
      where(is.numeric) & !contains("Pct") & !contains("Percent") & 
        !starts_with("N_") & !starts_with("Total_") & !ends_with("_N"),
      format_number
    ))
}

# Create a note about formatting for Excel sheets
create_formatting_note <- function() {
  data.frame(
    Note = c(
      "FORMATTING GUIDE:",
      "- All numbers shown to 1 decimal place",
      "- Percentages shown with % sign (e.g., 50.0%)",
      "- Recommended font: Helvetica Neue, 14pt, black",
      "- N columns show counts, separate from percentages"
    )
  )
}

# Example usage in scripts:
# Instead of:
#   mutate(Pct = round(100 * N / Total, 1))
#
# Use:
#   mutate(Pct = (N / Total) * 100) %>%
#   format_excel_output()
#
# Or for more control:
#   mutate(
#     N = format_count(N),
#     Pct = format_percent(Pct_numeric)
#   )


