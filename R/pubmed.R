#' PubMed Literature Scraper for TCMDATA
#'
#' @description 
#' Fetches and summarizes PubMed literature based on specific Traditional Chinese Medicine (TCM) 
#' and disease keywords within a defined time frame.
#'
#' @param tcm_name Character. Keywords for the TCM (e.g., "Ginseng").
#' @param disease_name Character. Keywords for the disease (e.g., "Fatigue").
#' @param year_range Integer vector of length 2. The publication year range (e.g., c(2015, 2026)).
#' @param email Character. User email address required by NCBI E-utils.
#' @param retmax Integer. Maximum number of PMIDs to retrieve. Default is 1000.
#' @param ... Additional arguments passed to other methods.
#' @importFrom dplyr bind_rows group_by summarise arrange n desc %>% filter
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @return A list of class 'tcm_pubmed' containing:
#' \item{data}{A data frame with PMID, Title, Year, Journal, and Abstract.}
#' \item{stats}{A summary data frame of publication counts by year.}
#' \item{query}{The constructed search query string.}
#' @export
#'
#' @examples
#' \dontrun{
#' result <- get_pubmed_data(tcm_name = "Artemisinin", 
#' disease_name = "Malaria", 
#' year_range = c(2020, 2025),
#' retmax = 100,
#' email = "your email")
#' 
#' print(str(result))
#' }
#' 
get_pubmed_data <- function(tcm_name, 
                            disease_name, 
                            year_range = NULL, 
                            email = NULL, 
                            retmax = 1000,
                            ...) {
  
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for get_pubmed_data(). Please install it.")
  }

  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required for get_pubmed_data(). Please install it.")
  }
  
  if (is.null(email)) {
    stop("NCBI requires an email address to use E-utils. Please provide one.")
  }
  
  if (is.null(year_range)) {
    current_year <- as.numeric(format(Sys.Date(), "%Y"))
    year_range <- c(current_year - 10, current_year)
    message("Time range not specified. Defaulting to the last 10 years: ", 
            year_range[1], " - ", year_range[2])
  }
  
  # Build Search Query
  query <- paste0("(", tcm_name, "[Title/Abstract]) AND (", disease_name, "[Title/Abstract]) ",
                  "AND ('", year_range[1], "/01/01'[Date - Publication] : '", 
                  year_range[2], "/12/31'[Date - Publication])")
  
  message("Search Query: ", query)
  
  # Get Total Count Only (check)
  check_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  check_res <- httr::GET(check_url, query = list(db = "pubmed", term = query, 
                                                 retmax = 0, retmode = "xml", email = email, ...))
  check_xml <- xml2::read_xml(httr::content(check_res, "text", encoding = "UTF-8"))
  total_count <- as.numeric(xml2::xml_text(xml2::xml_find_first(check_xml, "//Count")))
  
  if (length(total_count) == 0 || total_count == 0) {
    message("No literature found.")
    return(NULL)
  }
  
  if (identical(retmax, "all")) {
    final_retmax <- total_count
    message("User requested ALL records. Will retrieve ", total_count, " articles.")
    if (total_count > 5000) warning("Wait time may exceed 10 minutes...")
  } else {
    final_retmax <- as.numeric(retmax)
  }
  
  message("--------------------------------------------------")
  message("PubMed Total Matches: ", total_count, " articles found.")
  
  # Bias Warning
  if (total_count > final_retmax) {
    warning(paste0("\n[!!! BIAS WARNING !!!]\n",
                   "Total records found (", total_count, ") exceeds your limit (", final_retmax, ").\n",
                   "PubMed returns recent articles first. Older years may be truncated.\n",
                   "The 'statistical' plot will be BIASED.\n",
                   "Solution: Set retmax = ", total_count, " or retmax = 'all' to fix this.\n",
                   "Continuing in 3 seconds... (Press Esc to Stop)"), 
            immediate. = TRUE, call. = FALSE)
    Sys.sleep(3)
  } else {
    message("Coverage: 100% (", total_count, " records found). Trends will be accurate.")
    # If total found is less than limit, only fetch what exists
    final_retmax <- total_count
  }
  
  # Retrieve PMID List 
  search_res <- httr::GET(check_url, query = list(db = "pubmed", term = query, 
                                                  retmax = final_retmax, usehistory = "y", email = email, ...))
  
  search_xml <- xml2::read_xml(httr::content(search_res, "text", encoding = "UTF-8"))
  pmids <- xml2::xml_text(xml2::xml_find_all(search_xml, "//IdList/Id"))
  total_found <- length(pmids)
  
  if (total_found == 0) {
    message("No literature found for the given criteria.")
    return(NULL)
  }
  
  message(total_found, " records found. Starting batch download...")
  
  # Fetch Details
  batch_size <- 200
  all_records <- list()
  
  pb <- utils::txtProgressBar(min = 0, max = total_found, style = 3)
  
  for (start in seq(1, total_found, by = batch_size)) {
    end <- min(start + batch_size - 1, total_found)
    batch_pmids <- paste(pmids[start:end], collapse = ",")
    
    fetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_res <- httr::GET(fetch_url, query = list(db = "pubmed", id = batch_pmids, 
                                                   retmode = "xml", email = email, ...))
    
    Sys.sleep(0.4) 
    
    fetch_xml <- xml2::read_xml(httr::content(fetch_res, "text", encoding = "UTF-8"))
    articles <- xml2::xml_find_all(fetch_xml, "//PubmedArticle")
    
    batch_df <- lapply(articles, function(article) {
      year_node <- xml2::xml_find_first(article, ".//PubDate/Year")
      if (is.na(xml2::xml_text(year_node))) {
        medline_date <- xml2::xml_text(xml2::xml_find_first(article, ".//PubDate/MedlineDate"))
        pub_year <- regmatches(medline_date, regexpr("\\d{4}", medline_date))
      } else {
        pub_year <- xml2::xml_text(year_node)
      }
      
      data.frame(
        PMID = xml2::xml_text(xml2::xml_find_first(article, ".//PMID")),
        Title = xml2::xml_text(xml2::xml_find_first(article, ".//ArticleTitle")),
        Year = as.numeric(pub_year),
        Journal = xml2::xml_text(xml2::xml_find_first(article, ".//Journal/Title")),
        Abstract = xml2::xml_text(xml2::xml_find_first(article, ".//AbstractText")),
        stringsAsFactors = FALSE
      )
    })
    
    all_records <- c(all_records, batch_df)
    utils::setTxtProgressBar(pb, end)
  }
  close(pb) 
  
  # Summary and Statistics
  final_df <- bind_rows(all_records) %>% 
    filter(!is.na(.data$Year))
  
  summary_stats <- final_df %>%
    group_by(.data$Year) %>%
    summarise(Count = dplyr::n(), .groups = 'drop') %>%
    arrange(dplyr::desc(.data$Year))
  
  # Construct S3 Object
  result <- list(data = final_df, stats = summary_stats, query = query)
  class(result) <- "tcm_pubmed"
  
  return(result)
}


#' Plot PubMed Search Results
#'
#' @description 
#' Visualizes the results of a PubMed search, showing the annual publication trend or the top journals.
#'
#' @param x An object of class 'tcm_pubmed'.
#' @param type Character. Type of plot to return: "trend", or "journal". Default is "trend".
#' @param N Integer. The number of top journals to display in the bar chart. Defaults to 10.
#' @param bar_col Character. Hex code or name for the fill color of the bars in the journal chart. Defaults to "#3498db".
#' @param line_col Character. Hex code or name for the color of the trend line. Defaults to "#2c3e50".
#' @param point_col Character. Hex code or name for the color of the points on the trend line. Defaults to "#e74c3c".
#' @param ... Additional arguments passed to other methods.
#' 
#' @import ggplot2
#' @importFrom dplyr group_by summarise n slice_max %>%
#' @importFrom stats reorder
#' @importFrom rlang .data
#' 
#' @return A ggplot object or a list of plots.
#' @export
#' @examples
#' \dontrun{
#' ## get result from pubmed
#' result <- get_pubmed_data(tcm_name = "Artemisinin", 
#' disease_name = "Malaria", 
#' year_range = c(2020, 2025),
#' retmax = 100,
#' email = "your email")
#' 
#' ## plot
#' p1 <- plot(result, type = "trend")
#' p1
#' 
#' p2 <- plot(result, type = "journal")
#' p2
#' 
#' aplot::plot_list(p1, p2)
#' }
#' 
plot.tcm_pubmed <- function(x, 
                            type = c("trend", "journal"), 
                            N = 10,
                            bar_col = "#4E79A7",
                            line_col = "#2F4858",
                            point_col = "#B85C5C",
                            ...) {
  
  type <- match.arg(type)
  
  plot_data <- x[["stats"]]
  raw_data <- x[["data"]]
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    message("The object contains no stats data.")
    return(NULL)
  }
  
  p_final <- NULL
  
  ## trend plot
  if (type == "trend") {
    year_breaks <- if (nrow(plot_data) > 1) {
      seq(min(plot_data$Year), max(plot_data$Year), by = 1)
    } else {
      plot_data$Year
    }
    
    p_final <- ggplot(plot_data, aes(x = .data$Year, y = .data$Count)) +
      geom_line(color = line_col, linewidth = 0.55) +
      geom_point(color = point_col, size = 1.4) +
      scale_x_continuous(breaks = year_breaks) +
      labs(x = NULL, y = "Count") +
      .theme_tcm_pub(grid = "y")
  }
  
  ## journal plot
  if (type == "journal") {
    if (is.null(raw_data)) {
      message("The object contains no raw data for journal plot.")
      return(NULL)
    }
    
    journal_stats <- raw_data %>%
      group_by(.data$Journal) %>%
      summarise(Count = dplyr::n(), .groups = 'drop') %>%
      slice_max(.data$Count, n = N)
    
    p_final <- ggplot(journal_stats, 
                      aes(x = stats::reorder(.data$Journal, .data$Count), y = .data$Count)) +
      geom_col(fill = bar_col) +
      coord_flip() +
      labs(x = NULL, y = "Count", title = paste("Top", N, "Journals")) +
      .theme_tcm_pub(grid = "x")
  }
  return(p_final)
}


#' Extract and Sort PubMed Literature Table
#' @description 
#' Accesses the literature data within a 'tcm_pubmed' object and returns it as 
#' a data frame sorted by publication year in descending order.
#' @param x An object of class 'tcm_pubmed' generated by \code{get_pubmed_data}.
#' @param n Integer. Number of top records to return. If \code{NULL} (default), returns all records.
#' @param file Character. Optional file path (e.g., "results.csv"). If provided, the table will be saved to this location.
#' @importFrom dplyr arrange desc slice_head %>%
#' @importFrom rlang .data
#' @importFrom utils write.csv
#' @return A data frame containing the sorted literature records.
#' @export
#'
#' @examples
#' \dontrun{
#' ## get result from pubmed
#' result <- get_pubmed_data(tcm_name = "Artemisinin", 
#' disease_name = "Malaria", 
#' year_range = c(2020, 2025),
#' retmax = 100,
#' email = "your email")
#' 
#' # Get the sorted table
#' lit_table <- get_pubmed_table(result)
#' 
#' # Get top 10 records and save to CSV
#' lit_table <- get_pubmed_table(result, n = 10, file = "top10_lit.csv")
#' }
#' 
get_pubmed_table <- function(x, n = NULL, file = NULL) {
  
  if (!inherits(x, "tcm_pubmed")) {
    stop("Input must be an object of class 'tcm_pubmed'.")
  }
  
  df <- x[["data"]]
  if (is.null(df) || nrow(df) == 0) {
    message("The object contains no data.")
    return(NULL)
  }
  
  if (is.null(n)){
    n <- nrow(df)
  }
  
  sorted_df <- df %>%
    arrange(dplyr::desc(.data[["Year"]])) %>%
    slice_head(n = n)
  
  if (!is.null(file)) {
    tryCatch({
      utils::write.csv(sorted_df, file = file, row.names = FALSE)
      message("Table successfully exported to: ", file)
    }, error = function(e) {
      warning("Failed to export CSV: ", e$message)
    })
  }
  
  return(sorted_df)
}
