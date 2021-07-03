#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          June 09, 2021
# Purpose:       Scrape wikileaks "Iraq War Diaries" (aka SIGACTS III).
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#        - Wikileaks is missing sigacts data for May 2004 and March 2009... gr.
#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# ADMINISTRATIVE                                                          ----
#-----------------------------------------------------------------------------#

#---------------------------#
# Clear working environment
#---------------------------#
rm(list = ls())

#---------------------------#
# Load required packages
#---------------------------#
library(tidyverse)
library(rvest)

#---------------------------#
# Local functions
#---------------------------#

get_event_urls <- function(index_url){
  # This function extracts all unique event links from an index page.
  # These urls can be individually passed to get_event_data()

  index_url <- sprintf("https://wikileaks.org%s",index_url) %>% read_html()

  # Pull each event hyperlink from the index page
  # sub_links <- index %>% html_nodes(".c2") %>% html_nodes("a") %>% html_attr("href")
  sub_links <- index_url %>% html_nodes(".big") %>%  html_nodes("a") %>% html_attr("href")

  # Account for formatting differences per page - filter out non-event links.
  # sub_links <- sub_links[str_detect(sub_links, "event")]
  sub_links <- sub_links[str_detect(sub_links, "event") & !str_detect(sub_links, "type")]
  closeAllConnections()
  return(sub_links)
}

get_event_data <- function(event_url){
  # This function takes the url for each specific SIGACT event and
  # returns a tidied data frame organizing the data.

  # Load the html page
  event_url <- paste0("https://wikileaks.org", event_url) %>% read_html()

  dat_main <- event_url %>% html_nodes(".report") %>% html_table() %>% as.data.frame() %>%
    group_by(Reference.ID) %>%
    pivot_wider(names_from = Var.10,
                values_from = Enemy:Host.nation)

  dat_aux <- event_url %>% html_nodes("code") %>% html_text()
  othr    <- dat_aux[2]
  desc    <- dat_aux[1]

  pat <- c("Report key: ", "Tracking number: ", "Attack on: ", "Complex atack: ",
           "Reporting unit: ", "Unit name: ", "Type of unit: ", "Originator group: ",
           "Updated by group: ", "MGRS: ", "CCIR: ", "Sigact: ", "DColor: ")

  othr <- str_split(othr, pattern = paste(pat, collapse = "|"))[[1]][2:14]
  names(othr) <- str_remove_all(pat, ": ")

  othr <- bind_rows(othr)

  final <- bind_cols(dat_main, othr) %>% mutate(description = desc) %>%
    mutate(across(.cols = everything(), .fns = ~(na_if(., ""))))
  closeAllConnections()
  return(final)
}
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# ORGANIZATION                                                            ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# STEP 1 - Collect "Page 0" for all events organized "by date"
# ----------------------------------- #
index_zeros <- read_html("https://wikileaks.org/irq/#by_date") %>%
  html_nodes("#by_date") %>% html_nodes("a") %>% html_attr("href")
# ----------------------------------- #

# ----------------------------------- #
# STEP 2 - For each starting page, collect the number of event index pages
#          per year-month for later event extraction.
# ----------------------------------- #
index_pages <- c()

for(i in 1:length(index_zeros)){
  # Progress message
  cat(sprintf("\14Working on page %s of %s", i, length(index_zeros)))

  # Create full link
  index_tmp <- paste0("https://wikileaks.org", index_zeros[i]) %>% read_html()

  # Extract number of index pages per year-month date
  page_tmp <-
    # Page navigator appears twice (top and bottom), only need it onces
    (index_tmp %>% html_nodes(".paginator") %>% html_text())[1] %>%

    # It returns as a character vector, extract numbers
    str_extract_all(., "(\\d)+") %>%

    # The last number is the number of pages
    unlist() %>% tail(., 1) %>% as.numeric()

  # The url count starts at "0", so need to subtract 1
  page_tmp <- page_tmp - 1

  # Save to the storage vector
  index_pages <- c(index_pages, page_tmp)

  # Pause to give the server a break
  Sys.sleep(runif(1, 1,2))

};rm(i, index_tmp, page_tmp)
# ----------------------------------- #


# ----------------------------------- #
# STEP 3 - Organize data for later iteration
# ----------------------------------- #
seeds <- data.frame(index = index_zeros,
                    pages = index_pages,
                    month = str_extract(index_zeros, pattern = "_(\\d+)_") %>% str_remove_all(., "_"),
                    year  = str_extract(index_zeros, pattern = "(\\d+)"))

# To loop later I need a modified version of the URL where I can easly
# pass in a page number:
seeds$index_template <- str_replace(string      = seeds$index,
                                    pattern     = "_0[.]",
                                    replacement = "_%s.")

# Clean up and save
rm(index_pages, index_zeros)

save(seeds, file = "Scripts/sigacts-wikileaks/data/seeds.rdata")
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# COLLECT EVENT URL LINKS                                                 ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# "seeds" contains contains the urls for each month's "page 0" as well as the
# number of pages of events that occurred in that particular month.

# Here I will loop through each page containing events for every month and
# collect the urls to individual event pages. These urls can then be passed to
# get_event_data() for the final data collection pass.
# ----------------------------------- #

# Create a storage data frame for event urls
event_urls <- data.frame()

# Since there are so many pages [7799], I breakdown this scrape per year to
# perform in parts
seeds_subset <- seeds %>% filter(year == 2009)

for(i in 1:nrow(seeds_subset)){

  tmp_month <- seeds_subset$month[i]
  tmp_year  <- seeds_subset$year[i]
  tmp_pages <- seeds_subset$pages[i]

  for(z in 0:seeds_subset$pages[i]){
    # Progress message
    cat(sprintf("\14Date: %s-%s\nPage: %s of %s",
                tmp_month, tmp_year,
                z, tmp_pages))

    # Construct per-page index url
    index_tmp <- sprintf(seeds_subset$index_template[i], z)

    # Get event urls from each page
    event_urls_tmp <- tryCatch(
      expr = {get_event_urls(index = index_tmp)},
      error = function(e){NA},
      warning = function(w){NA}
    )

    # Tidy extract:
    event_urls_tmp <- data.frame(event_url = event_urls_tmp,
                                 month     = tmp_month,
                                 year      = tmp_year,
                                 url_page  = z)

    # Save urls
    event_urls <- bind_rows(event_urls, event_urls_tmp)

    # Pause to give the server a break before moving to next page
    Sys.sleep(runif(1, 0.5, 1.25))
  }
  # Pause to give the server a break before moving to next month-year date
  Sys.sleep(runif(1, 0.5, 1.25))
};rm(i, z, index_tmp, event_urls_tmp, tmp_month, tmp_year, tmp_pages)


# Save - overwrite for piecemeal extract
save(event_urls, file = "Scripts/sigacts-wikileaks/data/event_urls.rdata")

rm(seeds_subset)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# EVENT EXTRACTION                                                       ----
#-----------------------------------------------------------------------------#

# Year-wise extraction of events:
{event_urls_subset <- event_urls %>% filter(year == 2005, month == "01")
target <- nrow(event_urls_subset)}
# Storage vector for event data
# sigacts <- data.frame()

for(i in 1:target){
  # Progress message
  cat(sprintf("\14Percent complete: %s", round((i/target)*100,3)))

  # Get event data
  event <- tryCatch(
    expr    = {get_event_data(event_url = event_urls_subset$event_url[i])},
    error   = function(e){NA},
    warning = function(w){NA}
  )

  if(is.na(event)){
    event <- as.list(rep(NA, ncol(sigacts))) %>%
      as.data.frame()
    colnames(event) <- colnames(sigacts)

    event$event_url <- event_urls_subset$event_url[i]
    # bad_urls <- c(bad_urls, event_urls_subset$event_url[i])
  }

  # Save event data to main data frame
  sigacts <- bind_rows(sigacts, event)

  # Give the server a break
  Sys.sleep(runif(1, 0.5, 1.25))
};rm(i, event)

table(lubridate::month(sigacts$Date), lubridate::year(sigacts$Date))
table(event_urls$month, event_urls$year)

save(sigacts, file = "Scripts/sigacts-wikileaks/data/sigacts.rdata")
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
#save.image()
#rm(list = ls())
#-----------------------------------------------------------------------------#

rm(list=ls())
load("Scripts/sigacts-wikileaks/data/event_urls.rdata")
load("Scripts/sigacts-wikileaks/data/sigacts.rdata")
