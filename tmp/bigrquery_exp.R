
billing <- "sound-mantra-229923"

con <- connection_open(
  bigrquery::bigquery(),
  project = "publicdata",
  dataset = "samples",
  billing = billing
)


natality <- tbl(con, "natality")

natality %>%
  select(year, month, day, weight_pounds) %>% 
  head(10) %>%
  collect() %>%
  identity()
