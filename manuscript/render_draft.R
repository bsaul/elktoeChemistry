vers <- 'v016'
dir <- 'manuscript/'
rmarkdown::render(
  input       = paste0(dir, 'elktoe_shells.Rmd'),
  output_file = I(paste0('drafts/elktoe_shells_', vers)),
  output_format = "all"
)
