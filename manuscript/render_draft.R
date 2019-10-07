vers <- 'v010'
dir <- 'manuscript/'
rmarkdown::render(
  input       = paste0(dir, 'elktoe_shells.Rmd'),
  output_file = paste0('drafts/elktoe_shells_', vers, '.pdf' )
  )

