vers <- 'v008'
dir <- 'analysis/manuscripts/paper1/'
rmarkdown::render(
  input       = paste0(dir, 'elktoe_shells.Rmd'),
  output_file = paste0('drafts/elktoe_shells_', vers, '.pdf' )
  )

