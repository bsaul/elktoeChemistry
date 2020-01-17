# Render elktoe shells manuscript

vers <- 'v017'
dir <- 'manuscript/'
rmarkdown::render(
  input       = paste0(dir, 'elktoe_shells.Rmd'),
  output_file = I(paste0('drafts/elktoe_shells_', vers)),
  output_format = "all"
)

file.copy(from = "manuscript/elktoe_shells.docx", 
          to = sprintf("manuscript/drafts/elktoe_shells_%s.docx", vers))
unlink("manuscript/elktoe_shells.docx")

file.copy(from = sprintf("manuscript/elktoe_shells_%s.log", vers), 
          to = sprintf("manuscript/drafts/elktoe_shells_%s.log", vers))
unlink(sprintf("manuscript/elktoe_shells_%s.log", vers))