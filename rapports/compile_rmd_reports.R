
gh_render <- function(file) {
  
  cat("Caution notice :
  This script has a component that rewrites latex math formulas into
  GitHub-compatible URIs. To avoid any bug, you should not have $ followed or
  preceded by a space ( ) anywhere in markdown text parts of your .rmd file if 
  it's not as a latex formula marker like here : $\\pi \\in [1,4]$.\n\n")
  
  if(!grepl(".Rmd$", file)) cat("File must be a Rmarkdown file (.Rmd extension).") ; break
  
  needed_packages = list("parsermd","rmarkdown")
  missing_package = !all(needed_packages %in% rownames(installed.packages()))
  
  for (p in needed_packages) {
    
    if(!require(p)) cat(sprintf("The '%s' package has to be installed.\n"), p)
    
  }
  
  if(missing_package) break
  
  require(parsermd)
  require(rmarkdown)
  
  file = "C:/Users/guillert/Desktop/Projets/sourdough_evolution/rapports/2022-01-14_Data_phenot_overview_including_outliers.Rmd"
  
  rmd_file <- parsermd::parse_rmd(file)
  
  rmd_file
  
  
  # Encodes inline LaTeX formulas
  
  inline_encoding <- function(x, type) {
    
    if (type %in% c("rmd_markdown","rmd_heading")) {
      
      tmp <- gsub(x = x, "(?<=( )|(\\()|(^))(\\$)", 
                  "<img src='https://render.githubusercontent.com/render/math?math=", 
                  perl = T) %>%
        
        gsub(x = ., "(\\$)(?=( )|(\\))|($)|(\\.)|(\\,))", 
             "'>", 
             perl = T)[1]
      
      "(?<=\\$)[^\\$\n]+(?=\\$) < to match what is between the dollars. to modify"
      
    } else {
      
      tmp <- x
      
    }
    
    return(tmp)
    
  }
  
  as_tibble(rmd_file)$ast[[22]]
  inline_encoding(as_tibble(rmd_file)$ast[[22]]$name, "rmd_heading")
  
  to_edit <- rmd_file %>%
    as_tibble() %>%
    mutate(ast = map2(ast, type, ~ .x$code))
  

  parsermd::render(x = ., name = "test.md")
  
  
  # Compiles local version
  rmarkdown::render(input = file, output_format = "all")
  
  
  # TODO:
  # - verify that the corresponding .md file exists and is up-to-date, or break.
  # - import it with readLines
  # - list regexs
  # - replace $s by html balises
  # - encodeURIs
  # - export with _gh.md extension
  
  if(!grepl(".Rmd$", file)) cat("File must be a Rmarkdown file (.Rmd extension).") ; break
  
  
  tmp_rmd <- readLines(file)
  
  latex_formulas_patterns <- c(inline = "$.$",
                               multiline = "\$\$([^\$]+)\$\$")
  
  tmp_rmd <- gsub(perl = TRUE)
  
  
  # Compiles Github version
  write(tmp_rmd, "tmp.Rmd")
  rmarkdown::render(input = "tmp.Rmd", 
                    output_format = "md_document", 
                    output_file = sub(".Rmd$", "_gh.md$", file))
  file.remove("tmp.Rmd")
  
}

