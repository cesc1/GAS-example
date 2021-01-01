#### Web Scrapping de base de dades dbSNP
## L'objectiu és extreure el gen, cromosoma i posició d'un SNP de la base de dades dbSNP.
## S'ha treballat amb funcions, però pot ser millor crear una classe.



### Name: load_ncbi_datanode
## Desc: Llegeix la web en format html, i selecciona la taula d'informació que volem.
## Packages: rvest, magrittr
## Arguments: 
  # snp: nom de l'snp a trobar. character
## Retorna: Objecte classe xml_nodeset

load_ncbi_datanode <- function(snp){
  # Carrega html
  url <- sprintf("https://www.ncbi.nlm.nih.gov/snp/?term=%s", snp)
  web <- read_html(url)
  # Carrega taula
  datanode <- web %>%
    html_nodes(".snpsum_dl_left_align") %>%
    extract(1)
  # Mira si s'ha trobat la taula (snp)
  if(length(datanode) == 0){
    warning(sprintf("No s'ha trobat l'snp: %s", snp))
    return(NA)
  }
  return(datanode)
}

### Name: extractRow
## Desc: Mira quina de les files coincideix amb el pattern. Serveix de suport per les funcions
  # extractGen() i extractCh().
## Packages: rvest, magrittr. stringr
## Arguments: 
  # ncbi_datanode: objecte classe xml_nodeset
  # row_name: pattern to match. character
## Retorna: Fila que coincideix de la taula html. Integer

extractRow <- function(ncbi_datanode, row_name){
  if(is.na(ncbi_datanode))
    return(NA)
  fila <- ncbi_datanode %>%
    html_nodes("dt") %>%
    html_text() %>%
    str_remove("[: ]") %>% 
    str_extract(row_name) %>%
    is.na() %>%
    not()
  if(sum(fila) != 1){
    if(sum(fila) > 1) 
      warning(sprintf("2 files de %s, mirar-ho, molt extrany.", row_name))
    return(NA)
  }
  return(which(fila))
}


### Name: extractGen
## Desc: De la taula html carregada d'un snp, extreu el gen
## Packages: rvest, magrittr. stringr
## Arguments: 
  # ncbi_datanode: objecte classe xml_nodeset
## Retorna: character

extractGen <- function(ncbi_datanode) {
  if(is.na(ncbi_datanode)) 
    return(NA)
  # Get the gene row
  geneRow <- extractRow(ncbi_datanode, "Gene")
  if(is.na(geneRow)){
    return(NA)
  }
  # Extract gene
  gen <- ncbi_datanode %>%
    html_nodes("dd") %>%
    html_text %>%
    extract(geneRow) %>%
    str_split(" ") %>%
    extract2(1) %>%
    extract(1)
  return(gen)
}


### Name: extractCh
## Desc: De la taula html carregada, extreu el cromosoma i posició
## Packages: rvest, magrittr. stringr
## Arguments: 
  # ncbi_datanode: objecte classe xml_nodeset
## Retorna: chromosome (char) and position (char) in Mb. list

extractCh <- function(ncbi_datanode){
  if(is.na(ncbi_datanode)) 
    return(list(ch = NA, pos = NA))
  # Get the chromosome row
  chRow <- extractRow(ncbi_datanode, "Chromosome")
  if(is.na(chRow)){
    return(NA)
  }
  # Extract chromosome, position
  chromosome <- ncbi_datanode %>%
    html_nodes("dd") %>%
    html_text() %>%
    extract(chRow) %>%
    str_split("\\n") %>%
    extract2(1) %>%
    extract(1) %>%
    str_split(":") %>%
    extract2(1)
  if(length(chromosome) != 2)
    warning("Hi ha algo que no funciona al extraure cromosoma i posició (veure str_split(char, ":")).")
  
  return(list(ch = chromosome[1], pos = as.numeric(chromosome[2])/1e6))
}
