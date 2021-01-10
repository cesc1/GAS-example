##### Conjunt de funcions per facilitar GEA #####
## Principalment s'utilitza el paquet SNPassoc



#### Data Cleaning


### Name: remove_NA
## Desc: Elimina els SNPs que superen cert percentatge d'NA
## Packages: SNPassoc
## Arguments: 
  # data: objecte setupSNP
  # maxPerc: Numeric percentatge maxim d'NAs. 0_100 numeric
## Retorna: Objecte setupSNP

remove_NA <- function(data, maxPerc = 1){
  # Class check
  if(!is(data, "setupSNP")){
    stop("Class error. Object data needs to be from class setupSNP (library SNPassoc).")
  }
  # Evita l'autoprint de summary
  invisible(capture.output(
    res <- summary(data)))
  
  # Agafa els IDs que no volem, i es passa a bolean
  snpNamesOut <- rownames(res)[res$missing > maxPerc]
  columnOut <- colnames(data) %in% snpNamesOut
  
  return(data[, !columnOut])
}


### Name: remove_HWE
## Desc: Elimina els SNPs que no compleixen les condicions de HWE sobre la distribucio dels controls
## Packages: SNPassoc
## Arguments: 
  # data: objecte setupSNP
  # strat_name: Nom de la columna si volem strats. character
  # reflevel: Reference level name. character (pot ser int, ej: 2 per nivell ref)
  # signif: Nivell de significació del test. numeric. 
## Retorna: Objecte setupSNP


remove_HWE <- function(data, strat_name, reflevel = 2, signif = 0.05){
  # Class check
  if(!is(data, "setupSNP")){
    stop("Class error. Object data needs to be from class setupSNP (library SNPassoc).")
  }
  # Es decideix si volem stratificar.
  if(missing(strat_name)){
    res <- tableHWE(data)
    reflevel <- "HWE (p value)"
  } else{
    res <- tableHWE(data, data[, strat_name])
  }
  
  # Agafa els IDs que no volem, i es passa a bolean
  snpNamesOut <- rownames(res)[res[, reflevel] < signif]
  columnOut <- colnames(data) %in% snpNamesOut
  
  return(data[, !columnOut])
}


### Name: associationOut
## Desc: Funció que crea els models logistics per un SNP i extrau la informació del model mes significatiu.
## Packages: SNPassoc
## Arguments: 
  # base_formula: Formula de la part ajustada. character. ej: "caco ~ gender + sex". character
  # snpID: ID de l'snp
  # data: Objecte setupSNP.
  # ajust: si no volem ajustar, podem posar ajust = FALSE. La formula sera ej: "caco ~". boolean
## Retorna: data.frame

associationOut <- function(base_formula, snpID, data, ajust = TRUE){
  formula <- as.formula(paste(base_formula, 
                              ifelse(ajust, "+", ""), "get(snpID)"))
  mod_dom <- association(formula, data = data, model = "dominant")
  mod_rec <- association(formula, data = data, model = "recessive")
  mod_add <- association(formula, data = data, model = "log-additive")
  
  # Extraiem els valors del model
  matrx_res <- rbind(c(mod_dom[3, "OR"], mod_dom[2, "p-value"]), 
                     c(mod_rec[3, "OR"], mod_rec[2, "p-value"]), 
                     c(mod_add[2, c("OR", "p-value")]))
  # S'ajunta en un dataframe per treballar mes comodament
  data_res <- data.frame(snpID = snpID,
                         Model_Type = c("Dominant", "Recessiu", "Additiu"), 
                         matrx_res)
  colnames(data_res) <- c("ID", "Model_type", "OR", "P_value")
  
  # Seleccio del millor model
  min_pval <- which.min(data_res$P_value)
  return(data_res[min_pval, ])
}


### Name: associationAll
## Desc: Funció que crea els models logistics per tots els SNP, aplicant associationOut()
## Packages: SNPassoc
## Arguments: 
  # base_formula: Formula de la part ajustada. character. ej: "caco ~ gender + sex". character
  # snpIDs: IDs dels SNPs. character vector.
  # data: Objecte setupSNP.
  # ajust: si no volem ajustar, podem posar ajust = FALSE. La formula sera ej: "caco ~". boolean
## Retorna: data.frame

associationAll <- function(base_formula, snpIDs, data, ajust = TRUE){
  taula_res <- snpIDs %>%
    lapply(associationOut, 
           base_formula = base_formula, data = data, ajust = ajust) %>%
    do.call(what = rbind)
  taula_res$P_value_adj <- p.adjust(taula_res$P_value, method = "BH")
  return(taula_res)
}




### Name: rounding_Xtable
## Desc: Arrodoneix els digits de les columnes que volem, i els transforma tipus char (per Xtable).
  # També es pot ordenar una de les columnes a priori. 
## Packages: NA
## Arguments: 
  # data: Objecte data.frame
  # digits: Llista amb el nombre de digits que volem com en C. Ha de ser una llista anomenada
    # amb la columna corresponent: ej: list(OR = ".2", P_value = 0.3)
## Retorna: Taula amb les columnes arrodonides com a characters. data.frame.

rounding_Xtable <- function(data, digits, ordenar = NULL){
  # Datacheck
  if(!is.data.frame(data)){
    if(is.matrix(data)){
      data <- as.data.frame(data)
    }else{
      stop("La classe de data ha de ser matrix/data.frame")
    }
  }
  # Digit names check
  if(any(!names(digits) %in% names(data))){
    stop("Algun dels noms de digits no coincideix amb data")
  }
  # Ordenant, si ho volem
  if(!is.null(ordenar)){
    # Length ordenar check
    if(length(ordenar) != 1)
      stop("L'argument ordenar només pot tenir un valor.")
    #  Name ordenar check
    if(!ordenar %in% names(data))
      stop(sprintf("Variable %s no esta en data.", ordenar))
    # Ordena
    data <- data[order(data[, ordenar]), ]
  }
  # Arrodonint
  columnes_arrod <- sapply(names(digits), function(colNom){
    outString <- paste0("%", digits[[colNom]], "f")
    return(sprintf(outString, data[, colNom]))
  })
  # Substituint columnes
  data[, names(digits)] <- data.frame(columnes_arrod)
  return(data)
}


## PCA

gen_to_num <- function(data){
  numdata <- data %>% # Model aditiu
    sapply(as.numeric) %>%
    subtract(1) %>%
    data.frame() %>%
    extract(complete.cases(data), ) # De moment apliquem casos complets per no complicar
  return(numdata)
}

