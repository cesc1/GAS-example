##### Funcions per fer taula descriptiva univariant #####
#### Nomes ajuda una mica a crear una taula descriptiva univariant compacta, però no és essencial per
# l'estudi genètic. La funció principal és create_table(). Idea de posar les funcions en una clase.
# S'ha reaprofitat el codi que vaig crear en l'assignatura de dades transversals.


### Name: cat_table
## Desc: Crea una taula descriptiva per xtable de frequencies absolutes, relatives... 
   # d'una variable categorica
## Packages: NA
## Arguments: 
  # x: vector character/factor
  # var_name: El nom que donar a la variable
  # biskip: Si la variable és dicotomica, tenim la opció de printar nomes una categoria.
## Retorna: data.frame

cat_table <- function(x, var_name, biskip = FALSE) {
  data <- data.frame(table(x))
  colnames(data) <- c("Levels", "fA")
  data$fR <- data$fA / sum(data$fA) * 100
  data$fR <- sprintf("%.1f", data$fR)
  
  if(biskip) {
    if(nrow(data) == 2){
      data <- data[1, ]
    }
  }
  Variable <- c(var_name, rep(NA, nrow(data) - 1))
  data <- cbind(Variable, data)
  return(data)
}


# Funcio de suport per num_table
is_simetric <- function(x) {
  skewness <- sum((x - mean(x))^3) / length(x) / sd(x)^3
  return(abs(skewness) < 1)
}


### Name: num_table
## Desc: Crea una taula descriptiva per xtable per variable numerica
## Packages: NA
## Arguments: 
  # x: vector numeric/integer
  # var_name: El nom que donar a la variable
  # estadistic: l'estadistic que volem calcular. Si null, s'utilitzara mitjana o mediana 
  # depenent de la simetria
## Retorna: data.frame

num_table <- function(x, var_name, estadistic = NULL) {
  x <- x[complete.cases(x)] # de moment s'eliminen NAs, pero es pot millorar
  if(is.null(estadistic)) {
    if(is_simetric(x)) {
      estadistic <- mean
      Levels = "mean"      
    } else {
      estadistic <- median
      Levels = "median"      
    }
  }
  
  m <- estadistic(x)
  q <- quantile(x, probs = c(0.25, 0.75))
  data <- data.frame(Variable = var_name, Levels = Levels, 
                     fA = sum(complete.cases(x)), 
                     fR = sprintf("(%.1f, %.1f, %.1f)", q[1], m, q[2]))
  return(data)
}


### Name: is_categorical
## Desc: Diu si una variable es categorica/factor/logical.
## Packages: NA
## Arguments: 
  # x: vector
## Retorna: boolean

is_categorical <- function(x) {
  if(is.integer(x) | is.numeric(x)) {
    return(FALSE)
    
  } else if(is.character(x) | is.factor(x) | is.logical(x)) {
    return(TRUE)
    
  } else{stop(paste("No s'ha implementat l'us de la classe:", class(x)))}
}


### Name: create_table
## Desc: Crea descriptius univariants de les variables d'un dataframe
## Packages: purrr
## Arguments: 
  # data: data.frame
  # var_names: Nom o etiqueta de les variables (nom que sortira al printar)
## Retorna: data.frame. Pensat per printar amb xtable.

create_table <- function(data, var_names, fA = FALSE) {
  list_tables <- purrr::map2(data, var_names, function(x, var_name) {
    if(is_categorical(x)){
      cat_table(x, var_name)
    } else {
      num_table(x, var_name)
    }
  })
  full_table <- do.call("rbind", list_tables)
  
  rownames(full_table) <- NULL
  colnames(full_table) <- c("Variable", "Levels/Method", "n", "Percent/quartils")
  if(!fA){
    full_table$fA <- NULL
  }
  return(full_table)
}




#### Bivariate analisis

### Name: bi_cat_table
## Desc: Crea una taula de contingencia.
## Packages: NA
## Arguments: 
  # x: vector character/factor
  # var_name: El nom que donar a la variable
  # biskip: Si la variable és dicotomica, tenim la opció de printar nomes una categoria.
## Retorna: data.frame

bi_cat_table <- function(x, resposta, var_name) {
  data <- t(as.matrix(table(resposta, x)))
  
  Labels <- c(rownames(data))
  Variable <- c(var_name, rep(NA, nrow(data) - 1))
  data <- cbind(Variable, Labels, data)
  rownames(data) <- NULL
  return(data)
}
