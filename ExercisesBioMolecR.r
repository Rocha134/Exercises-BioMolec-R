#Francisco Rocha Juárez. A01730560

dna <- c("a","t","c","g")

# 1 Creación de una función para generar una secuencia aleatoria de acidos
nucleicos
sequence <-function(dna, n){
  return(sample(dna, n, replace = T))
}
#Aplicamos la función para crear una muestra aleatoria de DNA
sequence(dna, 30)

# 2 Creación de una función para calcular la longitud de una secuencia
size <-function(dna){
  print("Numero de bases: ")
  return(length(dna))
}

#Probando la función
gen_molde <- sequence(dna,30)
gen_molde
size(gen_molde)

# 3 Creación de una funcion que calcula el porcentaje de cada nucleotido
en la secuencia
base.percentage <-function(dna){
  a <- 0
  t <- 0
  g <- 0
  c <- 0
  for(i in 1:length(dna)){
    if(dna[i] == "a"){
      a <- a + 1
    } else if(dna[i] =="t"){
      t <- t + 1
    } else if(dna[i] == "c"){
      c <- c + 1
    } else if(dna[i] == "g"){
      g <- g + 1
    }
  }
  print("Porcentaje de Adenina: ")
  print((a/length(dna))*100)
  print("Porcentaje de Timina: ")
  print((t/length(dna))*100)
  print("Porcentaje de Citosina: ")
  print((c/length(dna))*100)
  print("Porcentaje de Guanina: ")
  print((g/length(dna))*100)
}

# Ejecuto la función
base.percentage(gen_molde)

# 4 Funcion que recibe una hebra directa y regresa la inversa
invert <-function(dna) {
  return(rev(dna))
}

#Ejecutamos
invert(gen_molde)

# 5 Funcion que recibe una hebra directa y obtiene la complementaria
complement <-function(dna){
  cdna <- c()
  for(i in 1:length(dna)){
    if(dna[i] == "a"){
      cdna[i] <- "t"
    } else if(dna[i] == "t"){
      cdna[i] <- "a"
    } else if(dna[i] == "c"){
      cdna[i] <- "g"
    } else if(dna[i] == "g"){
      cdna[i] <- "c"
    }
  }
  return(cdna)
}

# Ejecutamos
complement(gen_molde)
gen_molde

# 6 Programa una funcion que transcribe DNA a RNA
transcription <-function(dna){
  mrna <- c()
  for(i in 1:length(dna)){
    if(dna[i] == "a"){
      mrna[i] <- "u"
    } else if(dna[i] == "t"){
      mrna[i] <- "a"
    } else if(dna[i] == "c"){
      mrna[i] <- "g"
    } else if(dna[i] == "g"){
      mrna[i] <- "c"
    }
  }
  return(mrna)
}

# Ejecutamos
transcription(gen_molde)
gen_molde

# 7 Funcion que traduce una secuencia de RNA a PROTEINA
traduction <-function(mrna){
  codones <- c()
  i <- 1
  j <- 1
  while(i < length(mrna)){
    codones[j] <- sprintf("%s%s%s", mrna[i], mrna[i + 1], mrna[i + 2])
    j <- j + 1
    i <- i + 3
  }
   for (i in 1:length(codones)){
     if (codones[i] == "aug"){
       print("START")
     } else if (codones[i] == "uaa" || codones[i] == "uag" || codones[i]
== "uga"){
       print("STOP")
     } else if (codones[i] != "aug" || codones[i] != "uaa" || codones !=
"uag" || codones != "uga"){
       if(codones[i] == "uuu" || codones[i] == "uuc"){
         print("Phe")
       } else if(codones[i] == "uua" || codones[i] == "uug" || codones[i]
== "cuu" || codones[i] == "cuc" || codones[i] == "cua" || codones[i] ==
"cug" ){
         print("Leu")
       } else if(codones[i] == "auu" || codones[i] == "auc" || codones[i]
== "aua" || codones[i] == "aug"){
         print("lle")
       } else if(codones[i] == "guu" || codones[i] == "guc" || codones[i]
== "gua" || codones[i] == "gug"){
         print("Val")
       } else if(codones[i] == "ucu" || codones[i] == "ucc" || codones[i]
== "uca" || codones[i] == "ucg" || codones[i] == "agu" || codones[i] ==
"agc" || codones[i] == "aga" || codones[i] == "agg"){
         print("Ser")
       } else if(codones[i] == "ccu" || codones[i] == "ccc" || codones[i]
== "cca" || codones[i] == "ccg"){
         print("Pro")
       } else if(codones[i] == "acu" || codones[i] == "acc" || codones[i]
== "aca" || codones[i] == "acg"){
         print("Thr")
       } else if(codones[i] == "gcu" || codones[i] == "gcc" || codones[i]
== "gca" || codones[i] == "gcg"){
         print("Ala")
       } else if(codones[i] == "uau" || codones[i] == "uac"){
         print("Tyr")
       } else if(codones[i] == "cau" || codones[i] == "cac"){
         print("His")
       } else if(codones[i] == "caa" || codones[i] == "cag"){
         print("Gin")
       } else if(codones[i] == "aau" || codones[i] == "aac"){
         print("Asn")
       } else if(codones[i] == "aaa" || codones[i] == "aag"){
         print("Lys")
       } else if(codones[i] == "gau" || codones[i] == "gac"){
         print("Asp")
       } else if(codones[i] == "gaa" || codones[i] == "gag"){
         print("Glu")
       } else if(codones[i] == "cgu" || codones[i] == "cgc" || codones[i]
== "cga" || codones[i] == "cgg" || codones[i] == "aga" || codones[i] ==
"agg"){
         print("Arg")
       } else if(codones[i] == "ggu" || codones[i] == "ggc" || codones[i]
== "gga" || codones[i] == "ggg"){
         print("Gly")
       }
     }
   }
}
 
# Ejecutamos
prueba <-transcription(gen_molde)
prueba
traduction(prueba)
