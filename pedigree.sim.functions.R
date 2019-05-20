#CHANGE THIS TO YOUR PEDIGREESIM INSTALLATION PATH
pedigreesim.path = "/home/USERNAME/PEDIGREESIM/"






#####################################
#####################################
full.path = pedigreesim.path
setwd(full.path)
load(paste0(full.path, "pedsim.files/axc.chr1a.cm"))
load(paste0(full.path, "pedsim.files/cxa.6b.cm"))

#create directory structure
if(!dir.exists(paste0(full.path, "analyses"))) dir.create("analyses")
if(!dir.exists(paste0(full.path, "config.files"))) dir.create("config.files")
if(!dir.exists(paste0(full.path, "simresults"))) dir.create("simresults")
if(!dir.exists(paste0(full.path, "example.config"))) dir.create("example.config")
if(!dir.exists(paste0(full.path, "10k.sims"))) dir.create("10k.sims")

#### MISC FUNCTIONS ####

p = function(...){
  #quicker paste function. must supply vector of 
  #elements to paste together (i.e. with c())
  paste(..., collapse = "", sep = "")
}

convert.to.character.data.frame = function(df){
  df = as.data.frame(df)
  g = df
  g[] = lapply(df, as.character)
  return(g)
}

#### PEDSIM FUNCTIONS ####

run.1000.sims = function(project.name, f.gen, num.sims, w.selec, selec.str, selec.pos, pop.size, cm.pos, perform.analysis.flag){
  #generates multiple genotyping datasets by running pedigreesim with the specified parameters for the specified number of simulations. Performs gamete selection procedure from the specified number of filial generations.

  #args:
  #project.name - character sting, the name of the project  
  #f.gen - integer, the highest filial generation of the population
  #num.sims - integer the number of simulations to perform
  #w.selec - boolean flag, indicates whether selection is to applied
  #selec.str - integer, the strength of selection to be used (ranging from 0 - 1)
  #selec.pos - integer, the marker number against which to apply selection (e.g. 50 for the 50th marker)
  #pop.size - integer, the number of individuals in the population
  #cm.pos - numeric vector, the centimorgan positions of markers
  #perform.analysis.flag - boolean, indicates whether to analyse the resulting genotyping dataframes using the perform.analysis() function
  if(missing(num.sims)) num.sims = 1000
  if(missing(w.selec)) w.selec = F
  if(missing(selec.str)) selec.str = 10
  if(missing(selec.pos)) selec.pos = 200
  if(missing(cm.pos)) cm.pos = axc.chr1a.cm
  if(missing(perform.analysis.flag)) perform.analysis.flag = T
  
  t1 = mclapply(1:num.sims, function(x){
    if(w.selec == T){
      complete.f5.selection.procedure(length(cm.pos), pop.size, f.gen, cm.pos, full.path, paste0(project.name, x), "KOSAMBI", T, selec.str, selec.pos)  
    } else {
      complete.fgen.no.selec.procedure(length(cm.pos), pop.size, f.gen, cm.pos, full.path, paste0(project.name, x), "KOSAMBI", T)
    }
  }, mc.cores = 60)
  
  clean.files()  
  
  # browser()
  
  if(w.selec == T){
    files.to.combine = list.files(p(full.path, "simresults/"), 
                                  pattern = p(project.name, ".*output_genotypes.dat_wselec.f", f.gen, ".dat"))  
  } else {
    files.to.combine = list.files(p(full.path, "simresults/"), 
                                  pattern = p(project.name, ".*output_genotypes.dat"))
  }
  
  
  all.sim.files = mclapply(files.to.combine, function(x){
    
    
    
    g = read.ped(p(full.path, "simresults/", x))
    g = convert.ped.sim.to.zygote(g)
    g = reset.rownames(g)
    rownames(g) = paste0("marker", rownames(g))
    g = as.data.frame(t(g))
    g = g[grep(p("F", f.gen), rownames(g)), ]
  }, mc.cores = 60)
  
  if(class(all.sim.files[[1]]) != "data.frame"){
    all.sim.files = lapply(all.sim.files, as.data.frame)
  }
  
  # browser()
  
  if(f.gen == 2){
    f2flag = T
  } else {
    f2flag = F
  }
  
  # browser()
  
  s(all.sim.files, p(full.path, "simresults/", project.name), "NA")
  
  if(perform.analysis.flag == T){
    if(w.selec == F){
      # perform.analysis = function(geno.list, selection, selection.position, selection.strength, recombination.position, F2, fgen)
      analysis1 = perform.analysis(all.sim.files, recombination.position = "AxC1A", F2 = f2flag, fgen = f.gen)
    } else {
      analysis1 = perform.analysis(all.sim.files, selection = T, selection.position = selec.pos, selection.strength = selec.str, recombination.position = "AxC1A", F2 = f2flag, fgen = f.gen)
    }
    
    write.csv(analysis1, p(full.path, "analyses/", project.name, ".csv"), row.names = F)
  }
  

  
  
  files.to.rm = list.files(p(full.path, "simresults/"), pattern = "genotypes.dat")
  lapply(paste0(p(full.path, "simresults/"), files.to.rm), file.remove)
  
}


make.founder.genotypes = function(num.markers, file.path.to.write){
  #args:
  # num.markers - integer, the number of markers to make genotypes for
  # file.path.to.write - character string, the path to the file to write genotypes
  fg = c("marker\tP1_1\tP1_2\tP2_1\tP2_2", paste0("marker", 1:num.markers, "\tA\tA\tB\tB"))
  write(fg, file.path.to.write, ncolumns = 1)
}

make.pedigree.file = function(f.gen, number.individuals, file.path.to.write){
  #makes a pedigreesim pedigree file for a single seed descent population structure
  #args:
  # f.gen - integer, the filial generation 
  # number.individuals - integer, the number of individuals in the population
  # file.path.to.write - character string, the file path to write pedigree file
  if(f.gen < 2) return("f.gen must be at least 2")
  top.lines = list(c("Name\tParent1\tParent2", "P1\tNA\tNA", "P2\tNA\tNA", "F1_1\tP1\tP2"))
  f2.lines = list(paste0("F2_", 1:number.individuals, "\tF1_1\tF1_1"))
  
  if(f.gen > 2){
    more.lines = lapply(3:f.gen, function(x){
      paste0("F", x, "_", 1:number.individuals, "\tF", (x - 1), "_", 1:number.individuals, "\tF", (x - 1), "_", 1:number.individuals)
    })
    all.lines = c(top.lines, f2.lines, more.lines)
  } else {
    all.lines = c(top.lines, f2.lines)
  }
  
  
  all.lines2 = do.call(c, all.lines)
  write(all.lines2, file.path.to.write, ncolumns = 1)
  
}

make.map.file = function(cm.positions, file.path.to.write){
  #makes a pedigreesim map file
  #args:
  # cm.positions - a numeric vector of postitions of markers in centimorgans
  # file.path.to.write - character string, the file path to write the map file
  lines1 = c("marker\tchromosome\tposition", paste0("marker", 1:length(cm.positions), "\tA\t", cm.positions))
  write(lines1, file.path.to.write, ncolumns = 1)
}

make.chrom.file = function(cm.length, file.path.to.write){
  #make a pedigreesim chromosome file
  #args:
  # cm.length - the length of the chromosome in centimorgans
  # file.path.to.write - character string, the file path to write the chromosome file
  lines1 = c("chromosome\tlength\tcentromere\tprefPairing\tquadrivalents", paste0("A\t", cm.length, "\t", (cm.length / 2), "\t1.00\t0.0"))
  write(lines1, file.path.to.write, ncolumns = 1)
}

make.pedsim.par.file = function(project.name, file.path.to.write, map.function, gen.file.path, ped.file.path){
  #make a pedigree sim parameter file
  #args:
  # project.name - character string, the name of the project
  # file.path.to.write - character string, the file path to write the parameter file
  # map.function - character string, the map function to use (defaults to Kosambi)
  # gen.file.path - character string, the file path of the pedigreesim founder file
  # ped.file.path - character string, the file path of the pedigreesim pedigree file
  if(missing(map.function)) map.function = "KOSAMBI"
  if(missing(gen.file.path)) gen.file.path = p(full.path, "config.files/", project.name, ".gen")
  if(missing(ped.file.path)) ped.file.path = p(full.path, "config.files/", project.name, ".ped")
    
  lines1 = c("PLOIDY = 2", p("MAPFUNCTION = ", map.function), "MISSING = NA",
             p("CHROMFILE = ", full.path, "config.files/", project.name, ".chrom"),
             p("MAPFILE = ", full.path, "config.files/", project.name, ".map"),
             p("FOUNDERFILE = ", gen.file.path),
             p("PEDFILE = ", ped.file.path),
             p("OUTPUT = ", full.path, "simresults/", project.name, ".output"))
  write(lines1, file.path.to.write, ncolumns = 1)
}

make.pedsim.par.file.f3.selec = function(project.name, uniqueid, file.path.to.write, map.function){
  #makes a pedigreesim parameter file to be used in the f3 selection procedure
  #project.name - character sting, the name of the project
  #uniqueid - character string, a uniqueid to identify this parameter file (e.g. "DHF2840")
  #file.path.to.write - character string, the file path to write the parameter file
  #map.function - character string, the name of the map function to be used (defaults to Kosambi)
  if(missing(map.function)) map.function = "KOSAMBI"  
  lines1 = c("PLOIDY = 2", p("MAPFUNCTION = ", map.function), "MISSING = NA",
             p("CHROMFILE = ", full.path, "config.files/", project.name, ".chrom"),
             p("MAPFILE = ", full.path, "config.files/", project.name, ".map"),
             p("FOUNDERFILE = ", full.path, "config.files/", project.name, ".", uniqueid, "f3.selec.gen"),
             p("PEDFILE = ", full.path, "config.files/", project.name, ".", uniqueid, "f3.selec.ped"),
             p("OUTPUT = ", full.path, "simresults/", uniqueid, "_", project.name, ".output"))
  write(lines1, file.path.to.write, ncolumns = 1)
}

make.pedsim.config.files.and.run = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1){
  #makes pedigreesim configuration files and runs a simulation
  #args:
  #num.markers - integer, the number of markers to be used
  #num.individuals - integer, the number of individuals in the population
  #f.gen - integer, the highest filial generation of the population
  #cm.positions - numeric vector, the centimorgan positions of markers
  #write.path.folder - character string, the path of the pedigreesim installation
  #project.name - character string, the name of the project
  #map.function - character string, the map function to use (defaults to Kosambi)
  #run1 - boolean, indicates whether to run pedigreesim (if not, only generates config files)
  if(missing(run1)) run1 == F
  if(missing(map.function)) map.function = "KOSAMBI"
  make.founder.genotypes(num.markers, p(write.path.folder, "/config.files/", project.name, ".gen"))
  make.pedigree.file(f.gen, num.individuals, p(write.path.folder, "/config.files/", project.name, ".ped"))
  make.map.file(cm.positions, p(write.path.folder, "/config.files/", project.name, ".map"))
  make.chrom.file(max(cm.positions), p(write.path.folder, "/config.files/", project.name, ".chrom"))
  make.pedsim.par.file(project.name, p(write.path.folder, "/", project.name, ".par"), map.function = map.function)
  
  if(run1 == T){
    system(p("java -jar " full.path, "PedigreeSim.jar ", full.path, project.name, ".par"))  
  }  
}


make.pedsim.par.files.and.run.f3.selec = function(write.path.folder, project.name, uniqueid, map.function, run1){
  #args:
  #write.path.folder - character string, the path of the pedigreesim installation
  #project.name - character sting, the name of the project
  #uniqueid - character string, a uniqueid to identify this parameter file (e.g. "DHF2840")
  #map.function - character string, the name of the map function to be used (defaults to Kosambi)
  #run1 - boolean, indicates whether to run pedigreesim (if not, only generates config files)
  if(missing(run1)) run1 == F
  if(missing(map.function)) map.function = "KOSAMBI"
  make.pedsim.par.file.f3.selec(project.name, uniqueid, p(write.path.folder, "/", uniqueid, "_", project.name, ".par"), map.function = map.function)
  
  if(run1 == T){
    system(p("java -jar " full.path, "PedigreeSim.jar ", full.path, uniqueid, "_", project.name, ".par"))  
  }
}

convert.ped.sim.to.zygote = function(ped.sim){
  #converts a genotyping dataframe generated by pedigree sim from gametes to zygotes 
  #args:
  #ped.sim - dataframe, the genotyping dataframe produced by pedigreesim

  # ped.sim = read.delim(ped.sim.genotype.path, sep = "\t", stringsAsFactors = F)
  rownames(ped.sim) = ped.sim[, 1]
  #remove marker column
  ped.sim = ped.sim[, -1]
  ped.sim = convert.to.character.data.frame(ped.sim)
  
  #combine gametes
  ped.sim2 = lapply(seq(1, ncol(ped.sim), 2), function(q){
    combined.geno1 = unlist(Map(function(x, y){
      if(x == "A" & y == "A") z = "A"
      if(x == "A" & y == "B") z = "H"
      if(x == "B" & y == "A") z = "H"
      if(x == "B" & y == "B") z = "B"
      z
    }, ped.sim[, q], ped.sim[, (q + 1)]))  
    combined.geno1
  })
  
  ped.sim3 = do.call(cbind, ped.sim2)
  ped.sim3 = as.data.frame((ped.sim3))
  colnames(ped.sim3) = colnames(ped.sim)[seq(1, ncol(ped.sim), 2)]
  ped.sim3 = convert.to.character.data.frame(ped.sim3)
  ped.sim3
}

read.ped = function(ped.path) read.delim(ped.path, sep = "\t", stringsAsFactors = F)

f2.selection.procedure = function(project.name, selection.strength, selection.pos){
  #args:
  # selection.strength - integer, denominator of selection fraction
  # selection.pos - integer, locus at which to apply selection
  
  #if geno.10k doesn't already exists, load it and process it
  if(!exists("geno.10k")){
    geno.10k = read.ped(p(full.path, "10k.sim/axc1a.pop10k.output_genotypes.dat"))
    geno.10k = geno.10k[8:ncol(geno.10k)]
    geno.10k = geno.10k[, which(sapply(geno.10k, function(x) x[[selection.pos]] == "B"))]
  } 
  
  
  f2.genotypes = read.ped(p(full.path, "simresults/", project.name, ".output_genotypes.dat"))
  
  to.replace = names(which(sapply(f2.genotypes[4:ncol(f2.genotypes)], function(x){
    dice = sample(selection.strength, 1)
    if(dice == 1 & x[[selection.pos]] == "A"){
      return(T)
    } else {
      return(F)
    }
  })))
  
  replacements1 = geno.10k[, sample(ncol(geno.10k), length(to.replace))]
  
  colnames(replacements1) = to.replace
  
  f2.genotypes[, to.replace] = replacements1
  write.table(f2.genotypes,
              p(full.path, "simresults/", project.name, ".output_genotypes.dat_wselec.f2.dat"),
              sep = "\t",
              row.names = F,
              quote = F)
  return()
}

f3.selection.procedure = function(project.name, uniqueid, f.gen, selection.strength, selection.pos){
  #selection procedure for filial generations above F2
  #args:
  # output.genotypes.path - string, path to genotypes produced by PedSim
  # uniqueid - string, unique identifier for files produced by this procedure
  # f.gen - integer indicating filial generation
  # selection.strength - integer, denominator of selection pressure fraction
  # selection.pos - integer, locus at which  to apply selection
  
  # browser()
  
  f.gen2 = p("F", f.gen)
  prev.f.gen = p("F", (f.gen - 1))
  psr.path = p(full.path, "simresults/")
  output.genotypes.path = p(psr.path, project.name, ".output_genotypes.dat")
  prev.ped.file.path = p(full.path, "config.files/",
                         project.name, ".", f.gen, ".ped")
  
  
  
  ped1 = read.delim(output.genotypes.path, sep = "\t", stringsAsFactors = F)
  
  #first get individuals in which selection via gamete competition is possible (heterozygotes)
  ped.zygote = convert.ped.sim.to.zygote(ped1)
  
  ped.zygote2 = ped.zygote[, grep(f.gen2, colnames(ped.zygote))]
  
  
  selection.ind = which(sapply(ped.zygote2, function(x){
    dice1 = sample(selection.strength, 1)
    x[[selection.pos]] == "H" & dice1 == 1
  }))
  
  if(length(selection.ind) > 0){
    
    #### SELECTION PROCEDURE ####
    
    #get columns names of all gametes in which selection is possible
    snames = strsplit(names(selection.ind), "_")
    snames2 = unlist(lapply(snames, function(x){
      c(paste0(x[[1]], "_", x[[2]], "_", 1),
        paste0(x[[1]], "_", x[[2]], "_", 2))
    }))
    
    snames.prev.filial.gen = gsub(f.gen2, prev.f.gen, snames2)
    snames.prev.filial.gen.trunc = gsub("_[1-2]$", "", snames.prev.filial.gen)
    
    #prepare founder genotypes file
    ped.selec = ped1[, snames.prev.filial.gen]
    marker = ped1$marker
    ped.selec2 = as.data.frame(cbind(marker, ped.selec))
    
    #make new founder genotypes file with lines that are to undergo selection
    write.table(ped.selec2, p(full.path, "config.files/", project.name, ".", uniqueid, "f3.selec.gen"), sep = "\t", row.names = F, quote = F)
    
    #process pedigree file
    pedigree.file1 = read.delim(prev.ped.file.path, sep = "\t", stringsAsFactors = F)
    pedigree.file2 = pedigree.file1[which(pedigree.file1$Parent1 %in% snames.prev.filial.gen.trunc), ]
    pedigree.file3 = pedigree.file2[match(sort(rep(pedigree.file2$Parent1, 5)), pedigree.file2$Parent1), ]
    
    
    # if(nrow(pedigree.file3) == 0) browser()
    
    pedigree.file3 = reset.rownames(pedigree.file3)
    pedigree.file3$Name = paste0(pedigree.file3$Name, ".", rep(1:5, nrow(pedigree.file3) / 5))
    
    # pedigree.file3$Parent1 = gsub("_", ".", pedigree.file3$Parent1)
    # pedigree.file3$Parent2 = gsub("_", ".", pedigree.file3$Parent2)
    
    print("ped.parents")
    
    ped.parents = data.frame(unique(pedigree.file3$Parent1), rep(NA, length(unique(pedigree.file3$Parent1))), rep(NA, length(unique(pedigree.file3$Parent1))))
    # if(is.null(ncol(ped.parents))) browser()
    ped.parents = reset.colnames(ped.parents)
    # if(is.null(ncol(ped.parents))) browser()
    colnames(ped.parents) = c("Name", "Parent1", "Parent2")
    
    pedigree.file4 = rbind(ped.parents, pedigree.file3)
    
    #make new pedigree file with lines that are to undergo selection
    write.table(pedigree.file4, p(full.path, "config.files/", project.name, ".", uniqueid, "f3.selec.ped"), sep = "\t",
                row.names = F, quote = F)
    
    perform.again = T
    counter1 = 1
    while(perform.again == T){
      make.pedsim.par.files.and.run.f3.selec(full.path, project.name, uniqueid, "KOSAMBI", T)
      
      # if(counter1 > 6) browser()
      
      #run pedigreesim with new pedigree and founder files
      # browser()
      
      selection.test = read.ped(p(psr.path, uniqueid, "_", project.name, ".output_genotypes.dat"))
      
      selection.test.f3 = selection.test[, grep(f.gen2, colnames(selection.test))]
      
      selection.test.f3.2 = selection.test.f3[, which(sapply(selection.test.f3, function(x){
        x[[selection.pos]] == "B"
      }))]
      
      # if(is.null(colnames(selection.test.f3.2))) browser()
      
      
      
      colnames(selection.test.f3.2) = gsub("\\.[1-5]", "", colnames(selection.test.f3.2))
      trunc.colnames1 = gsub("_[1-2]$", "", colnames(selection.test.f3.2))
      #only need one of each gamete 
      selection.test.f3.3 = as.data.frame(selection.test.f3.2[, match(unique(trunc.colnames1), trunc.colnames1)])
      selection.test.f3.3 = convert.to.character.data.frame(selection.test.f3.3)
      colnames(selection.test.f3.3) = unique(trunc.colnames1)
      
      # q1 = gsub("_[1-2]$", "", colnames(selection.test.f3.3))
      
      #do we have gametes with a B at locus 200 for all of the lines we performed selection on ?
      perform.again = !all(unique(gsub("_[1-2]$", "", snames2)) %in% colnames(selection.test.f3.3))
      counter1 = counter1 + 1
    }
    
    #grab columns of ped1 (original genotype file) to replace
    original.colnames.to.replace = names(which(sapply(ped1[, snames2], function(x) x[[selection.pos]] == "A")))
    
    #check columns match, if so replace original genotyping data with new selected genotyping data.
    if(all(gsub("_[1-2]$", "", colnames(ped1[, original.colnames.to.replace])) == gsub("_[1-2]$", "", colnames(selection.test.f3.3)))){
      colnames(selection.test.f3.3) = colnames(ped1[, original.colnames.to.replace])
      #perform replacement
      ped1[, original.colnames.to.replace] = selection.test.f3.3
    }
    
    write.table(ped1, p(output.genotypes.path, "_wselec.f", f.gen, ".dat"), sep = "\t",
                row.names = F, quote = F)
    
    
    #### END OF SELECTION PROCEDURE ####
  } else {
    write.table(ped1, p(output.genotypes.path, "_wselec.f", f.gen, ".dat"), sep = "\t",
                row.names = F, quote = F)
  }
  
  
  
  # browser()
}

make.files.next.iteration = function(project.name, f.gen){
  #this function assumes one of the selection procedure functions has been run previously
  
  prev.fgen = f.gen - 1
  
  prev.geno = read.ped(p(full.path, "simresults/", 
                         project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".dat"))
  
  
  to.keep = colnames(prev.geno)[grep(p("F", prev.fgen), colnames(prev.geno))]
  
  parents1 = unique(gsub("_[1-2]$", "", to.keep))
  pedigree1 = data.frame(parents1, NA, NA)
  
  new.gen = gsub(p("F", prev.fgen), p("F", f.gen), parents1)
  
  ped2 = data.frame(new.gen, parents1, parents1)
  # if(is.null(colnames(pedigree1))) browser()
  colnames(pedigree1) = c("Name", "Parent1", "Parent2")
  # if(is.null(colnames(ped2))) browser()
  colnames(ped2) = c("Name", "Parent1", "Parent2")
  ped2 = ped2[grep(p("F", f.gen), ped2$Name), ]
  
  pedigree1 = convert.to.character.data.frame(pedigree1)
  ped2 = convert.to.character.data.frame(ped2)
  ped3 = rbind(pedigree1, ped2)
  # print("before row reset 1")
  # tryCatch(function(g) reset.rownames(ped3), error = function(x) browser())
  ped3 = reset.rownames(ped3)
  write.table(ped3,
              p(full.path, "config.files/", project.name, ".", f.gen, ".ped"),
              sep = "\t",
              row.names = F,
              quote = F)
  
  if(f.gen > 2){
    gen.file1 = read.ped(p(full.path, "simresults/",
                           project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".dat"))
    
    gen.file2 = gen.file1[, c(1, grep(p("F", prev.fgen), colnames(gen.file1)))]
    
    
    
    # newcolnames1 = strsplit(colnames(gen.file2), "_")
    # colnames(gen.file2) = sapply(newcolnames1, function(x) p(x[1], ".", x[2], "_", x[3]))
    # colnames(gen.file2)[1] = "marker"
    
    write.table(gen.file2,
                p(full.path, "simresults/",
                  project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".trunc.gen"),
                sep = "\t",
                row.names = F,
                quote = F)
    
    make.pedsim.par.file(project.name = project.name,
                         file.path.to.write = p(full.path, project.name, ".f", f.gen, ".par"),
                         map.function = "KOSAMBI",
                         gen.file.path = p(full.path, "simresults/", 
                                           project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".trunc.gen"),
                         ped.file.path = p(full.path, "config.files/", project.name, ".", f.gen, ".ped")
    )
    
  } else {
    make.pedsim.par.file(project.name = project.name,
                         file.path.to.write = p(full.path, project.name, ".f", f.gen, ".par"),
                         map.function = "KOSAMBI",
                         gen.file.path = p(full.path, "simresults/", 
                                           project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".dat"),
                         ped.file.path = p(full.path, "config.files/", project.name, ".", f.gen, ".ped")
    )
  }
  
  
  
  
  
  # browser()
  
  system(p("java -jar ", full.path, "PedigreeSim.jar ", full.path, project.name, ".f", f.gen, ".par"))
  
  
}

clean.files = function(){
  write.path.folder = full.path
  files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = "_founderalleles.dat")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = "_alleledose.dat")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = ".hsa")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = ".hsb")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".gen")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".map")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".ped")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
  
  files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".chrom")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
  
  files.to.rm = list.files(p(write.path.folder),  pattern = ".par")
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, x)))
}

complete.f5.selection.procedure = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder,
                                           project.name, map.function, run1, selection.strength, selection.pos){
  #args:
  #num.markers - integer, the number of markers to be used
  #num.individuals - integer, the number of individuals in the population
  #f.gen - integer, the highest filial generation of the population
  #cm.positions - numeric vector, the centimorgan positions of markers
  #write.path.folder - character string, the path of the pedigreesim installation
  #project.name - character sting, the name of the project  
  #map.function - character string, the map function to use (defaults to Kosambi)
  #run1 - boolean, indicates whether to run pedigreesim (if not, only generates config files)
  #selection.strength - integer, the strength of selection to be used (ranging from 0 - 1)
  #selection.pos - integer, the marker number against which to apply selection (e.g. 50 for the 50th marker)
  
  make.pedsim.config.files.and.run(num.markers, num.individuals, 2, cm.positions, write.path.folder,
                                   project.name, map.function, run1)
  f2.selection.procedure(project.name, selection.strength, selection.pos)
  
  if(f.gen > 2){
    lapply(3:f.gen, function(x){
      make.files.next.iteration(project.name, x)
      unique.code = paste0(LETTERS[sample(26, 5)], collapse = "", sample(9, 4))
      f3.selection.procedure(project.name, unique.code, x, selection.strength, selection.pos)
      print("done 1")
      
      files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = unique.code)
      lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
      
      files.to.rm2 = list.files(p(write.path.folder, "config.files/"), pattern = unique.code)
      lapply(files.to.rm2, function(x) file.remove(p(write.path.folder, "config.files/", x)))
      
      files.to.rm3 = list.files(p(write.path.folder), pattern = unique.code)
      lapply(files.to.rm3, function(x) file.remove(p(write.path.folder, x)))
      
    })
  }
  
  
  
  
  files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = p(project.name, "\\.output_genotypes\\.dat"))
  files.to.rm = files.to.rm[-grep(p("wselec", "\\.f", f.gen, "\\.dat"), files.to.rm)]
  lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
}


complete.fgen.no.selec.procedure = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1){
  #args:
  #num.markers - integer, the number of markers to be used
  #num.individuals - integer, the number of individuals in the population
  #f.gen - integer, the highest filial generation of the population
  #cm.positions - numeric vector, the centimorgan positions of markers
  #write.path.folder - character string, the path of the pedigreesim installation
  #project.name - character sting, the name of the project  
  #map.function - character string, the map function to use (defaults to Kosambi)
  #run1 - boolean, indicates whether to run pedigreesim (if not, only generates config files)
  
  make.pedsim.config.files.and.run(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1)
}


#### ANALYSIS FUNCTIONS 1 ####


convert.recomb.to.seg = function(recomb.data){
  #args: 
  # recomb.data - dataframe containing genotyping information
  seg.data = unlist(lapply(recomb.data, function(x){
    # browser()
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    a / (a + b)
  })) 
  seg.data
}

convert.recomb.to.seg2 = function(recomb.data, p.value1){
  #args:
  # recomb.data - a dataframe containing genotype information
  # p.value1 - boolean flag indicating whether to only return p-values from chi-square tests
  
  if(missing(p.value1)) p.value1 = F
  
  seg.data = lapply(recomb.data, function(x){
    # browser()
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    chisq.test(c(a, b))
  })
 
  if(p.value1 == T){
    return(unlist(lapply(seg.data, function(x) x[[3]])))
  } else {
    return(seg.data)
  }
  
}

convert.recomb.to.seg.f2 = function(recomb.data, p.value1){
  #args:
  # recomb.data - a dataframe containing genotype information
  # p.value1 - boolean flag indicating whether to only return p-values from chi-square tests
  
  if(missing(p.value1)) p.value1 = F
  
  seg.data = lapply(recomb.data, function(x){
    # browser()
    a = length(which(x == "A"))
    h = length(which(x == "H"))
    b = length(which(x == "B"))
    chisq.test(c(a, h, b), p = c(1/4, 1/2, 1/4))
  })
  
  if(p.value1 == T){
    return(unlist(lapply(seg.data, function(x) x[[3]])))
  } else {
    return(seg.data)
  }
}


# lapply(selected1, function(x){
#   # browser()
#   a = length(which(x == "A"))
#   b = length(which(x == "B"))
#   chisq.test(c(a, b))
# })



check.for.distorted.markers = function(recomb.data, thresholdlevel, F2){
  #returns the coordinates of markers that show significant seg. dist. according to chi-square test
  #args:
  # recomb.data - a dataframe containing genotyping data
  # thresholdlevel - an integer specifying the level of significance, defaults to 0.05
  if(missing(thresholdlevel)) thresholdlevel = 0.05
  if(missing(F2)) F2 = F
  
  if(F2 == T){
    g = convert.recomb.to.seg.f2(recomb.data)
  } else {
    g = convert.recomb.to.seg2(recomb.data)
  }
  
  g2 = which(unname(unlist(lapply(g, function(x) x[3]))) < thresholdlevel)
  return(g2)
}

check.for.distorted.markers.fdr = function(recomb.data, thresholdlevel, method1, F2){
  #returns the coordinates of markers that show significant seg. dist. according to chi-square test
  #args:
  # recomb.data - a dataframe containing genotyping data
  # thresholdlevel - an integer specifying the level of significance, defaults to 0.05
  if(missing(thresholdlevel)) thresholdlevel = 0.05
  if(missing(method1)) method1 = "BH"
  if(missing(F2)) F2 = F
  
  if(F2 == T){
    g = convert.recomb.to.seg.f2(recomb.data)
  } else {
    g = convert.recomb.to.seg2(recomb.data)
  }
  
  g2 = which(p.adjust(unname(unlist(lapply(g, function(x) x[3]))), method1) < thresholdlevel)
  return(g2)
}

prepare.gg.data = function(recomb.data){
  new.data = data.frame(convert.recomb.to.seg(recomb.data), "", stringsAsFactors = F)
  new.data[, 2][check.for.distorted.markers(recomb.data)] = "T"
  new.data[, 2][check.for.distorted.markers(recomb.data, 0.01)] = "T1"
  new.data[, 2][check.for.distorted.markers(recomb.data, 0.001)] = "T2"
  colnames(new.data) = c("seg.ratio", "sig")
  # new.data$chromo = unique(recomb.data[1, 3:ncol(recomb.data)])
  new.data
}

make.plots = function(list.of.recomb.data, plot.titles){
  #args:
  # list.of.recomb.data - a list of genotyping dataframes
  count1 = make.counter()
  if(missing(plot.titles)) plot.titles = 1:length(list.of.recomb.data)

  lapply(list.of.recomb.data, function(x){
    g = prepare.gg.data(x)
    # browser()
    ggplot(g, aes(x = 1:nrow(g), y = seg.ratio, group = sig, color = sig)) + geom_point() + 
      geom_hline(yintercept = 0.5) + ggtitle(plot.titles[count1()])# + coord_cartesian(ylim = c(0.5, 1.6))
  })
}

find.recomb.w.seg.dist = function(list.of.recomb.data){
  unlist(lapply(list.of.recomb.data, function(x){
    g = which(unlist(lapply(convert.recomb.to.seg2(x), function(x2) x2[3])) < 0.05)
    if(length(g) > 0) return(T)
    if(length(g) == 0) return(F)
  }))
}

#### ANALYSIS FUNCTIONS 2 ####

#functions to evaluate lists of genotype data simulations

number.within.10 = function(x, selection.position){
  #grab number of simulations in which the peak of segregation distortion is within 10 markers of the point of selection
  #args:
  #x - a list of genotyping dataframes 

  if(missing(selection.position)) selection.position = 200

  g = unlist(lapply(x, function(y){
    g = convert.recomb.to.seg(y)
    g1 = abs(g - 0.5)
    mean(which(g1 == max(g1)))
  }))

  length(which(g > (selection.position - 10) & g < (selection.position + 10)))
}

number.dist.markers = function(y, threshold, F2){
  #check number of simulations with distorted markers
  #args:
  #y - a list of genotyping dataframes 
  #threshold - integer, p-value threshold, defaults to 0.05
  if(missing(threshold)) threshold = 0.05
  if(missing(F2)) F2 = F
  
  length(which(unlist(lapply(y, function(x){
  g = check.for.distorted.markers(x, threshold, F2 = F2)
  if(length(g) > 0) return(T)
  return(F)
  }))))
}


number.dist.markers.fdr = function(y, threshold, method1, F2){
  #check number of simulations with distorted markers
  #args:
  #y - a list of genotyping dataframes 
  #threshold - integer, p-value threshold, defaults to 0.05
  if(missing(threshold)) threshold = 0.05
  if(missing(method1)) method1 = "BH"
  if(missing(F2)) F2 = F
  
  length(which(unlist(lapply(y, function(x){
  g = check.for.distorted.markers.fdr(x, threshold, method1, F2 = F2)
  if(length(g) > 0) return(T)
  return(F)
  }))))
}

mean.sd.num.distorted = function(y, threshold, F2){
  #grab the mean and standard deviations of number of distorted markers for each group of 1000 simulations
  #args:
  #y - a list of genotyping dataframes 

  if(missing(threshold)) threshold = 0.05
  if(missing(F2)) F2 = F
  q2 = unlist(lapply(y, function(x){
    g = check.for.distorted.markers(x, threshold, F2 = F2)
    length(g)
  }))
  
  g1 = list(mean(q2[which(q2 != 0)]), sd(q2[which(q2 != 0)]))
  names(g1) = c("mean", "sd")
  g1
}

mean.sd.num.distorted.fdr = function(y, threshold, method1, F2){
  #grab the mean and standard deviations of number of distorted markers for each group of 1000 simulations
  #args:
  #y - a list of genotyping dataframes 
  if(missing(threshold)) threshold = 0.05
  if(missing(method1)) method1 = "BH"
  if(missing(F2)) F2 = F

  q2 = unlist(lapply(y, function(x){
    g = check.for.distorted.markers.fdr(x, threshold, method1, F2 = F2)
    length(g)
  }))
  
  g1 = list(mean(q2[which(q2 != 0)]), sd(q2[which(q2 != 0)]))
  names(g1) = c("mean", "sd")
  g1
}


peak.dist = function(x){
  #what is the mean peak of distortion and what is the sd of the peak of distortion between simulations
  #args:
  #x - a list of genotyping dataframes 
  g = unlist(lapply(x, function(y){
    g = convert.recomb.to.seg(y)
    g1 = abs(g - 0.5)
    mean(which(g1 == max(g1)))
  }))

  g
}

check.magnitude.of.distortion = function(x){
  #what is the magnitude of distortion in the simulations at the peak
  #args:
  #x - a list of genotyping dataframes 
  g = unlist(lapply(x, function(y){
    g = convert.recomb.to.seg(y)
    g1 = abs(g - 0.5)
    max(g1)
  }))

  g
}



#ns = number.sims.w.dist.markers
#mn = mean.number.distorted.markers.in.sims.w.distortion
#sdn = sd.number.distorted.markers.in.sims.w.distortion
#nsp = number.of.simulations.where.peak.of.distortion.is.within.10.markers.of.selection.locus
#mmp = mean.magnitude.of.peak.distortion
#mmpsd = sd.magnitude.of.peak.distortion
#maxmag = max.magnitude.of.peak.distortion
#minmag = min.magnitude.of.peak.distortion


perform.analysis = function(geno.list, selection, selection.position, selection.strength, recombination.position, F2, fgen){
  #args:
  # geno.list - a list of genotyping dataframes
  # selection - boolean flag - was selection used in this simulation?
  # selection.position - integer, the position of selection
  # selection.strength - integer, the strength of selection
  # recombination.position - string, the recombination profile (if any) used
  # F2 - boolean flag indicating whether to use f2 type test of segregation distortion, otherwise assumes 1:1 ratio of homozygotes as distribution (ignores hets)
  # fgen - integer specifying the filial generation of the simulation
  

  if(missing(selection)) selection = "-"
  if(missing(selection.position)) selection.position = "-"
  if(missing(selection.strength)) selection.strength = "-"  
  
  if(missing(F2)) F2 = F
  if(missing(recombination.position)){
    return("need to enter recombination.position")
  } else {
  
  number.simulations = length(geno.list)
  pop.size = nrow(geno.list[[1]])
  num.markers = ncol(geno.list[[1]])  
  ns0.05 = number.dist.markers(geno.list, 0.05, F2 = F2)  
  ns0.01 = number.dist.markers(geno.list, 0.01, F2 = F2)
  ns0.001 = number.dist.markers(geno.list, 0.001, F2 = F2)
  ns.fdr.0.05 = number.dist.markers.fdr(geno.list, 0.05, F2 = F2)
  ns.bon.0.05 = number.dist.markers.fdr(geno.list, 0.05, "bonferroni", F2 = F2)
  mnsd = mean.sd.num.distorted(geno.list, 0.05)
  mn.0.05 = mnsd[[1]]
  sdn.0.05 = mnsd[[2]]
  mnsd = mean.sd.num.distorted(geno.list, 0.01)
  mn.0.01 = mnsd[[1]]
  sdn.0.01 = mnsd[[2]]
  mnsd = mean.sd.num.distorted(geno.list, 0.001)
  mn.0.001 = mnsd[[1]]
  sdn.0.001 = mnsd[[2]]
  mnsd = mean.sd.num.distorted.fdr(geno.list, 0.05)
  mn.fdr.0.05 = mnsd[[1]]
  sdn.fdr.0.05 = mnsd[[2]]
  mnsd = mean.sd.num.distorted.fdr(geno.list, 0.05, "bonferroni")

  mn.bon.0.05 = mnsd[[1]]
  sdn.bon.0.05 = mnsd[[2]]
  if(selection == "-" | selection.position == "-"){
    nsp = "-"
  } else {
    nsp = number.within.10(geno.list, selection.position)  
  }
  
  mag.distortion = check.magnitude.of.distortion(geno.list)
  mmp = mean(mag.distortion)
  mmpsd = sd(mag.distortion)
  maxmag = max(mag.distortion)
  minmag = min(mag.distortion)

  summarydf1 = newdf(c("R.object.name", "pop.size", "fgen", "selection", "selection.position", "selection.strength", "recombination.position", "number.simulations", "number.markers", "ns.0.05", "ns.0.01", "ns.0.001", "ns.fdr.0.05", "ns.bon.0.05", "mn.0.05", "mn.0.01", "mn.0.001", "mn.fdr.0.05", "mn.bon.0.05", "sdn.0.05", "sdn.0.01", "sdn.0.001", "sdn.fdr.0.05", "sdn.bon.0.05", "nsp", "mmp", "mmpsd", "maxmag", "minmag"))  


  summarydf1[1, ] = c(deparse(substitute(geno.list)), pop.size, fgen, selection, selection.position, selection.strength, recombination.position, number.simulations, num.markers, ns0.05, ns0.01, ns0.001, ns.fdr.0.05, ns.bon.0.05, mn.0.05, mn.0.01, mn.0.001, mn.fdr.0.05, mn.bon.0.05, sdn.0.05, sdn.0.01, sdn.0.001, sdn.fdr.0.05, sdn.bon.0.05, nsp, mmp, mmpsd, maxmag, minmag)

  summarydf1

  # summarydf1$R.object.name = 
  # summarydf1$pop.size = pop.size
  # summarydf1$selection = selection
  # summarydf1$selection.position = selection.position
  # summarydf1$selection.strength = selection.strength
  # summarydf1$recombination.position = recombination.position
  # summarydf1$

  }

  
}

