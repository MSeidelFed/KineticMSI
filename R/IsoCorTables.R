#' A function to produce the input files for Python-based correction of natural isotopic abundances (NIA)
#'
#' This function allows you to produce the input files needed for the Python tool [IsoCor](https://pypi.org/project/IsoCor/) that corrects isotopologue envelopes for NIA. The function takes a single csv file with the right format and transforms it. Subsequently, the function generates a matrix to the R environment that can be exported in any desired format to serve as input for IsoCor.
#' @param PathToCSV directory where the input file is stored.
#' @keywords KineticMSI NIA correction IsoCor
#' @export
#' @examples
#'
#' ...


IsoCorTables <- function(PathToCSV) {

  WT1 <- read.csv(file = PathToCSV, header = T)
  
  dir.create("IsoCorInput")
  
  ## metabolites table

  Metabolites <- as.factor(WT1[,"Molecular.formula"])

  name <- c()

  formula <- c()

  charge <- c()

  for (i in 1:length(levels(Metabolites))) {
  
    test <- length(strsplit(x = levels(Metabolites)[i], split = "")[[1]])
  
    if (test != 0) {
    
      index_met_i <- which(WT1[,"Molecular.formula"] == levels(Metabolites)[i])
    
      if (length(index_met_i) > 1) {
      
        cat("...\n")
        stop("ERROR: Duplicated features in the input table...")
      
      } else {
      
        name <- c(name, as.character(WT1[index_met_i ,"Metabolite"]))
    
        formula <- c(formula, as.character(WT1[index_met_i ,"Molecular.formula"]))
    
        charge <- c(charge, -1) 
      
      }
    }
  }

  Metabolites_out_table <- cbind(name, formula, charge)

  write.table(x = Metabolites_out_table, file = paste0("IsoCorInput","/Metabolites.dat"), quote = F, row.names = F, sep = "\t")
  
  ## data table

  unique_lipids <- as.character(unique(WT1[,"Metabolite"]))

  IsoCor_Mats <- matrix(NA, nrow = 1, ncol = 6)

  for (i in 1:length(unique_lipids)) {

    lipid_group1 <- WT1[which(WT1[,"Metabolite"] == unique_lipids[i]),]

    formula_subgroup <- as.character(lipid_group1[1,"Molecular.formula"])

    formula_subgroup_charachters <- unlist(strsplit(x = formula_subgroup, ""))

    H_position <- grep(pattern = "H", x = formula_subgroup_charachters)

    alpha_positions <- grep("[[:alpha:]]",formula_subgroup_charachters)

    nextTO_H_position <- alpha_positions[which(alpha_positions == H_position) + 1]

    H_digits <- (H_position + 1) : (H_position + nextTO_H_position - H_position - length(H_position))

    H_number <- paste0(formula_subgroup_charachters[H_digits], collapse = "")

    isotopologue = 0:H_number

    ## mat

    ToMelt_mat <- matrix(data = NA,
                         nrow = length(isotopologue),
                         ncol = dim(lipid_group1)[2] - 4)

    rownames(ToMelt_mat) <- isotopologue

    ToMelt_mat[1:dim(lipid_group1)[1], 1:(dim(lipid_group1)[2] - 4)] <- as.matrix(lipid_group1[,5:dim(lipid_group1)[2]])

    test <- reshape2::melt(data = ToMelt_mat, varnames = "isotopologue", value.name = "area")

    metabolite = rep(as.character(unique(lipid_group1[,"Metabolite"])), dim(test)[1])

    IsoCor_Mat <- cbind(sample = test[,2],
                        metabolite = metabolite,
                        derivative = rep(NA, dim(test)[1]),
                        isotopologue = test[,1],
                        area = test[,3],
                        resolution = rep(70000, dim(test)[1]))

    IsoCor_Mats <- rbind(IsoCor_Mats, IsoCor_Mat)
  }

  out_mat <- IsoCor_Mats[2:dim(IsoCor_Mats)[1],]

  replacement_AREA = vector(length = dim(out_mat)[1])

  replacement_DER = vector(length = dim(out_mat)[1])

  for (i in 1:dim(out_mat)[1]) {

    if (is.na(out_mat[i,"area"])) {runner = 0} else {runner = out_mat[i,"area"]}

    replacement_AREA[i] = runner

    if (is.na(out_mat[i,"derivative"])) {runner = "NULL"} else {runner = out_mat[i,"area"]}

    replacement_DER[i] = runner
  }

  out_mat[,5] = replacement_AREA
  out_mat[,3] = replacement_DER
  
  write.table(out_mat, file = paste0("IsoCorInput","/Data_example.tsv"), row.names = F, quote = F, sep = "\t")

  return(list(Metabolites_out_table, out_mat))

}
