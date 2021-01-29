




IsoCorTables <- function(PathToCSV_file) {
  
  WT1 <- read.csv(file = PathToCSV_file, header = T)
  
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
    
    
    library(reshape2)
    
    test <- melt(data = ToMelt_mat, varnames = "isotopologue", value.name = "area")
    
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
  
  return(out_mat)
  
}
