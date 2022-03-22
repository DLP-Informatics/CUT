

LeftNuc <- gene_arr[iter-1] #left neighbor nucleotide
MiddleNuc <- gene_arr[iter] #center nucleotide
RightNuc <- gene_arr[iter+1] #right neighbor nucelotide

LeftNuc <- toupper(LeftNuc) #left neighbor nucleotide
MiddleNuc <- toupper(MiddleNuc) #center nucleotide
RightNuc <- toupper(RightNuc) #right neighbor nucelotide

# DiagTriplet <- paste(LeftNuc, MiddleNuc, RightNuc, sep = "")
# DiagArr <- c(DiagArr, DiagTriplet



rowswitch <- paste(LeftNuc, RightNuc, sep = "") #row matrix assignment operator
colswitch <- MiddleNuc #column matrix assignment operator

#null reference checker. checking the 5mer block for null references
if(is.na(LeftNuc))
{
  j <- j+3
  next
}
if(is.na(MiddleNuc))
{
  j <- j+3
  next
}
if(is.na(RightNuc))
{
  j <- j+3
  next
}



#checks the nucleotide triplets for non-ambiguous entries, skips if not atcg
if(LeftNuc != 'A' & LeftNuc != 'T' & LeftNuc != 'C' & LeftNuc != 'G')
{
  j <- j+3
  next
}
if(MiddleNuc != 'A' & MiddleNuc != 'T' & MiddleNuc != 'C' & MiddleNuc != 'G')
{
  j <- j+3
  next
}
if(RightNuc != 'A' & RightNuc != 'T' & RightNuc != 'C' & RightNuc != 'G')
{
  j <- j+3
  next
}



#Row matrix assignment operator for 3mer
rowiter <- switch(rowswitch,
                  "TT" = 1,
                  "TG" = 2,
                  "TC" = 3,
                  "TA" = 4,
                  "GT" = 5,
                  "GG" = 6,
                  "GC" = 7,
                  "GA" = 8,
                  "CT" = 9,
                  "CG" = 10,
                  "CC" = 11,
                  "CA" = 12,
                  "AT" = 13,
                  "AG" = 14,
                  "AC" = 15,
                  "AA" = 16
)
#column assignment operator
coliter <- switch(colswitch,
                  "T" = 1,
                  "G" = 2,
                  "C" = 3,
                  "A" = 4
)

if(orientation == '+')
{
  increment <- as.integer(CU_Forward[rowiter,coliter])
  increment <- increment + 1
  CU_Forward[rowiter, coliter] <- increment 
  
}
if(orientation == '-')
{
  increment <- as.integer(CU_Reverse[rowiter,coliter])
  increment <- increment + 1
  CU_Reverse[rowiter, coliter] <- increment 
  
}


