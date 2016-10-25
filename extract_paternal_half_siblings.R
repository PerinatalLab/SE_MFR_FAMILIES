# Dominika Modzelewska
# 21 October 2016

# loading dataset (merged mfr and biofaldrar)
load("/home/dom/Desktop/mfr/dt.RData")
load("/home/dom/Desktop/mfr/paternal_half_siblings.Rda")

#############
############# CLEANING DATASET
# maternal id inocsistences between mfr and biofaldrar
# missing maternal id
# missing paternal id
# duplicated maternal and child entries
# no delation of missing child's id (since all missing are due death born)

dat = dt[, c("lpnr_BARN", "lpnr_mor", "LopnrMor", "LopnrFar", "AR", "GRDBS")]

#### CHECKING THE INCONSISTENCES IN THE MATERNAL ID
incon = sum(dat$lpnr_mor==dat$LopnrMor, na.rm=T)
cat("number of inconsistent maternal id's, between mfr and biofod:", nrow(dat)-incon) # 200

# checking why the inconsistences are observed: 1. MISSING VALUES
a = sum(is.na(dat$lpnr_mor[which(dat$lpnr_mor!=dat$LopnrMor)]) | dat$lpnr_mor[which(dat$lpnr_mor!=dat$LopnrMor)]=="0")
b = sum(is.na(dat$LopnrMor[which(dat$lpnr_mor!=dat$LopnrMor)]) | dat$LopnrMor[which(dat$lpnr_mor!=dat$LopnrMor)]=="0")
cat("NB: some missing values are indicated as 0:",
    a,  "of that inconsistences is due to missing maternal's id value in mfr ",
    b, "of that inconsistences is due to missing maternal's id value in biofod")

# those with missing id values will be replaced by maternal id from MFR
dat$LopnrMor[which((dat$lpnr_mor!=dat$LopnrMor) & dat$LopnrMor==0)] = "NA"
dat$LopnrMor[which(dat$LopnrMor=="NA")] = dat$lpnr_mor[which(dat$LopnrMor=="NA")]
cat(sum(dat$lpnr_mor!=dat$LopnrMor, na.rm=T), "inconsistences remain in the dataset")

# checking why the inconsistences are observed: 2. DIGIT PROBLEMS
# check in how many digits they are different
data.frame(dat$lpnr_mor[which(dat$lpnr_mor!=dat$LopnrMor)], dat$LopnrMor[which(dat$lpnr_mor!=dat$LopnrMor)]) # done "by eye", check it

#eliminating the inconsistences
bad_rows = which(dat$lpnr_mor!=dat$LopnrMor) # you can't do dat[dat$lpnr_mor!=dat$LopnrMor, ] as then we create a datest and we do not specify the numbers o # of the row to be delated, what is needed in next step
dat = dat[-bad_rows, ]

# removing unwanted column LopnrMor
dat = dat[, !names(dat)=="LopnrMor"]

# removing missing paternal id's 
bad_rows = which(dat$LopnrFar==0 | is.na(dat$LopnrFar))
dat = dat[-bad_rows, ]
cat("missing paternal id's:", length(bad_rows)) # 31200

# remove duplicated child and mother entries
tmp0 = dat[which( (!is.na(dat$lpnr_BARN)) & (!is.na(dat$lpnr_mor))),] # both mom and child have PersNumb

#  mother and child pair was entered twice:
morbar_unqcode = paste(tmp0$lpnr_BARN,tmp0$lpnr_mor,sep="_") # unique code for mother-child pair
dpl_codes = unique(morbar_unqcode[which(duplicated(morbar_unqcode))]) # codes that are duplicated
dpl_rows = which(morbar_unqcode %in% dpl_codes) # rows that contain duplicated pairs
unq_rows = which(! morbar_unqcode %in% dpl_codes) # rows that contain unique pairs

tmp0_dpl = tmp0[dpl_rows,] # duplicated mother-child pairs
tmp0_unq = tmp0[unq_rows,] # unique mother-child pairs, n = 4 052 959

dat = tmp0_unq # NOTE: at this moment we do not care about which row is preferred to remain in dataset

cat("remained number of individuals:" , nrow(dat)) # 4021525


#################
################# FINDING FAMILY STRUCTURES

colnames(dat)[1:3] = c("kid", "mom", "dad")

#######
####### FINDING CHLDREN WHO HAVE FULL SIBLING (they might also have half siblings, but in the created dataset it is not seen)

dat$parent_code = paste(dat$mom, dat$dad, sep="_") # all parent's codes in the dataset
parent_indic = unique(dat$parent_code[duplicated(dat$parent_code)]) # the codes of parents who have more than one child

dat$fsib = 0 # indicator if child who has a full siblings match
dat$fsib[which(dat$parent_code %in% parent_indic)]=1

sib_rows = which(dat$parent_code %in% parent_indic)
sib = dat[sib_rows, ] # dataset containing only children which have full siblings

#####
##### by the variable df$count there are indicated full siblings groups

# let's create dataset containg the parent_code and column with their children and GA
parent_indic = parent_indic[1:10] # to check

rows=NULL
datalist = list()
for(i in 1:length(parent_indic)){
  j = parent_indic[i]
  rows = grep(j, dat$parent_code) # pokazuje numery wierszy ciaz w ktorych byl ten sam ojciec i matka
  fsib = dat[rows, ]
  fsib$count = i
  datalist[[i]] = fsib
}
df = do.call(rbind, datalist)

cat("in a dataset we hava", mas(df$count), "families which have full siblings")
# to see how many full siblings each family has
table(no_family=df$count)
hist(df$count, breaks=100, freq=F) # to check which value (number of FS) occurs the most often in a dataset


######
######  FINDING PATERNAL HALF SIBLINGS

# find children who have the same father
father_indic = unique(dat$dad[duplicated(dat$dad)]) # the ids of fathers who have more than 1 child
cat("there are", length(father_indic), "fathers, who have more than one child") # 1342639

df_paternal_sib = dat[dat$dad %in% father_indic, ] # create a dataset which contains the multiple pregnacies for father (in other words, there are no fathers who have only one child)
head(df_paternal_sib)

# father_indic = father_indic[1:100]

# creating the matches of all paternal half siblings
ptm <- proc.time()
no_cheat= NULL
cb=list()

for(i in 1:length(father_indic)){
  datalist=list()
  comb_list=list()
  j = father_indic[i] # the id of father which we are checking
  parents = unique(df_paternal_sib$parent_code[which(df_paternal_sib$dad==j)])
  if (length(parents)==1){
    no_cheat = c(no_cheat, j)
    
  } else {
    
    for(z in 1:length(parents)){
      v = parents[z] # indicate which parther he had in a loop
      fisrt_child_id = df_paternal_sib$kid[which(df_paternal_sib$parent_code==v)] # getting know who is full who half sib for a given guy
      datalist[[z]] = fisrt_child_id # kids who are within a list are full siblings, the children from different lists are paternal half siblings
    }
    combination = combn(datalist, 2, simplify = F)
    
    for( q in 1:length(combination)){
      comb_list[[q]]= do.call(expand.grid, combination[[q]])
    }
  } 
  cb[[i]] = do.call(rbind, comb_list)
}
head(cb)
a = proc.time() - ptm
cat("time for", length(father_indic), "is:", a, "user, system, elapsed") # 11.002 0.163 11.172

#### create and save a dataset of paternal half siblings
paternal_half_siblings = do.call(rbind, cb) # create a dataset of paternal half siblings
colnames(paternal_half_siblings) = c("first_child", "his_phs") # renaming the columns: 1-first child, 2-his matched paternal half sibling
save("paternal_half_siblings", file="/home/dom/Desktop/mfr/paternal_half_siblings.RData")
save("paternal_half_siblings", file="/home/dom/Desktop/mfr/paternal_half_siblings.Rda")

### create and save a dataset of fathers who had only one sexual partner, and did not produce half siblings
no_half_sib= data.frame(no_cheat)
colnames(no_half_sib) = "father_id"
save("no_half_sib", file="/home/dom/Desktop/mfr/no_half_sib.Rdata" )

#### add to a dataset desired variables (kid is a key variable)
# we can merge easily the father's id to our paternal_half_siblings dataset since children from both columns are from the same father
# desired_data = df_paternal_sib[, c("dad", "kid")] 
# colnames(desired_data)[1]="first_child"
# data = merge(big_cb, sm_dat, "first_child")
colnames(df_paternal_sib)
# but to make sure that we do not have any error, and the two kids in paternal_half_siblings dataset are from the same father, we take the same approach as for mother:
desired_data = df_paternal_sib[ , c("dad", "mom", "kid", "AR", "GRDBS")]
colnames(desired_data) = c("dad_fist", "mom_first", "first_child", "ar_first", "ga_first") # we first merge based on the first kid
data_check = merge(desired_data, paternal_half_siblings, by= "first_child")

desired_data = df_paternal_sib[ , c("dad", "mom", "kid", "AR", "GRDBS")]
colnames(desired_data) = c("dad_hs", "mom_hs", "his_phs", "ar_hs", "ga_hs") # next we  merge based on the matched HS
data_check = merge(desired_data, data_check, by="his_phs")

sum((data_check$dad_first - data_check$data_hs)!=0) # fathers shouls be alwasy the same, thus the sum must be equal to 0
sum((data_check$mom_first - data_check$mom_hs)==0) # mothers should always be different, thus the sum must be equal to 0


# check if there are repeated pairs


########
######## ESTIMATING HERITABILITY BASED ON THE INTRACLASS §CORRELATION OF GESTATIONAL AGE BETWEEN PATERNAL HALF SIBLINGS

# STEP 1: estimating the correlation between paternal half siblings by:
  # 1.1 multilevel analysis
  # 1.2 taking randomly the half sibling pair from a one father, many times, calculate corr each time, and at the end take a mean of those values (?)
# STEP 2: estimating the heritability by multiplying corr coefficient by 4 (arrive from theory of heritability)


#######
####### CHECKING IF THE HERITABILITY ESTIMATES DIFFER ACCROSS THE YEARS AND AREAS OF SWEDEN















