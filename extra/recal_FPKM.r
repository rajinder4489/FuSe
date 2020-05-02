rm(list=ls())

#The colnames and rowwnames should be same foe both expr and effective_len file.

#Read the normalized read counts file
expr <- read.delim(file="", sep="\t", stringsAsFactors = F)

#Make sure that the transcript ids are rownames
expr_sort <- expr[ order(row.names(expr)), ]

#Read the effective length file, it can be obtained from one of the output files from the tool used to align the reads to genome. Refer the manual of the tool used.
effective_len <- read.delim(file="isoform_effective_len.txt", sep="\t", stringsAsFactors = F)

#Make sure that the transcript ids are rownames
effective_len_sort <- effective_len[ order(row.names(effective_len)), ]


if(setequal(colnames(expr_sort), colnames(effective_len_sort)) & setequal(rownames(expr_sort), rownames(effective_len_sort)))
{
	norm_fpkm <- expr*10^9/colSums(expr)/effective_len	#normalized read count used for calculating FPKM
}


is.na(norm_fpkm)<-sapply(norm_fpkm, is.infinite)		#converts all Inf to NA

norm_fpkm[is.na(norm_fpkm)] <- 0		#this is required to fix the cases where 0 divided by 0 happened resulting in NA and Inf converted to NA in the last step

write.table(norm_fpkm, file="", sep="\t")
