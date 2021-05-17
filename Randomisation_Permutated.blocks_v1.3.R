library(stringr)

# Setting the seed for the random sampling. Actual seed used to randomise is held by the unblinded statistician
set.seed(56482)

# NB for samples greater than 99, the function (allocation) below needs tweaking to correctly prepend the right number of zeros to the Sequence number (to be patched up)
 
n<-104  			# Total sample size
a<-c(4,6,8)			# Choose block sizes, it works for size 12 too
b<-c(4,8)			
s<-0				#initialise sum
my_vector = c()		#initialise vector

n_arms <- 2			# Number of arms
# Enter preferred arm labels here:
arm_labels <- c("Ketamine", "Midazolam")

# Loop to choose the appropriate last block size not to exceed the total sample size.

for (i in 1:n){
 if (s<n-8) {				
	i <- sample(a,1)
	s<-sum(i+s)
	my_vector=c(my_vector,i)
 } 
else if (s==n-8){				#if residual is 8, randomly allocate 4 or 8 block
	i <- sample(b,1)
	s<-sum(i+s)
	my_vector=c(my_vector,i)
 }
 else if (s==n-6){			#if residual is 6, allocate block size 6
	i=6
	s<-sum(i+s)
	my_vector=c(my_vector,i)
 } 
else if (s==n-4) {			#if residual is 4, allocate block size 4

	i=4
	s<-sum(i+s)
	my_vector=c(my_vector,i)
}}
s					# print the sum of block sizes			
block_sizes <- my_vector	#sequence of blocks


	block_number <- length(block_sizes)						# Number of blocks

	# Permute the blocks
	p_block_sizes <- sample(block_sizes, replace=F)


	# Initialise the data frame with the Randomisation ID (RID) in format 00x00yy
	# where x is the site number (78 in site 1, 26 in site 2) and yy is a patient sequence.
	# The allocation (arm) column is left blank for now:
	alloc <- data.frame(RID=paste(rep("001",n),str_pad(c(1:n),4,pad="0"), sep=""), arm=NA)

	current_pat <- 1
	# For each block:
	for (i in 1:block_number) {
	 # Set the allocation for the 
		alloc[current_pat:(current_pat+p_block_sizes[i]-1),]$arm <- c(sample(rep(seq(1:n_arms), (p_block_sizes[i]/n_arms)), p_block_sizes[i], replace=F))
		current_pat <- current_pat+p_block_sizes[i]
	}

	for (j in 1:n_arms) {
		alloc$arm<-replace(alloc$arm, alloc$arm==j, arm_labels[j]) }

# Some reporting to ensure all worked OK:
	print("Block size frequencies:")
	print(table(table(alloc$block)))
	print("Sample sizes:")
	print(table(alloc$arm))

	# Output the data to comma separated value file:
	write.csv(alloc, "allocation_list.csv", row.names=F)
	