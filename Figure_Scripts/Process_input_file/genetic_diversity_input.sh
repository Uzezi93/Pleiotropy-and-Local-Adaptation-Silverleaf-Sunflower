#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

# Measuring sfs for each population

module load angsd/0.935

# angsd -vcf-gl north_lfmm.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/north_lfmm.sfs

# angsd -vcf-gl coast_lfmm.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/coast_lfmm.sfs

# angsd -vcf-gl north_pcadapt.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/north_pcadapt.sfs

# angsd -vcf-gl coast_pcadapt.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/coast_pcadapt.sfs

angsd -vcf-gl coast_eQTL2.recode.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS3/coast_eQTL.sfs

angsd -vcf-gl north_eQTL2.recode.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS3/north_eQTL.sfs

angsd -vcf-gl coast_eGenes2.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS3/coast_eGenes.sfs

angsd -vcf-gl north_eGenes2.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS3/north_eGenes.sfs

angsd -vcf-gl coast_control2.recode.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS3/coast_control.sfs

angsd -vcf-gl north_control2.recode.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS3/north_control.sfs
#---------------------------------------------------------------------------------------------------------

# Measuring folded saf for each population

 for i in SFS3/*.idx
 do
  realSFS $i -fold 1 > $i.sfs
  done

#----------------------------------------------------------------------------------------------------------
# Measuring unfolded saf for each population

 for i in SFS3/*.idx
 do
 realSFS $i > ${i}_unfold.sfs
 done

#----------------------------------------------------------------------------------------------------------
# Calculate theta for each site

 for i in SFS3/*.idx
 do
 realSFS saf2theta $i -sfs ${i}_unfold.sfs -outname $i
 done

#----------------------------------------------------------------------------------------------------------
# Print theta statistics

 for i in SFS3/*.thetas.idx
 do
 thetaStat print $i > ${i}_theta_stat
 done

#----------------------------------------------------------------------------------------------------------
# Estimate Tajimas D and other statistics

 for i in SFS3/*.thetas.idx
 do
 thetaStat do_stat $i
 done

#---------------------------------------------------------------------------------------------------------
# Estimating genetic diversity for coastal populations shared loci between PCAdapt and LFMM
# Make seperate directories for seperate analysis

angsd -vcf-gl coast_shared.recode.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/coast_shared.sfs

realSFS SFS2/coast_shared.sfs.saf.idx -fold 1 > SFS2/coast_shared.sfs

realSFS SFS2/coast_shared.sfs.saf.idx  > SFS2/coast_shared_unfold.sfs

realSFS saf2theta SFS2/all_pop.sfs.saf.idx -sfs SFS2/all_pop_unfold.sfs -outname SFS2/coast_shared

cd SFS2/

thetaStat print coast_shared.thetas.idx  > coast_shared_theta_stat

thetaStat do_stat coast_shared.thetas.idx

#---------------------------------------------------------------------------------------------------------
# Estimating genetic diversity for northern populations shared loci between PCAdapt and LFMM

cd ../

angsd -vcf-gl north_shared.recode.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/north_shared.sfs

realSFS SFS2/north_shared.sfs.saf.idx -fold 1 > SFS2/north_shared.sfs

realSFS SFS2/north_shared.sfs.saf.idx  > SFS2/north_shared_unfold.sfs

realSFS saf2theta SFS2/all_pop.sfs.saf.idx -sfs SFS2/all_pop_unfold.sfs -outname SFS2/north_shared

cd SFS2/

thetaStat print north_shared.thetas.idx  > north_shared_theta_stat

thetaStat do_stat north_shared.thetas.idx

#-----------------------------------------------------------------------------------------------------------
# Estimating genetic diversity for northern populations

angsd -vcf-gl north_pop.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/north_pop.sfs

realSFS SFS2/north_pop.sfs.saf.idx -fold 1 > SFS2/north_pop.sfs

realSFS SFS2/north_pop.sfs.saf.idx  > SFS2/north_pop_unfold.sfs

realSFS saf2theta SFS2/north_pop.sfs.saf.idx -sfs SFS2/north_pop_unfold.sfs -outname SFS2/north_pop

cd SFS2/

thetaStat print north_pop.thetas.idx  > north_pop_theta_stat

thetaStat do_stat north_pop.thetas.idx

#-----------------------------------------------------------------------------------------------------------
# Estimating genetic diversity for coastal populations

cd ../

angsd -vcf-gl coast_pop.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/coast_pop.sfs

realSFS SFS2/coast_pop.sfs.saf.idx -fold 1 > SFS2/coast_pop.sfs

realSFS SFS2/coast_pop.sfs.saf.idx  > SFS2/coast_pop_unfold.sfs

realSFS saf2theta SFS2/coast_pop.sfs.saf.idx -sfs SFS2/coast_pop_unfold.sfs -outname SFS2/coast_pop

cd SFS2/

thetaStat print coast_pop.thetas.idx  > coast_pop_theta_stat

thetaStat do_stat coast_pop.thetas.idx
