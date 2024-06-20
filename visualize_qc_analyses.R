args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}







make_h2_comparison_scatter <- function(qc_h2_comparison_file) {
	df <- read.table(qc_h2_comparison_file, header=TRUE, sep="\t")


	pp <- ggplot(df, aes(x=h2_2021, y=h2_2024)) + geom_point(color='magenta3') +
		figure_theme() + 
		labs(x="h2 (2021)", y="h2 (2024)", title="h2") +
		geom_abline(slope=1, color='grey')

	return(pp)

}


make_h2_se_comparison_scatter <- function(qc_h2_comparison_file) {
	df <- read.table(qc_h2_comparison_file, header=TRUE, sep="\t")


	pp <- ggplot(df, aes(x=h2_se_2021, y=h2_se_2024)) + geom_point(color='magenta3') +
		figure_theme() + 
		labs(x="h2_se (2021)", y="h2_se (2024)", title="h2_se") +
		geom_abline(slope=1, color='grey')

	return(pp)

}




make_h2_z_comparison_scatter <- function(qc_h2_comparison_file) {
	df <- read.table(qc_h2_comparison_file, header=TRUE, sep="\t")


	pp <- ggplot(df, aes(x=h2_z_2021, y=h2_z_2024)) + geom_point(color='magenta3') +
		figure_theme() + 
		labs(x="h2_z (2021)", y="h2_z (2024)", title="h2_z") +
		geom_abline(slope=1, color='grey')

	return(pp)

}



make_h2_intercept_comparison_scatter <- function(qc_h2_comparison_file) {
	df <- read.table(qc_h2_comparison_file, header=TRUE, sep="\t")


	pp <- ggplot(df, aes(x=intercept_2021, y=intercept_2024)) + geom_point(color='magenta3') +
		figure_theme() + 
		labs(x="intercept (2021)", y="intercept (2024)", title="Intercept") +
		geom_abline(slope=1, color='grey')

	return(pp)

}


make_rg_comparison_scatter <- function(qc_rg_comparison_file) {
	df <- read.table(qc_rg_comparison_file, header=TRUE, sep="\t")

	pp <- ggplot(df, aes(x=rg_2021, y=rg_2024)) + geom_point(color='cyan4') +
		figure_theme() + 
		labs(x="rg (2021)", y="rg (2024)", title="rg") +
		geom_abline(slope=1, color='grey')

	return(pp)

}



make_h2_distribution_histogram <- function(alkesgrp_2024_EUR_h2_summary_file) {
	df <- read.table(alkesgrp_2024_EUR_h2_summary_file, header=TRUE, sep="\t")

	pp <- ggplot(df, aes(x=h2))+
 		 geom_histogram(color="magenta4", fill="magenta1", alpha=.14) +
 		 figure_theme() + 
 		 labs(x="SLDSC h2")

 	return(pp)
}


make_rg_distribution_histogram <- function(alkesgrp_2024_EUR_genetic_corr_summary_file) {
	df <- read.table(alkesgrp_2024_EUR_genetic_corr_summary_file, header=TRUE, sep="\t")

	pp <- ggplot(df, aes(x=genetic_correlation))+
 		 geom_histogram(color="cyan4", fill="cyan1", alpha=.14) +
 		 figure_theme() + 
 		 labs(x="cross-trait LDSC rg")

 	return(pp)
}



############################
# Command line args
############################
alkesgrp_2024_EUR_h2_summary_file = args[1]
alkesgrp_2024_EAS_h2_summary_file = args[2]
alkesgrp_2024_EUR_genetic_corr_summary_file = args[3]
alkesgrp_2024_EAS_genetic_corr_summary_file = args[4]
qc_comparison_to_2021_update_dir = args[5]


# Qc comparison summary files
qc_h2_comparison_file = paste0(qc_comparison_to_2021_update_dir, "h2_comparison.txt")
qc_rg_comparison_file = paste0(qc_comparison_to_2021_update_dir, "rg_comparison.txt")



################
# Scatter plot comparing h2 in 2024 and 2021
output_file <- paste0(qc_comparison_to_2021_update_dir, "h2_comparison_2024_2021.pdf")

h2_scatter <- make_h2_comparison_scatter(qc_h2_comparison_file)
h2_se_scatter <- make_h2_se_comparison_scatter(qc_h2_comparison_file)
h2_z_scatter <- make_h2_z_comparison_scatter(qc_h2_comparison_file)
h2_intercept_scatter <- make_h2_intercept_comparison_scatter(qc_h2_comparison_file)

joint_plot <- plot_grid(h2_scatter, h2_se_scatter, h2_z_scatter, h2_intercept_scatter, ncol=2)
ggsave(joint_plot, file=output_file, width=7.2, height=5.3, units="in")


################
# Scatter plot comparing rg in 2024 and 2021
output_file <- paste0(qc_comparison_to_2021_update_dir, "rg_comparison_2024_2021.pdf")
rg_scatter <- make_rg_comparison_scatter(qc_rg_comparison_file)
ggsave(rg_scatter, file=output_file, width=7.2, height=3.9, units="in")

################
# histogram plot showing distribution of heritabilities
output_file <- paste0(qc_comparison_to_2021_update_dir, "sldsc_h2_eur_distribution_histogram.pdf")
h2_histo <- make_h2_distribution_histogram(alkesgrp_2024_EUR_h2_summary_file)
ggsave(h2_histo, file=output_file, width=7.2, height=3.9, units="in")


################
# histogram plot showing distribution of genetic correlations
output_file <- paste0(qc_comparison_to_2021_update_dir, "ldsc_rg_eur_distribution_histogram.pdf")
rg_histo <- make_rg_distribution_histogram(alkesgrp_2024_EUR_genetic_corr_summary_file)
ggsave(rg_histo, file=output_file, width=7.2, height=3.9, units="in")







