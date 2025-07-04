#!/usr/bin/env Rscript

# Comprehensive SoupX Visualization Examples
# This script demonstrates all the new visualization and analysis capabilities

library(SoupX)
library(ggplot2)

cat("=== SoupX Enhanced Visualization Examples ===\n\n")

# Load test data
data(scToy)
cat("✓ Loaded scToy test data\n")

# Run basic SoupX analysis
cat("Running SoupX analysis...\n")
sc_est <- estimateSoup(scToy)
sc_cont <- autoEstCont(sc_est, verbose=FALSE, doPlot=FALSE)
adjusted <- adjustCounts(sc_cont, verbose=0)
cat("✓ Completed SoupX analysis\n\n")

# Example 1: Quality Control Dashboard
cat("1. Quality Control Dashboard\n")
cat("============================\n")
qc_plots <- plotQualityControl(sc_cont, adjusted)
cat("✓ Generated QC plots:\n")
cat("  - Top soup genes\n")
cat("  - Contamination distribution\n")
cat("  - UMI count distribution\n")
cat("  - Before/after comparison\n")
cat("  - Soup correlation\n\n")

# Display individual plots
print(qc_plots$soup_profile)
print(qc_plots$contamination_distribution)
print(qc_plots$umi_distribution)
if("before_after_comparison" %in% names(qc_plots)) {
  print(qc_plots$before_after_comparison)
}

# Example 2: Cluster-specific Contamination Analysis
cat("2. Cluster-specific Contamination Analysis\n")
cat("==========================================\n")
if("clusters" %in% colnames(sc_cont$metaData)) {
  cluster_plot <- plotContaminationByCluster(sc_cont, plotType = "boxplot")
  print(cluster_plot)
  
  # Try different plot types
  cluster_violin <- plotContaminationByCluster(sc_cont, plotType = "violin")
  cluster_bar <- plotContaminationByCluster(sc_cont, plotType = "bar")
  
  cat("✓ Generated cluster contamination plots (boxplot, violin, bar)\n\n")
} else {
  cat("⚠ No clustering information available for cluster analysis\n\n")
}

# Example 3: Gene-specific Contamination Analysis
cat("3. Gene-specific Contamination Analysis\n")
cat("=======================================\n")

# Get some example genes from the data
example_genes <- c("CD7", "LTB", "S100A9", "LYZ", "CD74")
available_genes <- intersect(example_genes, rownames(sc_cont$toc))

if(length(available_genes) > 0) {
  # Contamination ratio analysis
  gene_contamination <- plotGeneContamination(sc_cont, available_genes, 
                                             plotType = "contamination_ratio")
  print(gene_contamination)
  
  # Expression change analysis
  gene_change <- plotGeneContamination(sc_cont, available_genes, adjusted, 
                                      plotType = "expression_change")
  print(gene_change)
  
  # Soup contribution analysis
  gene_soup <- plotGeneContamination(sc_cont, available_genes, 
                                    plotType = "soup_contribution")
  print(gene_soup)
  
  cat("✓ Generated gene-specific contamination plots\n\n")
} else {
  cat("⚠ No example genes found in the data\n\n")
}

# Example 4: Comprehensive Report Generation
cat("4. Comprehensive Report Generation\n")
cat("==================================\n")
report <- generateSoupXReport(sc_cont, adjusted)
cat("✓ Generated comprehensive report with:\n")
cat("  - All QC plots\n")
cat("  - Cluster analysis (if available)\n")
cat("  - Soup correlation\n")
cat("  - Marker distribution (if clusters available)\n")
cat("  - Summary statistics\n\n")

# Display summary statistics
cat("Summary Statistics:\n")
cat("==================\n")
cat("Number of cells:", report$summary$n_cells, "\n")
cat("Number of genes:", report$summary$n_genes, "\n")
cat("Total UMIs:", report$summary$total_umis, "\n")
cat("Sparsity:", round(report$summary$sparsity * 100, 2), "%\n")
if("mean_contamination" %in% names(report$summary)) {
  cat("Mean contamination:", round(report$summary$mean_contamination, 4), "\n")
  cat("Contamination range:", round(report$summary$contamination_range, 4), "\n")
}
if("reduction_percent" %in% names(report$summary)) {
  cat("UMIs removed:", report$summary$umis_removed, "\n")
  cat("Reduction percentage:", round(report$summary$reduction_percent, 2), "%\n")
}
cat("\n")

# Example 5: Advanced Visualizations
cat("5. Advanced Visualizations\n")
cat("==========================\n")

# Create a custom visualization showing the relationship between UMI count and contamination
if("rho" %in% colnames(sc_cont$metaData)) {
  umi_vs_contamination <- ggplot(data.frame(
    UMIs = sc_cont$metaData$nUMIs,
    Contamination = sc_cont$metaData$rho
  ), aes(x = log10(UMIs), y = Contamination)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "loess", color = "red") +
    labs(title = "UMI Count vs Contamination Fraction",
         x = "log10(UMIs per Cell)",
         y = "Contamination Fraction") +
    theme_minimal()
  
  print(umi_vs_contamination)
  cat("✓ Generated UMI vs contamination correlation plot\n\n")
}

# Example 6: Performance Benchmarking
cat("6. Performance Benchmarking\n")
cat("===========================\n")
benchmark_results <- benchmark_soupx(sc_cont, iterations = 2, verbose = TRUE)
cat("✓ Completed performance benchmarking\n\n")

# Example 7: Interactive Analysis Suggestions
cat("7. Interactive Analysis Suggestions\n")
cat("==================================\n")
cat("For interactive analysis, consider:\n")
cat("- Using plotly::ggplotly() to make plots interactive\n")
cat("- Creating Shiny apps for dynamic parameter exploration\n")
cat("- Using DT::datatable() for interactive data tables\n")
cat("- Implementing custom color schemes for different cell types\n\n")

# Example 8: Export Options
cat("8. Export Options\n")
cat("=================\n")
cat("To save all plots to PDF:\n")
cat("report <- generateSoupXReport(sc_cont, adjusted, output_dir = './soupx_report')\n\n")

cat("To save individual plots:\n")
cat("ggsave('soup_profile.pdf', qc_plots$soup_profile, width = 10, height = 8)\n")
cat("ggsave('contamination_by_cluster.pdf', cluster_plot, width = 10, height = 8)\n\n")

cat("=== Visualization Examples Complete ===\n")
cat("✓ All new visualization functions demonstrated\n")
cat("✓ Quality control dashboard created\n")
cat("✓ Cluster and gene-specific analyses shown\n")
cat("✓ Comprehensive report generated\n")
cat("✓ Performance benchmarking completed\n\n")

cat("Key Benefits of New Visualizations:\n")
cat("===================================\n")
cat("1. Quality Control: Comprehensive QC dashboard for data assessment\n")
cat("2. Cluster Analysis: Identify cluster-specific contamination patterns\n")
cat("3. Gene Analysis: Understand which genes are most affected by contamination\n")
cat("4. Before/After Comparison: Visualize the impact of correction\n")
cat("5. Automated Reporting: Generate publication-ready reports\n")
cat("6. Performance Monitoring: Track computational performance\n")
cat("7. Interactive Options: Ready for interactive exploration\n\n") 