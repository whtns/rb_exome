##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vc_somatic_variants
##' @return
##' @author whtns
##' @export
plot_vaf_density <- function(vc_somatic_variants) {
    
    vaf_dist_plot <- ggplot(vc_somatic_variants, aes(x = AF.TUMOR, color = caller)) + 
        geom_density() + 
        theme_cowplot() +
        labs(title = "Variant allele frequency density by variant caller",
             x = "Variant allele frequency",
             y = "Density",
             color = "Variant Caller")
    
    # ------------------------------
    
    # fit <- euler(vc_somatic_variants)
    # 
    # venn_caller_overlap <- plot(fit)
    # 
    # out_plot <- vaf_dist_plot / 
    #     venn_caller_overlap + 
    #     plot_layout(heights = c(2, 1)) + 
    #     plot_annotation(tag_levels = 'A')
    
    # ggsave("~/rb_pipeline/doc/RB_exome_manuscript/stachelek_supplemental/fig_s08.pdf", out_plot)
    

}
