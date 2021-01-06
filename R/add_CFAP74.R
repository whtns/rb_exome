##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
add_CFAP74 <- function() {

    replace_c1orf <- tibble::tribble(
        ~Uploaded.variant,           ~Location, ~Allele,       ~Consequence,  ~Symbol,             ~Gene, ~Feature.type,            ~Feature,         ~Biotype,  ~Exon, ~cDNA.position, ~CDS.position, ~Protein.position, ~Amino.acids,   ~Codons, ~Existing.variant, ~Feature.strand, ~Transcript.support.level, ~APPRIS, ~SIFT, ~PolyPhen,
        "1_1985486_C/T", "1:1985486-1985486",    TRUE, "missense_variant", "CFAP74", "ENSG00000142609",  "Transcript", "ENST00000493964.5", "protein_coding", "6/38",           556L,          400L,              134L,        "A/T", "GCT/ACT",    "rs1446994577",             -1L,                        5L,    "P1",  0.05,     0.069
    )

}
