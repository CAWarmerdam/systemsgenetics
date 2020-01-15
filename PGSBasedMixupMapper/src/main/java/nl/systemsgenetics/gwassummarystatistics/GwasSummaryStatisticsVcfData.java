package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.vcf.VcfGenotypeData;

public class GwasSummaryStatisticsVcfData {
    public GwasSummaryStatisticsVcfData(VcfGenotypeData vcfGwasSummaryStatistics) {
        // Check the following
        // Variant IDs must be unique!
        // Reserved keys:
        // Field 	Description 	Required
        // ES 	Effect size estimate relative to the alternative allele 	YES!
        // SE 	Standard error of effect size estimate 	YES!
        // LP 	-log10 p-value for effect estimate 	NO
        // AF 	Alternate allele frequency in the association study 	NO
        // SS 	Sample size used to estimate genetic effect 	NO
        // EZ 	Z-score provided if it was used to derive the EFFECT and SE fields 	NO
        // SI 	Accuracy score of summary data imputation 	NO
        // NC 	Number of cases used to estimate genetic effect 	NO
        // ID 	Study variant identifier 	NO

    }
}
