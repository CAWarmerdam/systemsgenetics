package nl.systemsgenetics.gwassummarystatistics;

import com.google.common.collect.Iterators;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GeneticVariantMetaMap;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.vcf.VcfGeneticVariantMeta;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.molgenis.vcf.meta.VcfMeta;
import umcg.genetica.io.trityper.converters.VCFToTriTyper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.util.stream.StreamSupport;

public class VcfGwasSummaryStatistics extends VcfGenotypeData{
    public VcfGwasSummaryStatistics(File bzipVcfFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        super(bzipVcfFile, cacheSize, minimumPosteriorProbabilityToCall);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        super(bzipVcfFile, tabixIndexFile, minimumPosteriorProbabilityToCall);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        super(bzipVcfFile, tabixIndexFile, cacheSize, minimumPosteriorProbabilityToCall);
    }

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

        // Check if variant IDs are unique
//        long numberOfUniqueVariantIds = StreamSupport.stream(
//                Spliterators.spliteratorUnknownSize(vcfGwasSummaryStatistics.iterator(), Spliterator.ORDERED),
//                false).map(GeneticVariant::getPrimaryVariantId).distinct().count();
//        int numberOfVariants = Iterators.size(vcfGwasSummaryStatistics.iterator());
//        if (numberOfUniqueVariantIds != numberOfVariants) {
//            throw new GwasSummaryStatisticsVcfDataException(String.format("Encountered an error in the VCF file, " +
//                    "VCF file does not contain exclusively unique variant identifiers. (%d vs %d)",
//                    numberOfUniqueVariantIds, numberOfVariants));
//        }

    public double getStandardErrorOfEffectSizeEstimate(GeneticVariant variant) {
        GeneticVariantMeta variantMeta = variant.getVariantMeta();
        if (variantMeta instanceof VcfGeneticVariantMeta) {
            VcfGeneticVariantMeta variantMeta1 = (VcfGeneticVariantMeta) variantMeta;
        }
        return 0;
    }
}
