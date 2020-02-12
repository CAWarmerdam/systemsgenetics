package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;

import java.util.Iterator;

public class VariantFilterableGwasSummaryStatisticsDecorator extends VariantFilterableGenotypeDataDecorator implements MultiStudyGwasSummaryStatistics {
    public VariantFilterableGwasSummaryStatisticsDecorator(MultiStudyGwasSummaryStatistics originalGenotypeData, VariantFilter variantFilter) {
        super(originalGenotypeData, variantFilter);
    }

    @Override
    public float[][] getEffectSizeEstimates(GeneticVariant variant) {
        return ((MultiStudyGwasSummaryStatistics) originalGenotypeData).getEffectSizeEstimates(variant);
    }

    @Override
    public float[][] getTransformedPValues(GeneticVariant variant) {
        return ((MultiStudyGwasSummaryStatistics) originalGenotypeData).getTransformedPValues(variant);
    }

    @Override
    public Iterator<EffectAllele> effectAlleles(ReadOnlyGwasSummaryStatistics studyIndex) {
        return ((MultiStudyGwasSummaryStatistics) originalGenotypeData).effectAlleles(studyIndex);
    }

    @Override
    public float[][] getStandardErrorOfES(GeneticVariant variant) {
        return ((MultiStudyGwasSummaryStatistics) originalGenotypeData).getStandardErrorOfES(variant);
    }
}
