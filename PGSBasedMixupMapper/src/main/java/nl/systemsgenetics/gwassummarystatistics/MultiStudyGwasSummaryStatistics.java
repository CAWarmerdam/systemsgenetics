package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.Iterator;

public interface MultiStudyGwasSummaryStatistics extends RandomAccessGenotypeData {
    float[][] getEffectSizeEstimates(GeneticVariant variant);

    float[][] getTransformedPValues(GeneticVariant variant);

    Iterator<EffectAllele> effectAlleles(String studyName);

    float[][] getStandardErrorOfES(GeneticVariant variant);
}
