package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.GeneticVariantBackedEffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import org.molgenis.genotype.RandomAccessGenotypeData;

import java.util.ArrayList;
import java.util.Iterator;

public interface GeneticVariantBackedGwasSummaryStatistics extends GwasSummaryStatistics {

    double getEffectSizeEstimate(GeneticVariantBackedEffectAllele effectAllele);

    double getLogTransformedPValue(GeneticVariantBackedEffectAllele effectAllele);

    double getPValue(GeneticVariantBackedEffectAllele effectAllele);

    double getAlleleFrequency(GeneticVariantBackedEffectAllele effectAllele);
}
