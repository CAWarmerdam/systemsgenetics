package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.ArrayList;
import java.util.Iterator;

public interface GwasSummaryStatistics {
    public String getGwasId();

    public float[] getEffectSizeEstimates(GeneticVariant variant);

    public float[] getTransformedPValues(GeneticVariant variant);

    public Iterator<EffectAllele> effectAlleles();

    THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskEntries(
            RandomAccessGenotypeData genotypeData,
            double[] pValThres,
            String[] genomicRangesToExclude,
            boolean unweighted);
}
