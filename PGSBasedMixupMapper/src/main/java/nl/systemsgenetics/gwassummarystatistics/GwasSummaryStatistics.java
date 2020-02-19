package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.molgenis.genotype.RandomAccessGenotypeData;

import java.util.ArrayList;

public interface GwasSummaryStatistics extends Iterable<EffectAllele> {
    public String getGwasId();

    THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskEntries(
            RandomAccessGenotypeData genotypeData,
            double[] pValThres,
            String[] genomicRangesToExclude,
            boolean unweighted);
}
