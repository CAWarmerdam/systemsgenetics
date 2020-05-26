package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.ArrayList;

public interface GwasSummaryStatistics extends Iterable<EffectAllele> {
    public String getGwasId();

    /**
     * Get all effect alleles within the specified range
     *
     * @param seqName
     * @param rangeStart start of range, inclusive
     * @param rangeEnd end of range exclusive
     * @return
     */
    Iterable<EffectAllele> getEffectAllelesByRange(String seqName, int rangeStart, int rangeEnd);


}
