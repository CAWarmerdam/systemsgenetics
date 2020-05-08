package nl.systemsgenetics.gwassummarystatistics.effectAlleleFilter;

import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import org.molgenis.genotype.variant.GeneticVariant;

public interface EffectAlleleFilter {

    public boolean doesEffectAllelePassFilter(EffectAllele variant);

    /**
     * Most implementation always return true. Only implementation that only
     * filter on primary variant ID should return false if not pass
     *
     * @param id
     * @return
     */
    public boolean doesVariantIdPassFilter(String id);
}
