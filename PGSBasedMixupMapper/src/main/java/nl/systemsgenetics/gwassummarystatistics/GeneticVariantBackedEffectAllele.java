package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;

public class GeneticVariantBackedEffectAllele extends EffectAllele {
    private final int alleleIndex;
    private final GeneticVariant variant;
    private final ReadOnlyGwasSummaryStatistics summaryStatistics;

    public GeneticVariantBackedEffectAllele(GeneticVariant variant, ReadOnlyGwasSummaryStatistics summaryStatistics,
                                            int alleleIndex) {
        this.variant = variant;
        this.alleleIndex = alleleIndex;
        this.summaryStatistics = summaryStatistics;
    }

    public GeneticVariant getVariant() {
        return variant;
    }

    @Override
    public double getEffectSize() {
        return summaryStatistics.getEffectSizeEstimate(variant, alleleIndex);
    }

    @Override
    public double getLogTransformedPValue() {
        return summaryStatistics.getLogTransformedPValue(variant, alleleIndex);
    }

    @Override
    public double getPValue() {
        return summaryStatistics.getPValue(variant, alleleIndex);
    }

    @Override
    public Allele getAllele() {
        return variant.getAlternativeAlleles().get(alleleIndex);
    }

    @Override
    public String getSequenceName() {
        return variant.getSequenceName();
    }

    @Override
    public int getStartPos() {
        return variant.getStartPos();
    }

    @Override
    public String getPrimaryVariantId() {
        return variant.getPrimaryVariantId();
    }

    @Override
    public boolean matchesVariant(GeneticVariant variant) {
        return super.matchesVariant(variant)
                && matchesVariantAllelesOrComplement(variant);
    }

    public boolean matchesVariantAllelesOrComplement(GeneticVariant variant) {
        if (variant.getVariantAlleles().sameAlleles(this.variant.getVariantAlleles())) {
            return true;
        }
        Alleles complement = variant.getVariantAlleles().getComplement();

        return (complement.sameAlleles(this.variant.getVariantAlleles()));
    }
}
