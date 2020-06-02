package nl.systemsgenetics.gwassummarystatistics.effectAllele;

import nl.systemsgenetics.gwassummarystatistics.GeneticVariantBackedGwasSummaryStatistics;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Represents an effect allele, possibly originating from a genome wide association study.
 * This effect allele is backed by a genetic variant. This variant may or may not be able
 * to provide sample genotypes.
 * @author Robert Warmerdam
 */
public class GeneticVariantBackedEffectAllele extends EffectAllele {
    private final int alleleIndex;
    private final GeneticVariant variant;
    private final GeneticVariantBackedGwasSummaryStatistics summaryStatistics;

    public GeneticVariantBackedEffectAllele(GeneticVariant variant,
                                            GeneticVariantBackedGwasSummaryStatistics summaryStatistics,
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
        return summaryStatistics.getEffectSizeEstimate(this);
    }

    @Override
    public double getLogTransformedPValue() {
        return summaryStatistics.getLogTransformedPValue(this);
    }

    @Override
    public double getAlleleFrequency() {
        return summaryStatistics.getAlleleFrequency(this);
    }

    @Override
    public boolean variantIsSnp() {
        return variant.isSnp();
    }

    @Override
    public boolean variantIsBiallelic() {
        return variant.isBiallelic();
    }

    @Override
    public double getPValue() {
        return summaryStatistics.getPValue(this);
    }

    @Override
    public Allele getAllele() {
        return variant.getAlternativeAlleles().get(alleleIndex);
    }

    @Override
    public Allele getNonEffectAllele() {
        return variant.getRefAllele();
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

    /**
     * Returns whether or not the given variant corresponds to this effect allele.
     *
     * @param variant The genetic variant to compare with.
     * @return true if the variant matches the primary variant id,
     * sequence name and the starting position. Additionally, the given variant should
     * have the same alleles as the variant backing this effect allele, or their complement alleles,
     * for this method to return true.
     */
    @Override
    public boolean matchesVariant(GeneticVariant variant) {
        return super.matchesVariant(variant)
                && matchesVariantAllelesOrComplement(variant);
    }

    /**
     * Method that assesses whether or not the alleles, or their complements, of the specified genetic variant match
     * the alleles of the variant backing this effect allele.
     *
     * @param variant The variant to assess the alleles from.
     * @return true if the alleles of the specified variant, or their complement alleles,
     * match the alleles of the genetic variant backing this effect allele.
     */
    private boolean matchesVariantAllelesOrComplement(GeneticVariant variant) {
        if (variant.getVariantAlleles().sameAlleles(this.variant.getVariantAlleles())) {
            return true;
        }
        Alleles complement = variant.getVariantAlleles().getComplement();

        return (complement.sameAlleles(this.variant.getVariantAlleles()));
    }

    public int getAlleleIndex() {
        return alleleIndex;
    }
}
