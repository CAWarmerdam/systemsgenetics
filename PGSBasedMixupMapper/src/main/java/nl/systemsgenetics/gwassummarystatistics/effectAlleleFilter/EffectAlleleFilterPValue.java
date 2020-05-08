package nl.systemsgenetics.gwassummarystatistics.effectAlleleFilter;

import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;

public class EffectAlleleFilterPValue implements EffectAlleleFilter {

    private double minimumPValue;
    private double maximumPValue;

    public EffectAlleleFilterPValue(double minimumPValue, double maximumPValue) {
        this.minimumPValue = minimumPValue;
        this.maximumPValue = maximumPValue;
    }

    @Override
    public boolean doesEffectAllelePassFilter(EffectAllele variant) {
        double pValue = variant.getPValue();
        return pValue < maximumPValue && pValue > minimumPValue;
    }

    @Override
    public boolean doesVariantIdPassFilter(String id) {
        return true;
    }
}
