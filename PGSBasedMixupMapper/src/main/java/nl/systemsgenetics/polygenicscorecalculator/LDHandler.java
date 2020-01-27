package nl.systemsgenetics.polygenicscorecalculator;

import nl.systemsgenetics.gwassummarystatistics.EffectAllele;

import java.util.Iterator;

public interface LDHandler {
    public Iterator<EffectAllele> effectAlleleIterator(Iterator<EffectAllele> effectAlleleIterator);

    public double getRSquaredThreshold();

    public int getWindowSize();
}
