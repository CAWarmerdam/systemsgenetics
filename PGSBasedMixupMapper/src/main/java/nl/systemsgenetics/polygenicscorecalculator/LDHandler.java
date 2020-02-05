package nl.systemsgenetics.polygenicscorecalculator;

import nl.systemsgenetics.gwassummarystatistics.EffectAllele;

import java.util.Iterator;
import java.util.List;

public interface LDHandler {
    public Iterator<EffectAllele> effectAlleleIterator(Iterator<EffectAllele> effectAlleleIterator);

    public double getRSquaredThreshold();

    public List<Integer> getWindowSize();
}
