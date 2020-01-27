package nl.systemsgenetics.polygenicscorecalculator;

import nl.systemsgenetics.gwassummarystatistics.EffectAllele;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.*;

public class Clumper implements LDHandler {

    private Map<Set<GeneticVariant>, Double> ldMatrix;
    private int windowSize;
    private double rSquaredThreshold;
    private double indexVariantPValueThreshold;
    private List<GeneticVariant> indexVariants;

    public Clumper(RandomAccessGenotypeData genotypeData,
                   int windowSize, double rSquared, double indexVariantPValueThreshold) throws LdCalculatorException {

        this(calculateLDMatrix(genotypeData, windowSize, rSquared), windowSize, rSquared, indexVariantPValueThreshold);
    }

    public Clumper(Map<Set<GeneticVariant>, Double> ldMatrix, int windowSize, double rSquared, double indexVariantPValueThreshold) {
        this.ldMatrix = ldMatrix;
        this.windowSize = windowSize;
        this.rSquaredThreshold = rSquared;
        this.indexVariantPValueThreshold = indexVariantPValueThreshold;
    }

    public Iterator<EffectAllele> effectAlleleIterator(Iterator<EffectAllele> riskEntries) {
        indexVariants = new LinkedList<>();
        return new Iterator<EffectAllele>() {
            private EffectAllele nextEffectAllele = null;
            private EffectAllele lastReturnedEffectAllele = null;
            private boolean nextCalled = false;

            private boolean setNextRiskEntry() {
                // While there are more risk entries, loop through them to check if
                // one could be the next risk entry
                while (!riskEntries.hasNext()) {
                    nextEffectAllele = riskEntries.next();
                    if (nextEffectAllele == null) {
                        return false;
                    } else if (canFormChunk(nextEffectAllele.getVariant())) {
                        // If the a chunk can be formed, return true.
                        indexVariants.add(nextEffectAllele.getVariant());
                        return true;
                    }
                }
                // Return false and set the next risk entry to null if all risk entries have been checked.
                nextEffectAllele = null;
                return false;
            }

            @Override
            public boolean hasNext() {
                // Check if there if the next riskEntry exists and if it can form its own chunk
                if (!nextCalled && !riskEntries.hasNext()) {
                    return true;
                }
                return setNextRiskEntry();
            }

            @Override
            public EffectAllele next() {
                if (lastReturnedEffectAllele.equals(nextEffectAllele)) {
                    if (!hasNext()) {
                        throw new NoSuchElementException();
                    }
                }
                if (nextEffectAllele == null) {
                    throw new NoSuchElementException();
                }
                nextCalled = true;
                lastReturnedEffectAllele = nextEffectAllele;
                return nextEffectAllele;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private boolean canFormChunk(GeneticVariant newVariant) {

        for (GeneticVariant indexVariant :
                indexVariants) {
            // Check if the LD should be calculated on the spot
            // If the LD should be calculated on the spot, to this.
            Set<GeneticVariant> geneticVariantPair = Collections.unmodifiableSet(
                    new HashSet<>(Arrays.asList(indexVariant, newVariant)));
            if (ldMatrix.containsKey(geneticVariantPair)) {
                return false;
            }
        }
        return true;
    }

    private static Map<Set<GeneticVariant>, Double> calculateLDMatrix(
            RandomAccessGenotypeData genotypeData, int windowSize, double rSquaredThreshold) throws LdCalculatorException {
        HashMap<Set<GeneticVariant>, Double> ldMatrix = new HashMap<>();
        return ldMatrix;

//        // Get the chromosomes / sequences in the genotype data
//        for (String sequence : genotypeData.getSeqNames()) {
//
//            // Loop through the variants in this sequence in order (low bp to high bp)
//            Iterable<GeneticVariant> sequenceVariants = genotypeData.getSequenceGeneticVariants(sequence);
//            // TODO: Make sure that these iterables go through the variants ordered by bp position
//
//            for (GeneticVariant variant : sequenceVariants) {
//                int startPos = variant.getStartPos();
//                // Now loop through the variants starting from the current variant stopping at the variant
//                // for which the starting position difference with the current startpos does not exceed the window size.
//                for (GeneticVariant otherVariant : genotypeData.getVariantsByRange(
//                        sequence, startPos + 1, startPos + windowSize + 1)) {
//                    if (!variant.equals(otherVariant)) {
//                        // Calculate the R2
//                        double rSquared = LdCalculator.calculateRsquare(variant, otherVariant);
//                        // If the R2 is above the required threshold, put this pair in the map / matrix
//                        if (rSquared >= rSquaredThreshold) {
//                            ldMatrix.put( // Create a set as a key because order is not important this way
//                                    Collections.unmodifiableSet(new HashSet<>(Arrays.asList(variant, otherVariant))),
//                                    rSquared);
//                        }
//                    }
//                }
//            }
//        }
//        return ldMatrix;
    }

    public Map<Set<GeneticVariant>, Double> getLdMatrix() {
        return ldMatrix;
    }

    public double getIndexVariantPValueThreshold() {
        return indexVariantPValueThreshold;
    }

    public List<GeneticVariant> getIndexVariants() {
        return indexVariants;
    }

    @Override
    public int getWindowSize() {
        return windowSize;
    }

    @Override
    public double getRSquaredThreshold() {
        return rSquaredThreshold;
    }
}
