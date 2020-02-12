package nl.systemsgenetics.polygenicscorecalculator;

import nl.systemsgenetics.gwassummarystatistics.EffectAllele;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;

import javax.management.RuntimeErrorException;
import java.util.*;

public class Clumper implements LDHandler {

    private Map<Set<ComparableGeneticVariant>, Double> ldMatrix;
    private List<Integer> windowSize;
    private double rSquaredThreshold;
    private double indexVariantPValueThreshold;
    private List<ComparableGeneticVariant> indexVariants;

    public Clumper(RandomAccessGenotypeData genotypeData,
                   List<Integer> windowSize,
                   double rSquared,
                   double indexVariantPValueThreshold) throws LdCalculatorException, LDHandlerException {
        this.ldMatrix = calculateLDMatrix(genotypeData, windowSize, rSquaredThreshold);
        this.windowSize = windowSize;
        this.rSquaredThreshold = rSquared;
        this.indexVariantPValueThreshold = indexVariantPValueThreshold;
    }

    public Clumper(Map<Set<ComparableGeneticVariant>, Double> ldMatrix,
                   List<Integer> windowSize,
                   double rSquared,
                   double indexVariantPValueThreshold) {
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

            private boolean setNextRiskEntry() {
                // While there are more risk entries, loop through them to check if
                // one could be the next risk entry
                while (riskEntries.hasNext()) {
                    nextEffectAllele = riskEntries.next();
                    if (nextEffectAllele == null) {
                        return false;
                    } else if (canFormChunk(nextEffectAllele)) {
                        // If the a chunk can be formed, return true.
                        indexVariants.add(new ComparableGeneticVariant(nextEffectAllele.getVariant()));
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
                if (nextEffectAllele != null && !lastReturnedEffectAllele.equals(nextEffectAllele)) {
                    return true;
                }
                return setNextRiskEntry();
            }

            @Override
            public EffectAllele next() {
                if (nextEffectAllele == null ||
                        lastReturnedEffectAllele == null ||
                        lastReturnedEffectAllele.equals(nextEffectAllele)) {
                    if (!setNextRiskEntry()) {
                        throw new NoSuchElementException();
                    }
                }
                lastReturnedEffectAllele = nextEffectAllele;
                return nextEffectAllele;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private boolean canFormChunk(EffectAllele effectAllele) {
        if (effectAllele.getLogTransformedPValue() < indexVariantPValueThreshold) {
            return false;
        }

        ComparableGeneticVariant newVariant = new ComparableGeneticVariant(effectAllele.getVariant());

        for (ComparableGeneticVariant indexVariant :
                indexVariants) {
            // Check if the LD should be calculated on the spot
            // If the LD should be calculated on the spot, to this.
            Set<ComparableGeneticVariant> geneticVariantPair = Collections.unmodifiableSet(
                    new HashSet<>(Arrays.asList(indexVariant, newVariant)));
            if (doesPairExceedRSquaredThreshold(geneticVariantPair) &&
                    areLocatedWithinOrEqualToWindowSize(newVariant, indexVariant)) {
                return false;
            }
        }
        return true;
    }

    private boolean doesPairExceedRSquaredThreshold(Set<ComparableGeneticVariant> geneticVariantPair) {
        return ldMatrix.containsKey(geneticVariantPair) && ldMatrix.get(geneticVariantPair) >= rSquaredThreshold;
//        } else {
//            ArrayList<ComparableGeneticVariant> comparableGeneticVariants = new ArrayList<>(geneticVariantPair);
//            GeneticVariant firstMatchingVariant = getMatchingVariant(comparableGeneticVariants.get(0));
//            GeneticVariant otherMatchingVariant = getMatchingVariant(comparableGeneticVariants.get(1));
//
//            try {
//                return firstMatchingVariant != null && otherMatchingVariant != null
//                        && LdCalculator.calculateRsquare(firstMatchingVariant, otherMatchingVariant) >= rSquaredThreshold;
//            } catch (LdCalculatorException e) {
//                throw new RuntimeException("Encountered unexpected exception when calculating R-squared", e);
//            }
//        }
    }

//    private GeneticVariant getMatchingVariant(ComparableGeneticVariant comparableGeneticVariant) {
//        GeneticVariant firstVariant = comparableGeneticVariant.originalVariant;
//        Iterable<GeneticVariant> candidateVariants = genotypeData.getVariantsByPos(
//                firstVariant.getSequenceName(), firstVariant.getStartPos());
//
//        GeneticVariant firstMatchingVariant = null;
//        for (GeneticVariant candidateVariant : candidateVariants) {
//            if (new ComparableGeneticVariant(candidateVariant).equals(comparableGeneticVariant)) {
//                firstMatchingVariant = candidateVariant;
//                break;
//            }
//        }
//        return firstMatchingVariant;
//    }

    private boolean areLocatedWithinOrEqualToWindowSize(ComparableGeneticVariant firstVariant,
                                                        ComparableGeneticVariant secondVariant) {
        int startPositionOfFirstVariant = secondVariant.getOriginalVariant().getStartPos();
        int startPositionOfSecondVariant = firstVariant.getOriginalVariant().getStartPos();
        return Math.abs(startPositionOfFirstVariant - startPositionOfSecondVariant) <= windowSize.get(0);
    }

    private static Map<Set<ComparableGeneticVariant>, Double> calculateLDMatrix(
            RandomAccessGenotypeData genotypeData, List<Integer> windowSize, double rSquaredThreshold) throws LdCalculatorException, LDHandlerException {
        HashMap<Set<ComparableGeneticVariant>, Double> ldMatrix = new HashMap<>();

        // Get the chromosomes / sequences in the genotype data
        for (String sequence : genotypeData.getSeqNames()) {
            System.out.println("sequence = " + sequence);
            int lastStartPos = 0;

            // Loop through the variants in this sequence in order (low bp to high bp)
            Iterable<GeneticVariant> sequenceVariants = genotypeData.getSequenceGeneticVariants(sequence);

            // Outer loop through variants: loops from first to last variant
            for (GeneticVariant variant : sequenceVariants) {
                // First check if the variants indeed are processed in order
                int startPos = variant.getStartPos();
                if (lastStartPos > startPos) {
                    throw new LDHandlerException(
                            "Genetic variants in the genotype data should be ordered by position");
                }
                lastStartPos = startPos;
                // Now loop through the variants starting from the current variant stopping at the variant
                // for which the starting position difference with the current startpos does not exceed the window size.
                for (GeneticVariant otherVariant : genotypeData.getVariantsByRange(
                        sequence, startPos + 1, startPos + windowSize.get(0) + 1)) {
                    if (!variant.equals(otherVariant)) {
                        // Calculate the R2
                        double rSquared = LdCalculator.calculateRsquare(variant, otherVariant);
                        // If the R2 is above the required threshold, put this pair in the map / matrix
                        if (rSquared >= rSquaredThreshold) {
                            ldMatrix.put( // Create a set as a key because order is not important this way
                                    Collections.unmodifiableSet(new HashSet<>(Arrays.asList(
                                            new ComparableGeneticVariant(variant),
                                            new ComparableGeneticVariant(otherVariant)))),
                                    rSquared);
                        }
                    }
                }
            }
        }
        return ldMatrix;
    }

    public Map<Set<ComparableGeneticVariant>, Double> getLdMatrix() {
        return ldMatrix;
    }

    public double getIndexVariantPValueThreshold() {
        return indexVariantPValueThreshold;
    }

    public List<ComparableGeneticVariant> getIndexVariants() {
        return indexVariants;
    }

    @Override
    public List<Integer> getWindowSize() {
        return windowSize;
    }

    @Override
    public double getRSquaredThreshold() {
        return rSquaredThreshold;
    }

    private static class ComparableGeneticVariant {
        private final GeneticVariant originalVariant;

        private ComparableGeneticVariant(GeneticVariant variant) {
            this.originalVariant = variant;
        }

        GeneticVariant getOriginalVariant() {
            return originalVariant;
        }

        /*
         * (non-Javadoc)
         *
         * @see java.lang.Object#hashCode()
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 0;
            result = prime * result + ((originalVariant.getSequenceName() == null) ? 0 : originalVariant.getSequenceName().hashCode());
            result = prime * result + originalVariant.getStartPos();
            return result;
        }

        /*
         * (non-Javadoc)
         *
         * @see java.lang.Object#equals(java.lang.Object)
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof ComparableGeneticVariant)) {
                return false;
            }
            ComparableGeneticVariant otherComparableGeneticVariant = (ComparableGeneticVariant) obj;
            GeneticVariant other = otherComparableGeneticVariant.getOriginalVariant();

            if (originalVariant.getSequenceName() == null) {
                if (other.getSequenceName() != null) {
                    return false;
                }
            } else if (!originalVariant.getSequenceName().equals(other.getSequenceName())) {
                return false;
            }
            if (originalVariant.getStartPos() != other.getStartPos()) {
                return false;
            }

            // If we get here pos and sequence are identical

            if (originalVariant.getVariantAlleles() == null) {
                if (other.getVariantAlleles() != null) {
                    return false;
                }
            } else if (!originalVariant.getVariantAlleles().equals(other.getVariantAlleles())) {
                return false;
            }
            return true;
        }
    }
}
