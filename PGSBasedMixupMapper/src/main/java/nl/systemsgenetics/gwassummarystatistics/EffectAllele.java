package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.variant.GeneticVariant;

import java.util.*;

public class EffectAllele implements Comparable<EffectAllele> {
    private final int alleleIndex;
    private final GeneticVariant variant;
    private final MultiStudyGwasSummaryStatistics summaryStatistics;
    private final int studyIndex;

    private EffectAllele(GeneticVariant variant, MultiStudyGwasSummaryStatistics summaryStatistics,
                         int alleleIndex, int studyIndex) {
        this.variant = variant;
        this.alleleIndex = alleleIndex;
        this.studyIndex = studyIndex;
        this.summaryStatistics = summaryStatistics;
    }

    private static EffectAllele fromVariant(GeneticVariant variant, MultiStudyGwasSummaryStatistics summaryStatistics,
                                            int alleleIndex, int studyIndex) {
        return new EffectAllele(
                variant,
                summaryStatistics,
                alleleIndex,
                studyIndex);
    }

    public static Iterator<EffectAllele> effectAlleles(Iterator<GeneticVariant> variantIterator, MultiStudyGwasSummaryStatistics summaryStatistics, int studyIndex) {
        return new Iterator<EffectAllele>() {
            GeneticVariant variant = variantIterator.hasNext() ? variantIterator.next() : null;
            int currentAlleleIndex = 0;

            @Override
            public boolean hasNext() {
                if (variant == null) return false;
                while (currentAlleleIndex == variant.getAlternativeAlleles().getAlleleCount()) {
                    if (!variantIterator.hasNext()) {
                        return false;
                    }
                    variant = variantIterator.next();
                    currentAlleleIndex = 0;
                }
                return true;
            }

            @Override
            public EffectAllele next() {
                while (currentAlleleIndex == variant.getAlternativeAlleles().getAlleleCount()) {
                    if (!variantIterator.hasNext()) {
                        throw new NoSuchElementException();
                    }
                    variant = variantIterator.next();
                    currentAlleleIndex = 0;
                }
                return EffectAllele.fromVariant(
                        variant, summaryStatistics, currentAlleleIndex++, studyIndex);
            }
        };
    }

    public static Iterator<EffectAllele> sortEffectAllelesPerSequence(Iterator<EffectAllele> unorderedIterator) {
        return new Iterator<EffectAllele>() {
            TreeSet<EffectAllele> cachedSortedEffectAlleles = new TreeSet<>();
            List<String> previousSequenceNames = new ArrayList<>();
            EffectAllele futureAllele = unorderedIterator.hasNext() ? unorderedIterator.next() : null;

            private TreeSet<EffectAllele> getCachedSortedEffectAlleles(String currentSequenceName) {

                // Initialize a TreeSet for storing effect alleles for a particular sequence
                TreeSet<EffectAllele> sortedEffectAlleles = new TreeSet<>();

                // If the future allele was set sometime in the past, add this to the TreeSet
                if (futureAllele != null) {
                    sortedEffectAlleles.add(futureAllele);
                    futureAllele = null;
                }

                // While the unordered iterator has a next item,
                // loop through the effect alleles.
                while (unorderedIterator.hasNext()) {
                    // Get the sequence name for every effect allele
                    EffectAllele allele = unorderedIterator.next();
                    String sequenceName = allele.variant.getSequenceName();

                    if (sequenceName.equals(currentSequenceName)) {
                        sortedEffectAlleles.add(allele);

                    }
                    // If the effect allele corresponds to a different sequence than the previous
                    // effect alleles do something else.
                    else if (!previousSequenceNames.contains(sequenceName)) {
                        // If the new sequence name is not already encountered,
                        // assume that the previous sequence is done, add this current (now previous)
                        // sequence name to the list of previous sequence names.
                        previousSequenceNames.add(currentSequenceName);
                        // This allele still has to be processed. Do that in a future call to this method
                        // by setting it as the future allele.
                        futureAllele = allele;
                        // Break out of the while loop to return the sorted TreeSet of effect alleles for this sequence
                        break;
                    } else {
                        // If the new sequence name has been observed previously,
                        // this sorting method is not s
                        throw new IllegalArgumentException("Expected the input iterator to be grouped by sequence.");
                    }
                }
                return sortedEffectAlleles;
            }

            @Override
            public boolean hasNext() {
                if (!cachedSortedEffectAlleles.isEmpty()) {
                    return true;
                }
                return unorderedIterator.hasNext() || futureAllele != null;
            }

            @Override
            public EffectAllele next() {

                if (!cachedSortedEffectAlleles.isEmpty()) {
                    // There is still at least one effect allele in the cache,
                    // We can just return this now
                    return cachedSortedEffectAlleles.first();

                } else if (unorderedIterator.hasNext() || futureAllele != null) {
                    // The cache is empty!
                    // Fill the cache for the upcoming sequence.
                    cachedSortedEffectAlleles = getCachedSortedEffectAlleles(
                            futureAllele.variant.getSequenceName());
                    return cachedSortedEffectAlleles.first();
                } else {
                    throw new NoSuchElementException();
                }
            }
        };
    }

    public GeneticVariant getVariant() {
        return variant;
    }

    public double getEffectSize() {
        return summaryStatistics.getEffectSizeEstimates(variant)[studyIndex][alleleIndex];
    }

    public double getStandardErrorOfES() {
        return summaryStatistics.getStandardErrorOfES(variant)[studyIndex][alleleIndex];
    }

    public double getLogTransformedPValue() {
        return summaryStatistics.getTransformedPValues(variant)[studyIndex][alleleIndex];
    }

    @Override
    public int compareTo(EffectAllele other) {
        if (this.getLogTransformedPValue() > other.getLogTransformedPValue()){
            return 1;
        } else if (this.getLogTransformedPValue() == other.getLogTransformedPValue()) {
            return Double.compare(Math.abs(this.getEffectSize()), Math.abs(other.getEffectSize()));
        } else {
            return -1;
        }
    }
}
