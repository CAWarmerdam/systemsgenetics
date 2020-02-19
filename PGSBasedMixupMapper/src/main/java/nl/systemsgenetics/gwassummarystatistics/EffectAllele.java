package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.*;

public abstract class EffectAllele implements Comparable<EffectAllele> {

    public abstract String getSequenceName();

    public abstract int getStartPos();

    public abstract double getEffectSize();

    public abstract double getPValue();

    public abstract double getLogTransformedPValue();

    public abstract Allele getAllele();

    public abstract String getPrimaryVariantId();

    public boolean matchesVariant(GeneticVariant variant) {
        return variant != null
                && this.getPrimaryVariantId().equals(variant.getPrimaryVariantId())
                && this.getSequenceName().equals(variant.getSequenceName())
                && this.getStartPos() == variant.getStartPos();
    }

    @Override
    public int compareTo(EffectAllele other) {
        if (this.getPValue() < other.getPValue()){
            return -1;
        } else if (this.getPValue() == other.getPValue()) {
            return Double.compare(
                    Math.abs(other.getEffectSize()),
                    Math.abs(this.getEffectSize()));
        } else {
            return 1;
        }
    }

    @Override
    public String toString() {
        return String.join("\t",
                getPrimaryVariantId(),
                getSequenceName(),
                String.valueOf(getStartPos()),
                String.valueOf(getAllele()),
                String.valueOf(getEffectSize()),
                String.valueOf(getPValue()));

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
                    String sequenceName = allele.getSequenceName();

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
                    return cachedSortedEffectAlleles.pollFirst();

                } else if (unorderedIterator.hasNext() || futureAllele != null) {
                    // The cache is empty!
                    // Fill the cache for the upcoming sequence.
                    cachedSortedEffectAlleles = getCachedSortedEffectAlleles(
                            futureAllele.getSequenceName());
                    return cachedSortedEffectAlleles.pollFirst();
                } else {
                    throw new NoSuchElementException();
                }
            }
        };
    }
}
