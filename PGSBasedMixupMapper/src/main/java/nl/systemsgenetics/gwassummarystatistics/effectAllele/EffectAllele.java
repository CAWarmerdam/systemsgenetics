package nl.systemsgenetics.gwassummarystatistics.effectAllele;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.*;

/**
 * Represents an effect allele, possibly originating from a genome wide association study.
 * @author Robert Warmerdam
 */
public abstract class EffectAllele implements Comparable<EffectAllele> {

    /**
     * Get the primary variant ID
     *
     * @return String
     */
    public abstract String getPrimaryVariantId();

    /**
     * Get the Sequence this variant is located on
     *
     * @return the Sequence
     */
    public abstract String getSequenceName();

    /**
     * Gets the starting position on the sequence
     *
     * @return int
     */
    public abstract int getStartPos();

    /**
     * Gets the effect size of this allele.
     * Either a GWAS beta coefficient or log odds.
     *
     * @return the effect size as a double
     */
    public abstract double getEffectSize();

    /**
     * Gets the P-value corresponding to the effect.
     *
     * @return the P-value.
     */
    public abstract double getPValue();

    /**
     * Gets the -log10 transformed P-value corresponding to the effect.
     *
     * @return The -log10 transformed P-value.
     */
    public abstract double getLogTransformedPValue();

    /**
     * Gets the allele frequency of this effect allele.
     *
     * @return the frequency of this effect allele.
     */
    public abstract double getAlleleFrequency();

    /**
     * Method that returns if the variant is a SNP or not.
     *
     * @return true if the variant is a snp.
     */
    public abstract boolean variantIsSnp();

    /**
     * Method that indicates if the variant is biallelic.
     *
     * @return true if the variant is biallelic.
     */
    public abstract boolean variantIsBiallelic();

    /**
     * Gets the allele to which the effect corresponds.
     *
     * @return the allele.
     */
    public abstract Allele getAllele();

    /**
     * Returns whether or not the given variant corresponds to this effect allele.
     *
     * @param variant The genetic variant to compare with.
     * @return true if the variant matches at least the primary variant id,
     * sequence name and the starting position.
     */
    public boolean matchesVariant(GeneticVariant variant) {
        return variant != null
                && this.getPrimaryVariantId().equals(variant.getPrimaryVariantId())
                && this.getSequenceName().equals(variant.getSequenceName())
                && this.getStartPos() == variant.getStartPos();
    }

    /**
     * <p>Compares this effect allele with the specified effect allele for order.
     * Returns a negative integer, zero, or a positive integer as this object is less than,
     * equal to, or greater than the specified object.</p>
     *
     * <p>This effect allele is considered less than the specified effect allele if this effect allele's <i>p</i>-value
     * is less than the specified effect allele's <i>p</i>-value, or when these values are equal,
     * if this effect allele's effect size is greater than the specified effect allele's effect size.</p>
     *
     * <p>The opposite is of also true:
     * This effect allele is considered greater than the specified effect allele if this effect allele's <i>p</i>-value
     * is greater than the specified effect allele's <i>p</i>-value, or when these values are equal,
     * if this effect allele's effect size is less than the specified effect allele's effect size.</p>
     *
     * <p>When ordering effect alleles this definition consequently will result in a more significant
     * (or a more affecting) effect allele being ordered prior to a less significant (or a less affecting)
     * effect allele.</p>
     *
     * @param other the effect allele being compared.
     * @return a negative integer, zero, or a positive integer as this effect allele is less than, equal to or greater
     * than the specified object as defined above.
     */
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

    /**
     * Method that sorts an iterator of unordered effect alleles according to the implemented compareTo method,
     * for every sequence.
     *
     * For every sequence the method fills a sorted cache of effect alleles, and empties this one by one.
     * This is done for every consecutive sequence input sequence, and thus, the unordered iterator is at least
     * expected to be grouped by sequence. If this is not the case an exception will be thrown.
     *
     * @param unorderedIterator An iterator of effect alleles grouped by sequence.
     * @return An iterator of effect alleles grouped by sequence.
     * Within every sequence group the effect alleles are sorted.
     */
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
