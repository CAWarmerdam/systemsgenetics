package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import nl.systemsgenetics.gwassummarystatistics.effectAlleleFilter.EffectAlleleFilter;
import org.molgenis.genotype.RandomAccessGenotypeData;

import javax.annotation.Nonnull;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class FilterableGwasSummaryStatistics implements GwasSummaryStatistics {

    private GwasSummaryStatistics originalGwasSummaryStatistics;
    private EffectAlleleFilter effectAlleleFilter;

    public FilterableGwasSummaryStatistics(GwasSummaryStatistics originalGwasSummaryStatistics,
                                           EffectAlleleFilter effectAlleleFilter) {
        this.originalGwasSummaryStatistics = originalGwasSummaryStatistics;
        this.effectAlleleFilter = effectAlleleFilter;
    }

    @Override
    public String getGwasId() {
        return originalGwasSummaryStatistics.getGwasId();
    }

    @Override
    public Iterable<EffectAllele> getEffectAllelesByRange(String seqName, int rangeStart, int rangeEnd) {
        return filteredIterable(
                originalGwasSummaryStatistics.getEffectAllelesByRange(seqName, rangeStart, rangeEnd));
    }

    @Override
    @Nonnull
    public Iterator<EffectAllele> iterator() {
        return filteredIterable(originalGwasSummaryStatistics).iterator();
    }

    private Iterable<EffectAllele> filteredIterable(Iterable<EffectAllele> originalIterable){
        Iterator<EffectAllele> filteredEffectAlleleIterator = new Iterator<EffectAllele>() {
            private Iterator<EffectAllele> originalIterator = originalIterable.iterator();
            private EffectAllele next;

            @Override
            public boolean hasNext()
            {
                return next != null;
            }

            @Override
            public EffectAllele next() {

                if (next == null)
                {
                    throw new NoSuchElementException();
                }

                EffectAllele currentNext = next;

                // prepare next next
                goToNext();

                return currentNext;
            }

            private void goToNext()
            {
                while (originalIterator.hasNext())
                {
                    EffectAllele originalNext = originalIterator.next();

                    if (!effectAlleleFilter.doesEffectAllelePassFilter(originalNext))
                    {
                        // skip variants on exclude list
                        continue;
                    }

                    next = originalNext;
                    return;
                }
                // We do a return if we find a non excluded next. So if we get here it
                // is the end of the original iterator. Setting next to null so hasNext
                // knows it is the end.
                next = null;
            }
        };

        return () -> filteredEffectAlleleIterator;
    }
}
