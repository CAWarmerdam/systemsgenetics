package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.GeneticVariantBackedEffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;

import javax.annotation.Nonnull;
import java.util.*;

public class ReadOnlyGwasSummaryStatistics implements GeneticVariantBackedGwasSummaryStatistics {

    private static final Logger LOGGER = Logger.getLogger(ReadOnlyGwasSummaryStatistics.class);
    private VcfGwasSummaryStatistics vcfGwasSummaryStatistics;
    private int studyIndex;
    private String studyName;

    public ReadOnlyGwasSummaryStatistics(VcfGwasSummaryStatistics vcfGwasSummaryStatistics, String studyName) {

        List<String> sampleNames = Arrays.asList(vcfGwasSummaryStatistics.getStudyNames());
        if (!sampleNames.contains(studyName)) {
            throw new GwasSummaryStatisticsException(
                    String.format("Study name %s does not exist in the given summary statistics", studyName));
        }

        this.vcfGwasSummaryStatistics = vcfGwasSummaryStatistics;
        this.studyName = studyName;
        this.studyIndex = sampleNames.indexOf(studyName);
    }

    @Override
    public String getGwasId() {
        return studyName;
    }

    @Override
    public Iterable<EffectAllele> getEffectAllelesByRange(String seqName, int rangeStart, int rangeEnd) {
        ReadOnlyGwasSummaryStatistics thisSummaryStatistics = this;

        return () -> vcfGwasSummaryStatistics.effectAlleles(
                vcfGwasSummaryStatistics.getVariantsByRange(seqName, rangeStart, rangeEnd),
                thisSummaryStatistics);
    }

    @Override
    public double getEffectSizeEstimate(GeneticVariantBackedEffectAllele effectAllele) {
        return vcfGwasSummaryStatistics.getEffectSizeEstimates(effectAllele.getVariant())[studyIndex][effectAllele.getAlleleIndex()];
    }

    @Override
    public double getLogTransformedPValue(GeneticVariantBackedEffectAllele effectAllele) {
        return vcfGwasSummaryStatistics.getTransformedPValues(effectAllele.getVariant())[studyIndex][effectAllele.getAlleleIndex()];
    }

    @Override
    public double getPValue(GeneticVariantBackedEffectAllele effectAllele) {
        return Math.pow(10, -this.getLogTransformedPValue(effectAllele));
    }

    @Override
    public double getAlleleFrequency(GeneticVariantBackedEffectAllele effectAllele) {
        return vcfGwasSummaryStatistics.getAlleleFrequency(effectAllele.getVariant())[studyIndex][effectAllele.getAlleleIndex()];
    }

    @Override
    @Nonnull
    public Iterator<EffectAllele> iterator() {
        return vcfGwasSummaryStatistics.effectAlleles(
                vcfGwasSummaryStatistics.variantIterator(), this);
    }

}
