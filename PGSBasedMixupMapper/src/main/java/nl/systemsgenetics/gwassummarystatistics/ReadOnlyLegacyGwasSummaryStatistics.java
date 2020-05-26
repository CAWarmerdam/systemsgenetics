package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.polygenicscorecalculator.Main;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.regex.Pattern;

public class ReadOnlyLegacyGwasSummaryStatistics implements GwasSummaryStatistics {

    private static final Logger LOGGER = Logger.getLogger(ReadOnlyLegacyGwasSummaryStatistics.class);
    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
    private final TextFile readFiles;
    private final File riskFilePath;
    private String gwasId;
    private RandomAccessGenotypeData genotypeData;

    public ReadOnlyLegacyGwasSummaryStatistics(String gwasId, File riskFilePath,
                                               RandomAccessGenotypeData genotypeData) throws IOException {
        this.gwasId = gwasId;
        this.riskFilePath = riskFilePath;
        this.readFiles = new TextFile(riskFilePath.getAbsolutePath(), TextFile.R);
        this.genotypeData = genotypeData;
    }

    @Override
    public String getGwasId() {
        return gwasId;
    }

    @Override
    public Iterable<EffectAllele> getEffectAllelesByRange(String seqName, int rangeStart, int rangeEnd) {
        throw new UnsupportedOperationException("Not yet supported.");
    }

    @Override
    public Iterator<EffectAllele> iterator() {

        Map<String, GeneticVariant> variantIdMap = genotypeData.getVariantIdMap();

        return new Iterator<EffectAllele>() {
            private EffectAllele next;

            @Override
            public boolean hasNext() {
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
                try {
                    goToNext();
                } catch (IOException e) {
                    next = null;
                }

                return currentNext;
            }

            private void goToNext() throws IOException {

                readFiles.readLine();
                String s;
                while ((s = readFiles.readLine()) != null) {
                    String[] parts = TAB_PATTERN.split(s);

                    String variantId = parts[0];
                    if (variantIdMap.containsKey(variantId))
                    {
                        // skip variants on exclude list
                        continue;
                    }

                    double pValue = Double.parseDouble(parts[3]);
                    double effectSize = Double.parseDouble(parts[2]);

                    GeneticVariant variant = variantIdMap.get(variantId);
                    next = new RiskEntry(variantId, variant.getSequenceName(), variant.getStartPos(),
                            variant.getVariantAlleles().get(variant.getAlleleCount() - 1).getAlleleAsSnp(),
                            effectSize, pValue);
                    return;
                }
                // We do a return if we find a non excluded next. So if we get here it
                // is the end of the original iterator. Setting next to null so hasNext
                // knows it is the end.
                next = null;
            }
        };
    }
}
