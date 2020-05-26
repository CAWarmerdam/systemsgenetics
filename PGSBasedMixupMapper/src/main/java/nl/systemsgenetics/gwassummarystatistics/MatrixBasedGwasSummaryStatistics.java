package nl.systemsgenetics.gwassummarystatistics;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import com.opencsv.CSVWriter;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.GeneticVariantBackedEffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

import javax.annotation.Nonnull;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class MatrixBasedGwasSummaryStatistics implements GwasSummaryStatistics {

    private static final Logger LOGGER = Logger.getLogger(MatrixBasedGwasSummaryStatistics.class);
    private String gwasId;
    private final LinkedHashMap<String, Integer> variantIds;
    private DoubleMatrix1D betaCoefficients;
    private DoubleMatrix1D pValues;
    private final List<String> alleles;
    private final RandomAccessGenotypeData genotypeData;

    public MatrixBasedGwasSummaryStatistics(String gwasId, LinkedHashMap<String, Integer> variantIds,
                                            List<String> alleles,
                                            DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues,
                                            RandomAccessGenotypeData genotypeData) {
        // Set the gwas identifier
        this.gwasId = gwasId;
        this.genotypeData = genotypeData;

        // Make sure that the sizes of the three collections are equal.
        if (variantIds.size() != betaCoefficients.size() || variantIds.size() != pValues.size()) {
            throw new GwasSummaryStatisticsException(
                    String.format("Trying to add summary statistics with unequal number of variants, " +
                                    "betas or p-values (%d vs %d, %d)",
                            variantIds.size(), betaCoefficients.size(), pValues.size()));
        }

        this.variantIds = variantIds;
        this.betaCoefficients = betaCoefficients;
        this.pValues = pValues;
        this.alleles = alleles;
    }

    public MatrixBasedGwasSummaryStatistics(String gwasId, RandomAccessGenotypeData genotypeData) {
        this.gwasId = gwasId;
        this.variantIds = new LinkedHashMap<>();
        this.betaCoefficients = new DenseDoubleMatrix1D(0);
        this.pValues = new DenseDoubleMatrix1D(0);
        this.alleles = new ArrayList<>();
        this.genotypeData = genotypeData;
    }

    public void add(LinkedHashMap<String, Integer> variantIds, List<String> alleles,
                    DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues,
                    Double leastStringentPValueThreshold) {

        assertEqualLengths(variantIds, alleles, betaCoefficients, pValues);

        List<Integer> significantAssociations = new ArrayList<>();
        LinkedHashMap<String, Integer> updatedVariantIds = new LinkedHashMap<>();

        double[] pValueArray = pValues.toArray();

        int updatedIndex = 0;
        for (String variantId : variantIds.keySet()) {
            Integer index = variantIds.get(variantId);
            double pValue = pValueArray[index];
            if (pValue < leastStringentPValueThreshold) {
                significantAssociations.add(index);
                updatedVariantIds.put(variantId, updatedIndex++);
            }
        }

        int[] indicesArray = significantAssociations.stream().mapToInt(Integer::intValue).toArray();

        this.betaCoefficients = DoubleFactory1D.dense.append(
                this.betaCoefficients, betaCoefficients.viewSelection(indicesArray));
        this.pValues = DoubleFactory1D.dense.append(
                this.pValues, pValues.viewSelection(indicesArray));
        this.alleles.addAll(significantAssociations.stream().map(alleles::get)
                .collect(Collectors.toList()));

        int originalSize = this.variantIds.size();

        for (String variantId : updatedVariantIds.keySet()) {
            Integer variantIndex = updatedVariantIds.get(variantId);
            this.variantIds.put(
                    variantId,
                    variantIndex + originalSize);
        }
    }

    public void add(LinkedHashMap<String, Integer> variantIds, List<String> alleles,
                    DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues) {
        assertEqualLengths(variantIds, alleles, betaCoefficients, pValues);

        this.betaCoefficients = DoubleFactory1D.dense.append(this.betaCoefficients, betaCoefficients);
        this.pValues = DoubleFactory1D.dense.append(this.pValues, pValues);
        this.alleles.addAll(alleles);

        int originalSize = this.variantIds.size();

        for (String variantId : variantIds.keySet()) {
            this.variantIds.put(variantId, variantIds.get(variantId) + originalSize);
        }
    }

    private void assertEqualLengths(LinkedHashMap<String, Integer> variantIds, List<String> alleles, DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues) {
        if (variantIds.size() != betaCoefficients.size() || variantIds.size() != pValues.size()) {
            throw new GwasSummaryStatisticsException(
                    String.format("Trying to add summary statistics with unequal number of variants, " +
                            "betas or p-values (%d vs %d, %d)",
                            variantIds.size(), betaCoefficients.size(), pValues.size()));
        }

        assert alleles.size() == variantIds.size();
    }

    @Override
    public String getGwasId() {
        return gwasId;
    }

    @Override
    public Iterable<EffectAllele> getEffectAllelesByRange(String seqName, int rangeStart, int rangeEnd) {
        throw new UnsupportedOperationException("Not yet supported");
    }

    @Override
    @Nonnull
    public Iterator<EffectAllele> iterator() {
        Iterator<Map.Entry<String, Integer>> riskAlleleIterator = variantIds.entrySet().iterator();
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
                goToNext();

                return currentNext;
            }

            private void goToNext()
            {
                while (riskAlleleIterator.hasNext())
                {
                    Map.Entry<String, Integer> provisionalNext = riskAlleleIterator.next();

                    String variantId = provisionalNext.getKey();
                    if (variantIdMap.containsKey(variantId))
                    {
                        // skip variants on exclude list
                        continue;
                    }

                    double pValue = pValues.getQuick(variantIds.get(variantId));
                    double effectSize = betaCoefficients.getQuick(variantIds.get(variantId));

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

    public void save(Path pathPrefix, Double pValueThreshold) throws IOException {
        // Format the output path.
        Path outputPath = Paths.get(String.format("%s.txt",
                pathPrefix));

        LOGGER.info(String.format("Writing summary statistics to '%s'", outputPath));

        if (!outputPath.getParent().toFile().isDirectory() && !outputPath.getParent().toFile().mkdirs()) {
            throw new IOException(String.format("Could not create directory '%s'", outputPath.getParent()));
        }

        // Create a CSV writer without quotation characters and with tabs as separators to mimic
        // the legacy gwas summary statistics files.
        CSVWriter writer = new CSVWriter(new FileWriter(outputPath.toFile()),
                '\t', CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        // Write a header (is skipped when reading again)
        writer.writeNext(new String[]{"variant", "AssessedAllele", "ES", "pvalue"});

        // Write risk entries
        for (Map.Entry<String, Integer> variantAssociation : variantIds.entrySet()) {

            // Check if the pvalue for this association is below the given threshold.
            if (pValues.getQuick(variantAssociation.getValue()) < pValueThreshold) {
                // Only write the association whenever the p-value is below the threshold
                writer.writeNext(new String[]{
                        variantAssociation.getKey(),
                        alleles.get(variantAssociation.getValue()),
                        String.valueOf(betaCoefficients.getQuick(variantAssociation.getValue())),
                        String.valueOf(pValues.getQuick(variantAssociation.getValue()))});
            }
        }

        writer.close();
    }

    public int size() {
        return variantIds.size();
    }
}
