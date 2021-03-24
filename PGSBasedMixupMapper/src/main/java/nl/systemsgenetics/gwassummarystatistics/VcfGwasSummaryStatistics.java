package nl.systemsgenetics.gwassummarystatistics;

import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.GeneticVariantBackedEffectAllele;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMetaFormat;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Object representing a VCF format as defined by https://gwas.mrcieu.ac.uk/about/
 * @author Robert Warmerdam
 */
public class VcfGwasSummaryStatistics implements Closeable {
    private static final Map<String, VcfMetaFormat> RESERVED_KEYS = getReservedKeys();
    private VcfGenotypeData vcfGenotypeData;
    private RandomAccessGenotypeData genotypeData;

    private static Map<String, VcfMetaFormat> getReservedKeys() {
        HashMap<String, VcfMetaFormat> reservedKeys = new HashMap<>();
        reservedKeys.put("ES", createFieldProperties(
                "ES", "A", "Float", "Effect size estimate relative to the alternative allele"));
        reservedKeys.put("SE", createFieldProperties(
                "SE", "A", "Float", "Standard error of effect size estimate"));
        reservedKeys.put("LP", createFieldProperties(
                "LP", "A", "Float", "-log10 p-value for effect estimate"));
        reservedKeys.put("AF", createFieldProperties(
                "AF", "A", "Float", "Alternate allele frequency in the association study"));
        reservedKeys.put("SS", createFieldProperties(
                "SS", "A", "Float", "Sample size used to estimate genetic effect"));
        reservedKeys.put("EZ", createFieldProperties(
                "EZ", "A", "Float", "Z-score provided if it was used to derive the EFFECT and SE fields"));
        reservedKeys.put("SI", createFieldProperties(
                "SI", "A", "Float", "Accuracy score of summary data imputation"));
        reservedKeys.put("NC", createFieldProperties(
                "NC", "A", "Float", "Number of cases used to estimate genetic effect"));
        reservedKeys.put("ID", createFieldProperties(
                "ID", "1", "String", "Study variant identifier"));
        return reservedKeys;
    }

    private static VcfMetaFormat createFieldProperties(String id, String number, String type, String description) {
        LinkedHashMap<String, String> properties = new LinkedHashMap<>();
        properties.put("ID", id);
        properties.put("Number", number);
        properties.put("Type", type);
        properties.put("Description", description);
        return new VcfMetaFormat(properties);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, int cacheSize) throws IOException, GwasSummaryStatisticsException {
        this(bzipVcfFile, new File(bzipVcfFile.getAbsolutePath() + ".tbi"), cacheSize);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile) throws IOException, GwasSummaryStatisticsException {
        this(bzipVcfFile, tabixIndexFile, 100);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, int cacheSize) throws IOException, GwasSummaryStatisticsException {
        vcfGenotypeData = new VcfGenotypeData(bzipVcfFile, tabixIndexFile, cacheSize,
                0.4);// Represents the minimum posterior probability to call,
        // 0.4 is generally the default value);
        genotypeData = vcfGenotypeData;

        // Check the following
        // Variant IDs must be unique!
        // Reserved keys:
        // Field 	Description 	Required
        // ES 	Effect size estimate relative to the alternative allele 	YES!
        // SE 	Standard error of effect size estimate 	YES!
        // LP 	-log10 p-value for effect estimate 	NO
        // AF 	Alternate allele frequency in the association study 	NO
        // SS 	Sample size used to estimate genetic effect 	NO
        // EZ 	Z-score provided if it was used to derive the EFFECT and SE fields 	NO
        // SI 	Accuracy score of summary data imputation 	NO
        // NC 	Number of cases used to estimate genetic effect 	NO
        // ID 	Study variant identifier 	NO

        // Check if required fields ES and SE are in the format column
        assertReservedKeyValidity();
    }

    public void applyVariantFilter(VariantFilter variantFilter) {
        genotypeData = new VariantFilterableGenotypeDataDecorator(vcfGenotypeData, variantFilter);
    }

    private void assertReservedKeyValidity() throws GwasSummaryStatisticsException {
        for (String reservedKey : RESERVED_KEYS.keySet()) {
            assertVcfMetaValidity(reservedKey);
        }
    }

    private void assertVcfMetaValidity(String formatFieldKey) throws GwasSummaryStatisticsException {
        if (vcfGenotypeData.getFormatMeta(formatFieldKey) == null) {
            throw new GwasSummaryStatisticsException(String.format(
                    "Reserved format field '%s' not in format declerations", formatFieldKey));
        }
        if (!vcfGenotypeData.getFormatMeta(formatFieldKey).equals(RESERVED_KEYS.get(formatFieldKey))) {
            throw new GwasSummaryStatisticsException(String.format(
                    "Reserved format field '%s' not according to specifications", formatFieldKey));
        }
    }

    public float[][] getEffectSizeEstimates(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("ES"));
    }

    public float[][] getStandardErrorOfES(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("SE"));
    }

    public float[][] getTransformedPValues(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("LP"));
    }

    public float[][] getAlleleFrequency(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("AF"));
    }

    public String[] getStudyNames() {
        return vcfGenotypeData.getSampleNames();
    }

    public float[][] getSummaryStatisticsPerAlternativeAllele(GeneticVariant variant, VcfMetaFormat fieldFormat)
            throws GwasSummaryStatisticsException {
        if (!fieldFormat.getType().equals(VcfMetaFormat.Type.FLOAT)) {
            throw new IllegalArgumentException(String.format("The '%s' field's type '%s' is not 'Float'",
                    fieldFormat.getName(), fieldFormat.getType()));
        }
        if (!fieldFormat.getNumber().equals("A")) {
            throw new IllegalArgumentException(String.format("The '%s' field's number of values '%s' is not 'A'",
                    fieldFormat.getName(), fieldFormat.getNumber()));
        }

        // Get the vcf record
        VcfRecord vcfRecord = vcfGenotypeData.getVcfRecord(variant);

        // Get the number of studies
        final int nrStudies = vcfRecord.getNrSamples();
        if (nrStudies == 0) {
            return new float[0][];
        }

        // Get the number of alleles
        int alternativeAlleleCount = variant.getAlternativeAlleles().getAlleleCount();

        // Initialize values with zeros
        float[][] values = new float[nrStudies][alternativeAlleleCount];

        // Get the index of the required field
        int idx = vcfRecord.getFormatIndex(fieldFormat.getId());
        // Check if the field is present
        if (idx != -1) {
            // retrieve values from sample info
            int i = 0;
            // Get the value for every sample
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String valueString = vcfSample.getData(idx);
                if (valueString != null) {
                    // There should be the same number of values as alleles. These should be split by a comma (",")
                    String[] splitValuesString = StringUtils.split(valueString, ',');
                    // Check if the expected number of values corresponds to the actual number of values,
                    // and throw an exception if this is not the case.
                    if (splitValuesString.length != alternativeAlleleCount) {
                        throw new GwasSummaryStatisticsException(String.format(
                                "Error in '%s' value for study [%s], found %d value(s) (%s), " +
                                        "while %d were expected based on the alternative allele count",
                                fieldFormat.getId(), vcfGenotypeData.getSampleNames()[i], splitValuesString.length,
                                valueString, alternativeAlleleCount));
                    }

                    // For every value, try to convert it to a float
                    for (int j = 0; j < splitValuesString.length; j++) {
                        try {
                            // A dot (".") represents a missing value, we replace this with zero.
                            if (splitValuesString[j].equals(".")) {
                                values[i][j] = 0f;
                            } else if (splitValuesString[j].equals("inf")) {
                                values[i][j] = Float.POSITIVE_INFINITY;
                            } else {
                                values[i][j] = Float.parseFloat(splitValuesString[j]);
                            }
                        } catch (NumberFormatException e) {
                            List<String> variantIdentifiers = vcfRecord.getIdentifiers();
                            String variantIdentifier = variantIdentifiers != null ? variantIdentifiers.get(0) : "-";
                            variantIdentifier = StringUtils.abbreviate(variantIdentifier, 64);

                            throw new GwasSummaryStatisticsException(String.format(
                                    "Error in '%s' value for study [%s], variant '%s', found value: %s",
                                    fieldFormat.getId(), vcfGenotypeData.getSampleNames()[i],
                                    variantIdentifier, valueString));
                        }
                    }
                }
                ++i;
            }
        }
        return values;
    }

    static VcfMetaFormat getReservedKeyFormat(String reservedKey) {
        return RESERVED_KEYS.get(reservedKey);
    }

    public Iterator<EffectAllele> effectAlleles(Iterator<GeneticVariant> variantIterator,
                                                ReadOnlyGwasSummaryStatistics summaryStatistics) {

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
            public GeneticVariantBackedEffectAllele next() {
                while (currentAlleleIndex == variant.getAlternativeAlleles().getAlleleCount()) {
                    if (!variantIterator.hasNext()) {
                        throw new NoSuchElementException();
                    }
                    variant = variantIterator.next();
                    currentAlleleIndex = 0;
                }
                return new GeneticVariantBackedEffectAllele(variant, summaryStatistics, currentAlleleIndex++);
            }
        };
    }

    public Map<String, GeneticVariant> getVariantIdMap() {
        return genotypeData.getVariantIdMap();
    }

    public Iterator<GeneticVariant> variantIterator() {
        return genotypeData.iterator();
    }

    public Iterator<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
        return vcfGenotypeData.getVariantsByRange(seqName, rangeStart, rangeEnd).iterator();
    }

    @Override
    public void close() throws IOException {
        genotypeData.close();
    }
}
