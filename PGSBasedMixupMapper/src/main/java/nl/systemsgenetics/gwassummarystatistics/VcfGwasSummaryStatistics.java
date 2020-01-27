package nl.systemsgenetics.gwassummarystatistics;

import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMetaFormat;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class VcfGwasSummaryStatistics extends VcfGenotypeData implements MultiStudyGwasSummaryStatistics{
    private static final Map<String, VcfMetaFormat> RESERVED_KEYS = getReservedKeys();

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

    public VcfGwasSummaryStatistics(File bzipVcfFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException, VcfGwasSummaryStatisticsException {
        this(bzipVcfFile, new File(bzipVcfFile.getAbsolutePath() + ".tbi"), cacheSize, minimumPosteriorProbabilityToCall);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, double minimumPosteriorProbabilityToCall) throws IOException, VcfGwasSummaryStatisticsException {
        this(bzipVcfFile, tabixIndexFile, 100, minimumPosteriorProbabilityToCall);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException, VcfGwasSummaryStatisticsException {
        super(bzipVcfFile, tabixIndexFile, cacheSize, minimumPosteriorProbabilityToCall);

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

    private void assertReservedKeyValidity() throws VcfGwasSummaryStatisticsException {
        for (String reservedKey : RESERVED_KEYS.keySet()) {
            assertVcfMetaValidity(reservedKey);
        }
    }

    private void assertVcfMetaValidity(String formatFieldKey) throws VcfGwasSummaryStatisticsException {
        if (this.vcfMeta.getFormatMeta(formatFieldKey) == null) {
            throw new VcfGwasSummaryStatisticsException(String.format(
                    "Reserved format field '%s' not in format declerations", formatFieldKey));
        }
        if (!this.vcfMeta.getFormatMeta(formatFieldKey).equals(RESERVED_KEYS.get(formatFieldKey))) {
            throw new VcfGwasSummaryStatisticsException(String.format(
                    "Reserved format field '%s' not according to specifications", formatFieldKey));
        }
    }

    @Override
    public float[][] getEffectSizeEstimates(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("ES"));
    }

    public float[][] getStandardErrorOfES(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("SE"));
    }

    @Override
    public float[][] getTransformedPValues(GeneticVariant variant) {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("LP"));
    }

    public float[][] getSummaryStatisticsPerAlternativeAllele(GeneticVariant variant, VcfMetaFormat fieldFormat)
            throws VcfGwasSummaryStatisticsException {
        if (!fieldFormat.getType().equals(VcfMetaFormat.Type.FLOAT)) {
            throw new IllegalArgumentException(String.format("The '%s' field's type '%s' is not 'Float'",
                    fieldFormat.getName(), fieldFormat.getType()));
        }
        if (!fieldFormat.getNumber().equals("A")) {
            throw new IllegalArgumentException(String.format("The '%s' field's number of values '%s' is not 'A'",
                    fieldFormat.getName(), fieldFormat.getNumber()));
        }

        // Get the vcf record
        VcfRecord vcfRecord = getVcfRecord(variant);

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
                        throw new VcfGwasSummaryStatisticsException(String.format(
                                "Error in '%s' value for study [%s], found %d value(s) (%s), " +
                                        "while %d were expected based on the alternative allele count",
                                fieldFormat.getId(), vcfMeta.getSampleName(i), splitValuesString.length,
                                valueString, alternativeAlleleCount));
                    }

                    // For every value, try to convert it to a float
                    for (int j = 0; j < splitValuesString.length; j++) {
                        try {
                            // A dot (".") represents a missing value, we replace this with zero.
                            if (splitValuesString[j].equals(".")) {
                                values[i][j] = 0f;
                            } else {
                                values[i][j] = Float.parseFloat(splitValuesString[j]);
                            }
                        } catch (NumberFormatException e) {
                            throw new VcfGwasSummaryStatisticsException(String.format(
                                    "Error in '%s' value for study [%s], found value: %s",
                                    fieldFormat.getId(), vcfMeta.getSampleName(i), valueString));
                        }
                    }
                }
                ++i;
            }
        }
        return values;
    }

    public static VcfMetaFormat getReservedKeyFormat(String reservedKey) {
        return RESERVED_KEYS.get(reservedKey);
    }

    @Override
    public Iterator<EffectAllele> effectAlleles(String studyName) {
        return EffectAllele.effectAlleles(this.iterator(), this, Arrays.asList(getSampleNames()).indexOf(studyName));
    }
}
