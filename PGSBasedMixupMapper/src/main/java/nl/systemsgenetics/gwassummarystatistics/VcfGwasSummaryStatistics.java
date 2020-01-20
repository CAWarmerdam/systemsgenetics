package nl.systemsgenetics.gwassummarystatistics;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMetaFormat;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class VcfGwasSummaryStatistics extends VcfGenotypeData{
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

    public VcfGwasSummaryStatistics(File bzipVcfFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException, VcfGwasSummaryStatisticsException {
        this(bzipVcfFile, new File(bzipVcfFile.getAbsolutePath() + ".tbi"), cacheSize, minimumPosteriorProbabilityToCall);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException, VcfGwasSummaryStatisticsException {
        this(bzipVcfFile, tabixIndexFile, 100, minimumPosteriorProbabilityToCall);
    }

    public VcfGwasSummaryStatistics(File bzipVcfFile, File tabixIndexFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException, VcfGwasSummaryStatisticsException {
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

        // Check if variant IDs are unique
//        long numberOfUniqueVariantIds = StreamSupport.stream(
//                Spliterators.spliteratorUnknownSize(vcfGwasSummaryStatistics.iterator(), Spliterator.ORDERED),
//                false).map(GeneticVariant::getPrimaryVariantId).distinct().count();
//        int numberOfVariants = Iterators.size(vcfGwasSummaryStatistics.iterator());
//        if (numberOfUniqueVariantIds != numberOfVariants) {
//            throw new GwasSummaryStatisticsVcfDataException(String.format("Encountered an error in the VCF file, " +
//                    "VCF file does not contain exclusively unique variant identifiers. (%d vs %d)",
//                    numberOfUniqueVariantIds, numberOfVariants));
//        }

        // Check if required fields ES and SE are in the format column
        assertReservedKeyValidity();
    }

    private void assertReservedKeyValidity() throws VcfGwasSummaryStatisticsException {
        if (this.vcfMeta.getFormatMeta("ES") == null) {
            throw new VcfGwasSummaryStatisticsException("Required field ES not present");
        }
        if (!this.vcfMeta.getFormatMeta("ES").equals(RESERVED_KEYS.get("ES"))) {
            throw new VcfGwasSummaryStatisticsException("Required format field ES not according to specifications");
        }
        if (this.vcfMeta.getFormatMeta("SE") == null) {
            throw new VcfGwasSummaryStatisticsException("Required field SE not present");
        }
        if (!this.vcfMeta.getFormatMeta("SE").equals(RESERVED_KEYS.get("SE"))) {
            throw new VcfGwasSummaryStatisticsException("Requried format field SE not according to specifications");
        }
        assertVcfMetaValidity("LP");
        assertVcfMetaValidity("AF");
        assertVcfMetaValidity("SS");
        assertVcfMetaValidity("EZ");
        assertVcfMetaValidity("SI");
        assertVcfMetaValidity("NC");
        assertVcfMetaValidity("ID");
    }

    private void assertVcfMetaValidity(String formatFieldKey) throws VcfGwasSummaryStatisticsException {
        if (this.vcfMeta.getFormatMeta(formatFieldKey) != null &&
                !this.vcfMeta.getFormatMeta(formatFieldKey).equals(RESERVED_KEYS.get(formatFieldKey))) {
            throw new VcfGwasSummaryStatisticsException(String.format(
                    "Reserved format field '%s' not according to specifications", formatFieldKey));
        }
    }

    public float[][] getEffectSizeEstimates(GeneticVariant variant) throws VcfGwasSummaryStatisticsException {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("ES"));
    }

    public float[][] getStandardErrorOfES(GeneticVariant variant) throws VcfGwasSummaryStatisticsException {
        return getSummaryStatisticsPerAlternativeAllele(variant, getReservedKeyFormat("SE"));
    }

    public float[][] getPValues(GeneticVariant variant) throws VcfGwasSummaryStatisticsException {
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
        if (idx != -1) {
            // retrieve values from sample info
            int i = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String valueString = vcfSample.getData(idx);
                if (valueString != null) {
                    String[] splitValuesString = StringUtils.split(valueString, ',');
                    if (splitValuesString.length != alternativeAlleleCount) {
                        throw new VcfGwasSummaryStatisticsException(String.format(
                                "Error in '%s' value for study [%s], found %d value(s) (%s), " +
                                        "while %d were expected based on the alternative allele count",
                                fieldFormat.getId(), vcfMeta.getSampleName(i), splitValuesString.length,
                                valueString, alternativeAlleleCount));
                    }

                    for (int j = 0; j < splitValuesString.length; j++) {
                        try {
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
}
