package nl.systemsgenetics.gwassummarystatistics;

import org.molgenis.genotype.Sample;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.molgenis.vcf.meta.VcfMetaFormat;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.*;

import static org.testng.Assert.*;

public class VcfGwasSummaryStatisticsTest {

    private File exampleGwasVcfFile = new File(this.getClass()
            .getResource("/summarystatistics/gwas_vcf_sample.vcf.gz").toURI());

    public VcfGwasSummaryStatisticsTest() throws URISyntaxException {
    }


    @org.testng.annotations.Test
    public void testGetSummaryStatistics() throws IOException, VcfGwasSummaryStatisticsException {
        List<String> expectedVariantIds = new ArrayList<>(Arrays.asList(
                "rs10399793", "rs2462492", "rs114608975", "rs6702460", "rs8179466",
                "rs6680723", "rs12025928", "rs12238997", "rs72631875", "rs12029736", "rs9442385"));


        VcfGwasSummaryStatistics summaryStatistics =
                new VcfGwasSummaryStatistics(exampleGwasVcfFile,
                750000, // Represents the cache size in number of variants,
                // value copied from the GeneticRiskScoreCalculator module
                0.4);// Represents the minimum posterior probability to call,
                // 0.4 is generally the default value);

        assertEquals(summaryStatistics.getSamples(),
                Collections.singletonList(new Sample("IEU-b-1", null, null)));

        List<GeneticVariant> variants = new ArrayList<>();

        int variantIndex = 0;
        for (GeneticVariant variant : summaryStatistics) {
            assertEquals(variant.getPrimaryVariantId(), expectedVariantIds.get(variantIndex++));
            variants.add(variant);
        }

        float[][] effectSizeEstimates = summaryStatistics.getEffectSizeEstimates(
                variants.get(0));
        assertEquals(effectSizeEstimates, new float[][]{{-0.00457955f}});
        float[][] standardErrorOfES = summaryStatistics.getStandardErrorOfES(
                variants.get(0));
        assertEquals(standardErrorOfES, new float[][]{{0.00443769f}});
        float[][] pValues = summaryStatistics.getTransformedPValues(
                variants.get(0));
        assertEquals(pValues, new float[][]{{0.522879f}});
        float[][] ncs = summaryStatistics.getSummaryStatisticsPerAlternativeAllele(
                variants.get(0), VcfGwasSummaryStatistics.getReservedKeyFormat("NC"));
        assertEquals(ncs, new float[1][1]);
        float[][] effectSizeEstimatesLastVar = summaryStatistics.getEffectSizeEstimates(
                variants.get(10));
        assertEquals(effectSizeEstimatesLastVar, new float[][]{{0.00308509f, 0.00224527f}});
        float[][] standardErrorOfESLastVar = summaryStatistics.getStandardErrorOfES(
                variants.get(10));
        assertEquals(standardErrorOfESLastVar, new float[][]{{0.0152548f, 0.00492936f}});
        float[][] pValuesLastVar = summaryStatistics.getTransformedPValues(
                variants.get(10));
        assertEquals(pValuesLastVar, new float[][]{{0.0757207f, 0.187087f}});
        float[][] accuracyScore = summaryStatistics.getSummaryStatisticsPerAlternativeAllele(
                variants.get(10), VcfGwasSummaryStatistics.getReservedKeyFormat("SI"));
        assertEquals(accuracyScore, new float[][]{{0f, 0f}});

        // Getting a string should throw an exception
        try {
            summaryStatistics.getSummaryStatisticsPerAlternativeAllele(
                    variants.get(10), VcfGwasSummaryStatistics.getReservedKeyFormat("ID"));
            fail("getSummaryStatisticsPerAlternativeAllele() did not raise an " +
                    "error while a String was requested");
        } catch (IllegalArgumentException e) {
            assertEquals(e.getMessage(), "The 'FORMAT' field's type 'String' is not 'Float'");
        }

        // In the 10th variant (index 9) the allele count does not correspond to the numbers given
        // this should therefore give an exception
        try {
            summaryStatistics.getEffectSizeEstimates(
                    variants.get(9));
            fail("getEffectSizeEstimates() did not raise an " +
                    "exception while the alternative allele count " +
                    "does not correspond to the number of summary statistic values");
        } catch (VcfGwasSummaryStatisticsException e) {
            assertEquals(e.getMessage(),
                    "Error in 'ES' value for study [IEU-b-1], " +
                            "found 1 value(s) (-0.0042266), while 2 were " +
                            "expected based on the alternative allele count");
        }
    }

    @org.testng.annotations.Test
    public void testGetReservedKeyFormat() {
        LinkedHashMap<String, String> propertiesEffectSize = new LinkedHashMap<>();
        propertiesEffectSize.put("ID", "ES");
        propertiesEffectSize.put("Number", "A");
        propertiesEffectSize.put("Type", "Float");
        propertiesEffectSize.put("Description", "Effect size estimate relative to the alternative allele");
        VcfMetaFormat expectedFormatEffectSize = new VcfMetaFormat(propertiesEffectSize);

        LinkedHashMap<String, String> propertiesStandardError = new LinkedHashMap<>();
        propertiesStandardError.put("ID", "SE");
        propertiesStandardError.put("Number", "A");
        propertiesStandardError.put("Type", "Float");
        propertiesStandardError.put("Description", "Standard error of effect size estimate");
        VcfMetaFormat expectedFormatStandardError = new VcfMetaFormat(propertiesStandardError);

        VcfMetaFormat actualFormatEffectSize = VcfGwasSummaryStatistics.getReservedKeyFormat("ES");
        assertEquals(actualFormatEffectSize, expectedFormatEffectSize);

        VcfMetaFormat actualFormatStandardError = VcfGwasSummaryStatistics.getReservedKeyFormat("SE");
        assertEquals(actualFormatStandardError, expectedFormatStandardError);
    }

    @org.testng.annotations.Test
    public void testVariantFilterableGwasSummaryStatisticsDecorator() throws IOException {
        Set<String> variantIdsToInclude = new HashSet<>(Arrays.asList(
                "rs10399793", "rs8179466", "rs9442385"));
        float[][][] effectSizes = {{{-0.00457955f}}, {{0.00336766f}}, {{0.00308509f, 0.00224527f}}};

        VcfGwasSummaryStatistics summaryStatistics =
                new VcfGwasSummaryStatistics(exampleGwasVcfFile,
                        750000, // Represents the cache size in number of variants,
                        // value copied from the GeneticRiskScoreCalculator module
                        0.4);// Represents the minimum posterior probability to call,
        // 0.4 is generally the default value);

        MultiStudyGwasSummaryStatistics filteredVariants = new VariantFilterableGwasSummaryStatisticsDecorator(
                summaryStatistics,
                new VariantIdIncludeFilter(variantIdsToInclude));

        int filteredVariantIndex = 0;
        for (GeneticVariant variant : filteredVariants) {
            assertTrue(variantIdsToInclude.contains(variant.getPrimaryVariantId()));
            float[][] effectSizeEstimates = summaryStatistics.getEffectSizeEstimates(variant);
            float[][] effectSizeEstimatesfiltered = filteredVariants.getEffectSizeEstimates(variant);
            assertEquals(effectSizeEstimatesfiltered, effectSizeEstimates);
            assertEquals(effectSizeEstimatesfiltered, effectSizes[filteredVariantIndex++]);
        }
    }
}