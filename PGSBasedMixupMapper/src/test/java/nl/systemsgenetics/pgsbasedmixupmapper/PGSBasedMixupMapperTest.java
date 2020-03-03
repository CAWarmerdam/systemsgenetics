package nl.systemsgenetics.pgsbasedmixupmapper;

import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatistics;
import nl.systemsgenetics.polygenicscorecalculator.SimplePolygenicScoreCalculator;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.sampleFilter.SampleFilterableGenotypeDataDecorator;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.File;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

import static nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapper.*;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class PGSBasedMixupMapperTest {

    private File exampleGwasCouplingFile = new File(this.getClass()
            .getResource("/summarystatistics/gwasstats_phen.coupling.txt").toURI());
    private File exampleGenotypeData = new File(this.getClass()
            .getResource("/genotypedata/1000G.gwasstats_filtered-biallelic-403.vcf.gz").toURI());
    private File exampleSampleCouplingFile = new File(this.getClass()
            .getResource("/demoFiles/dummySampleCoupling.txt").toURI());
    private File examplePhenotypeFile = new File(this.getClass()
            .getResource("/phenotypedata/dummyphenotypes.tsv").toURI());
    private File expectedZScores = new File(this.getClass()
            .getResource("/expected/pgs_z-score.txt").toURI());
    private Path gwasSummaryStatisticsPath = new File(this.getClass()
            .getResource("/summarystatistics/").toURI()).toPath();
    private String[] defaultGenomicRangesToExclude = new String[]{"6:25000000-35000000"};
    private List<String> expectedSampleNames = Arrays.asList(
            "NA06986", "NA06985", "NA06989", "NA06994", "NA07000",
            "NA07037", "NA07048", "NA07051", "NA07056", "NA07347");

    public PGSBasedMixupMapperTest() throws URISyntaxException {
    }

    private static final char CSV_DELIMITER = ',';

    @org.testng.annotations.Test
    public void testCalculateZScoreMatrix() throws Exception {

        // Load the genotype to phenotype sample coupling map
        Map<String, String> genotypeToPhenotypeSampleCoupling = loadGenotypeToPhenotypeSampleCoupling(
                exampleSampleCouplingFile.toString(),
                CSV_DELIMITER);

        // Load the gwas to phenotype coupling map
        Map<String, String> gwasPhenotypeCoupling = loadGwasSummaryStatisticsPhenotypeCouplings(
                exampleGwasCouplingFile.toString(),
                CSV_DELIMITER);

        // Load trait data, only including the samples specified in the coupling map.
        DoubleMatrixDataset<String, String> phenotypeData = loadPhenotypeData(
                new HashSet<>(genotypeToPhenotypeSampleCoupling.values()),
                new HashSet<>(gwasPhenotypeCoupling.values()),
                examplePhenotypeFile);

        // Get the filter out the samples from the coupling file that could not be found.
        genotypeToPhenotypeSampleCoupling = genotypeToPhenotypeSampleCoupling.entrySet()
                .stream()
                .filter(map -> phenotypeData.getRowObjects().contains(map.getValue()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        Map<String, String[]> inputGenotypePaths = new HashMap<>();
        inputGenotypePaths.put(null, new String[]{exampleGenotypeData.getPath()
                .replace(".vcf.gz", "")});

        // Load Genotype data, only including the samples specified in the coupling map.
        RandomAccessGenotypeData genotypeData = loadGenotypeData(
                inputGenotypePaths, RandomAccessGenotypeDataReaderFormats.VCF,
                genotypeToPhenotypeSampleCoupling.keySet(), defaultGenomicRangesToExclude,
                0,
                false);

        // Initialize an Simple polygenic score calculator
        SimplePolygenicScoreCalculator polygenicScoreCalculator = new SimplePolygenicScoreCalculator(
                genotypeData,
                new ArrayList<>(Arrays.asList(5000, 1000)),
                new ArrayList<>(Collections.singletonList(0.01)),
                0.2,
                false,
                defaultGenomicRangesToExclude);

        try {
            // Get the gwas summary statistics map
            Map<String, GwasSummaryStatistics> gwasSummaryStatisticsMap = loadFilteredGwasSummaryStatisticsMap(
                    gwasPhenotypeCoupling,
                    genotypeData,
                    gwasSummaryStatisticsPath);

            // Initialize the Mix-up mapper
            PGSBasedMixupMapper pgsBasedMixupMapper = new PGSBasedMixupMapper(
                    genotypeData, phenotypeData.duplicate(), genotypeToPhenotypeSampleCoupling,
                    polygenicScoreCalculator);

            pgsBasedMixupMapper.calculatePolygenicScores(gwasSummaryStatisticsMap);

            pgsBasedMixupMapper.calculateZScoreMatrix();
            DoubleMatrixDataset<String, String> actualPolygenicScores = pgsBasedMixupMapper
                    .getPolygenicScores("HDL cholesterol");

            assertEquals(actualPolygenicScores.getRowObjects(), Collections.singletonList("IEU-a-780_P0.01"));
            assertEquals(actualPolygenicScores.getColObjects(), expectedSampleNames);

            assertEquals(actualPolygenicScores.getRow(0).toArray(),
                    new double[]{-0.657700, 2.15940, 1.28790, 0.0356000, 0.743200,
                            2.39350, 2.37500, -0.685600, 2.25310, 1.67860}, 1e-5);

            DoubleMatrixDataset<String, String> zScoreMatrix = pgsBasedMixupMapper.getZScoreMatrix();

            DoubleMatrixDataset<String, String> expectedZScores = DoubleMatrixDataset.loadDoubleTextData(
                    this.expectedZScores.toString(), '\t');

            assertEquals(
                    Arrays.stream(zScoreMatrix.getMatrixAs2dDoubleArray()).flatMapToDouble(Arrays::stream).toArray(),
                    Arrays.stream(expectedZScores.getMatrixAs2dDoubleArray()).flatMapToDouble(Arrays::stream).toArray(),
                    1e-5);

        } catch (PGSBasedMixupMapperException e) {
            System.err.println("Error running PGS Based Mixup Mapper: " + e.getMessage());
            System.err.println("See log file for stack trace");
            System.exit(1);
        }
    }

    @org.testng.annotations.Test
    public void testCalculateSplitPolygenicScores() {
        // Load the genotype to phenotype sample coupling map
        Map<String, String> genotypeToPhenotypeSampleCoupling = loadGenotypeToPhenotypeSampleCoupling(
                exampleSampleCouplingFile.toString(),
                CSV_DELIMITER);

        // Load the gwas to phenotype coupling map
        Map<String, String> gwasPhenotypeCoupling = loadGwasSummaryStatisticsPhenotypeCouplings(
                exampleGwasCouplingFile.toString(),
                CSV_DELIMITER);

        // Load trait data, only including the samples specified in the coupling map.
        DoubleMatrixDataset<String, String> phenotypeData = loadPhenotypeData(
                new HashSet<>(genotypeToPhenotypeSampleCoupling.values()),
                new HashSet<>(gwasPhenotypeCoupling.values()),
                examplePhenotypeFile);

        // Get the filter out the samples from the coupling file that could not be found.
        genotypeToPhenotypeSampleCoupling = genotypeToPhenotypeSampleCoupling.entrySet()
                .stream()
                .filter(map -> phenotypeData.getRowObjects().contains(map.getValue()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        Map<String, String[]> inputGenotypePaths = new HashMap<>();
        inputGenotypePaths.put(null, new String[]{exampleGenotypeData.getPath()
                .replace(".vcf.gz", "")});

        // Load Genotype data, only including the samples specified in the coupling map.
        RandomAccessGenotypeData genotypeData = loadGenotypeData(
                inputGenotypePaths, RandomAccessGenotypeDataReaderFormats.VCF,
                genotypeToPhenotypeSampleCoupling.keySet(), defaultGenomicRangesToExclude,
                0,
                false);

        // Initialize an Simple polygenic score calculator
        SimplePolygenicScoreCalculator polygenicScoreCalculator = new SimplePolygenicScoreCalculator(
                genotypeData,
                new ArrayList<>(Arrays.asList(5000, 1000)),
                new ArrayList<>(Collections.singletonList(0.01)),
                0.2,
                false,
                defaultGenomicRangesToExclude);

        // Get the gwas summary statistics map
        Map<String, GwasSummaryStatistics> gwasSummaryStatisticsMap = loadFilteredGwasSummaryStatisticsMap(
                gwasPhenotypeCoupling,
                genotypeData,
                gwasSummaryStatisticsPath);

        SampleIdIncludeFilter referenceSampleFilter = new SampleIdIncludeFilter(expectedSampleNames.subList(0, 5));
        List<String> expectedSampleSubset = expectedSampleNames.subList(0, 3);
        SampleIdIncludeFilter responseSampleFilter = new SampleIdIncludeFilter(expectedSampleSubset);

        DoubleMatrixDataset<String, String> actualPolygenicScores = polygenicScoreCalculator.calculate(
                gwasSummaryStatisticsMap.get("HDL cholesterol"), referenceSampleFilter, responseSampleFilter);

        assertEquals(actualPolygenicScores.getRowObjects(), Collections.singletonList("IEU-a-780_P0.01"));
        assertEquals(actualPolygenicScores.getColObjects(), expectedSampleSubset);

        // Calculate polygenic scores without using the subsetting method.
        RandomAccessGenotypeData filteredGenotypeData = new SampleFilterableGenotypeDataDecorator(
                genotypeData, referenceSampleFilter);

        // Initialize an Simple polygenic score calculator
        SimplePolygenicScoreCalculator filteredPolygenicScoreCalculator = new SimplePolygenicScoreCalculator(
                filteredGenotypeData,
                new ArrayList<>(Arrays.asList(5000, 1000)),
                new ArrayList<>(Collections.singletonList(0.01)),
                0.2,
                false,
                defaultGenomicRangesToExclude);

        DoubleMatrixDataset<String, String> expectedPolygenicScores = filteredPolygenicScoreCalculator.calculate(
                gwasSummaryStatisticsMap.get("HDL cholesterol"));

        assertEquals(actualPolygenicScores.getRow(0).toArray(),
                expectedPolygenicScores.viewColSelection(expectedSampleSubset).getRow(0).toArray(), 1e-5);
    }
}
