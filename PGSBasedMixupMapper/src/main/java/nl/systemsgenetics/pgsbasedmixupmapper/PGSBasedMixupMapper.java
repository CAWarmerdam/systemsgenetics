package nl.systemsgenetics.pgsbasedmixupmapper;

import cern.colt.function.tdouble.DoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.gwassummarystatistics.MultiStudyGwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.VariantFilterableGwasSummaryStatisticsDecorator;
import nl.systemsgenetics.gwassummarystatistics.VcfGwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.VcfGwasSummaryStatisticsException;
import nl.systemsgenetics.polygenicscorecalculator.*;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.*;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Robert Warmerdam
 */
public class PGSBasedMixupMapper {

    private static final char CSV_DELIMITER = ',';
    private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static final Logger LOGGER = Logger.getLogger(PGSBasedMixupMapper.class);
    private static final String HEADER
            = "  /---------------------------------------\\\n"
            + "  |          PGSBasedMixupMapper          |\n"
            + "  |                                       |\n"
            + "  |  University Medical Center Groningen  |\n"
            + "  \\---------------------------------------/";
    private final RandomAccessGenotypeData genotypeData;
    private Map<String, List<MultiStudyGwasSummaryStatistics>> gwasSummaryStatisticsMap;
    private final Map<String, String> genotypeSampleToPhenotypeSampleCoupling;
    private PolyGenicScoreCalculator polyGenicScoreCalculator;
    private final DoubleMatrixDataset<String, String> phenotypeMatrix;
    private final DoubleMatrixDataset<String, String> zScoreMatrix;
    private final Map<String, DoubleMatrixDataset<String, String>> zScoresMap = new HashMap<>();
    private final Map<String, DoubleMatrixDataset<String, String>> polyGenicScoresMap = new HashMap<>();
    private List<String> genotypeSampleIdentifiers;
    private List<String> phenotypeSampleIdentifiers;
    private Map<String, Integer> minimumAbsoluteResidualsIndices = new HashMap<>();
    private Set<Pair<String, String>> missingPhenotypes = new HashSet<>();
    private final DenseDoubleMatrix1D presenceArray;

    /**
     * Constructor method for resolving mix-ups based on polygenic scores.
     *
     * @param genotypeData The genotype data to resolve sample mix-ups in.
     * @param phenotypeMatrix The individual's phenotype data to resolve sample mix-ups in.
     * @param genotypeSampleToPhenotypeSampleCoupling Map with genotype sample identifiers as the keys
     *                                                with the original phenotype sample identifiers as values.
     * @param gwasSummaryStatisticsMap A map with GwasSummaryStatistics data per trait or phenotype.
     * @param polyGenicScoreCalculator An object with which polygenic scores are calculated.
     * @throws PGSBasedMixupMapperException if
     */
    public PGSBasedMixupMapper(RandomAccessGenotypeData genotypeData,
                               DoubleMatrixDataset<String, String> phenotypeMatrix,
                               Map<String, String> genotypeSampleToPhenotypeSampleCoupling,
                               Map<String, List<MultiStudyGwasSummaryStatistics>> gwasSummaryStatisticsMap,
                               PolyGenicScoreCalculator polyGenicScoreCalculator)
            throws PGSBasedMixupMapperException, PolyGenicScoreCalculatorException {
        this.genotypeData = genotypeData;
        this.gwasSummaryStatisticsMap = gwasSummaryStatisticsMap;
        this.genotypeSampleToPhenotypeSampleCoupling = genotypeSampleToPhenotypeSampleCoupling;
        this.polyGenicScoreCalculator = polyGenicScoreCalculator;
        this.phenotypeMatrix = phenotypeMatrix;

        // Get the genotype sample identifiers, and check if all samples from the coupling file are present.
        this.genotypeSampleIdentifiers = this.getGenotypeSampleIdentifiers();
        // Get the phenotype sample identifiers, and check if all samples from the coupling file are present.
        this.phenotypeSampleIdentifiers = this.getPhenotypeSampleIdentifiers();

        // Check if all phenotype samples are unique.
        assertUniquePhenotypeSamples(genotypeSampleToPhenotypeSampleCoupling);

        // Check if there are more than 0 phenotype samples.
        if (phenotypeMatrix.getRowObjects().size() < 1) {
            throw new IllegalArgumentException("The number of samples cannot be zero (0)");
        }

        // Scale phenotypes
        rankPhenotypes();
        zScoreMatrix = initializeGenotypePhenotypeMatrix();
        presenceArray = new DenseDoubleMatrix1D(phenotypeSampleIdentifiers.size());

        for (String phenotype : this.phenotypeMatrix.getColObjects()) {
            System.out.println(String.format("Calculating PGSs and corresponding Z-scores for trait '%s'", phenotype));
            // Calculate the Z scores for every phenotype
            try {
                DoubleMatrixDataset<String, String> phenotypeSpecificZScoreMatrix = calculateZScoreMatrix(phenotype);
                zScoresMap.put(phenotype, phenotypeSpecificZScoreMatrix);
                // Sum the Z score.
                zScoreMatrix.getMatrix()
                        .assign(phenotypeSpecificZScoreMatrix.getMatrix(), DoubleFunctions.plus);

                presenceArray.assign(getPhenotypePresence(phenotype), DoubleFunctions.plus);

            } catch (PGSBasedMixupMapperException e) {
                System.err.println("Skipping trait, Z-scores could not be calculated: " + e.getMessage());
                System.err.println("See log file for stack trace");
                LOGGER.warn("Skipping trait, Z-scores could not be calculated: " + e.getMessage(), e);
            }
        }

        // Create a 2D matrix from the 1D matrix representing a counter for the number of traits to divide the total
        // Z-scores with for every phenotype sample.
        DoubleMatrix2D otherMatrix = DoubleFactory2D.dense.repeat(
                DoubleFactory2D.dense.make(presenceArray.toArray(), 1),
                zScoreMatrix.rows(), 1);

        // Divide the Z-scores by the number of traits used for getting this Z-score
        zScoreMatrix.getMatrix().assign(otherMatrix, DoubleFunctions.div);

        // Check if there is a phenotype sample for which no presence was detected
        if (presenceArray.getMinLocation()[0] == 0) {
            LOGGER.warn("For at least one phenotype sample, no traits were available for comparison");
            System.err.println("Warning, for at least one phenotype sample, no traits were available for comparison");
        }

        // Check if there are Z-scores with
        if (Arrays.asList(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY)
                .contains(zScoreMatrix.getMatrix().getMaxLocation()[0])) {
            LOGGER.warn("The Z-score for at least one phenotype sample is infinite");
            System.err.println("Warning, the Z-score for at least one phenotype sample is infinite");
            zScoreMatrix.getMatrix().assign(v -> v == Double.NEGATIVE_INFINITY ? Double.POSITIVE_INFINITY : v);
        }
        if (zScoreMatrix.getMatrix().equals(zScoreMatrix.getElementQuick(0, 0))) {
            throw new PGSBasedMixupMapperException(String.format(
                    "All Z-scores are equal (%f). Cannot determine best matching samples.",
                    zScoreMatrix.getElementQuick(0, 0)));
        }

        // Normalize Z-scores to have mean at 0 along rows and cols
        // Normalize Z-scores to have an SD of 1
        zScoreMatrix.normalizeRows();
        zScoreMatrix.normalizeColumns();

        // Get the samples
        Map<String, String> bestMatchingPhenotypeSamplePerGenotype = new HashMap<>();
        for (String genotypeSample : genotypeSampleIdentifiers) {
            determineBestMatchingPhenotypeSample(genotypeSample, bestMatchingPhenotypeSamplePerGenotype);
        }

        resolveMixUps(bestMatchingPhenotypeSamplePerGenotype);
    }

    private DoubleMatrix1D getPhenotypePresence(String phenotype) {
        double[] phenotypePresence = new double[phenotypeSampleIdentifiers.size()];

        // Calculate the polygenic score for every genotype sample
        for (int i = 0; i < phenotypeSampleIdentifiers.size(); i++) {
            String phenotypeSampleIdentifier = phenotypeSampleIdentifiers.get(i);

            // Check if the Phenotype Score is missing
            // If this is missing, skip this phenotype score
            if (missingPhenotypes.contains(new ImmutablePair<>(phenotypeSampleIdentifier, phenotype))) {
                phenotypePresence[i] = 0;
            } else {
                phenotypePresence[i] = 1;
            }
        }
        return new DenseDoubleMatrix1D(phenotypePresence);
    }

    /**
     * Calculates standard (Z) scores for every combination of phenotype and genotype samples for the given phenotype.
     * Per combination, this indicates how many standard deviations (SDs) the actual phenotype minus the PGS,
     * deviates from the calculated residuals calculated between the original link of genotype and phenotype sample.
     *
     * @param phenotype The phenotype to calculate the standard (Z) score for.
     * @return The
     */
    private DoubleMatrixDataset<String, String> calculateZscoreMatrixNew(String phenotype) throws PolyGenicScoreCalculatorException, PGSBasedMixupMapperException {
        DoubleMatrixDataset<String, String> zScoreMatrixOfPolyGenicScoreDeviations =
                new DoubleMatrixDataset<>(genotypeSampleIdentifiers, phenotypeSampleIdentifiers);

        // Initialize polygenic scores
        DoubleMatrixDataset<String, Double> polyGenicScores = polyGenicScoreCalculator.calculate(genotypeData,
                gwasSummaryStatisticsMap.get(phenotype).get(0));

        // Initialize a matrix of residuals
        DoubleMatrixDataset<String, Double> residualsMatrix = new DoubleMatrixDataset<>(
                polyGenicScores.getRowObjects(), polyGenicScores.getColObjects());

        // Calculate the polygenic score for every genotype sample
        for (int i = 0; i < genotypeSampleIdentifiers.size(); i++) {
            String genotypeSample = genotypeSampleIdentifiers.get(i);
            // Calculate polygenic score
            // Obtain the actual phenotype score
            double actualPhenotypeScore = phenotypeMatrix.getElement(
                    genotypeSampleToPhenotypeSampleCoupling.get(genotypeSample), phenotype);
            // Calculate the residual
            residualsMatrix.viewRow(i).assign(actualPhenotypeScore);
        }

        // Index corresponding to the index of the residuals array with the lowest sum,
        // indicating that the p-value is the best
        int minimumResidualsIndex = -1;
        // Initialize value to compare residuals sum with.
        double absoluteMinimumResiduals = Double.MAX_VALUE;
        // Calculate for every row the residuals and get the row with the lowest sum of residuals
        for (int colIndex = 0; colIndex < residualsMatrix.getColObjects().size(); colIndex++) {
            double[] values = rankArray(polyGenicScores.viewCol(colIndex).toArray());
            polyGenicScores.viewCol(colIndex).assign(
                    values);

            // Calculate the residuals
            residualsMatrix.viewCol(colIndex).assign(polyGenicScores.viewCol(colIndex),
                    DoubleFunctions.minus);

            // Determine if the current row is the lowest
            double absoluteSumOfResiduals = Math.abs(residualsMatrix.viewCol(colIndex)
                    .copy().assign(DoubleFunctions.abs).zSum());
            if (absoluteSumOfResiduals < absoluteMinimumResiduals) {
                // If so, update the minimum
                minimumResidualsIndex = colIndex;
                absoluteMinimumResiduals = absoluteSumOfResiduals;
            }
        }

        // Get the residuals corresponding to the best P-value
        double[] residuals = residualsMatrix.viewCol(minimumResidualsIndex).toArray();

        // Calculate the standard deviation and the mean of the residuals.
        double variance = StatUtils.populationVariance(residuals);
        double sd = Math.sqrt(variance);
        double mean = StatUtils.mean(residuals);

        if (variance == 0) {
            throw new PGSBasedMixupMapperException("No variance in residuals");
        }

        // Calculate Z-score for every sample combination
        for (int sampleIndex = 0; sampleIndex < genotypeSampleIdentifiers.size(); sampleIndex++) {
            String genotypeSample = genotypeSampleIdentifiers.get(sampleIndex);
            for (String phenotypeSample : phenotypeSampleIdentifiers) {
                // For the current combination, calculate the residual.
                double actualTraitScores = phenotypeMatrix.getElement(phenotypeSample, phenotype);
                // Calculate the difference between the actual trait score and the PGS
                double polyGenicScoreDeviation = actualTraitScores - polyGenicScores
                        .getElementQuick(sampleIndex, minimumResidualsIndex);

                // Divide every standard deviation with the standard deviations within the population of (...)
                zScoreMatrixOfPolyGenicScoreDeviations.setElement(
                        genotypeSample, phenotypeSample, Math.abs(polyGenicScoreDeviation - mean) / sd);
            }
        }
        // Return the Z scores of polygenic score deviations.
        return zScoreMatrixOfPolyGenicScoreDeviations;
    }

    /**
     * Calculates standard (Z) scores for every combination of phenotype and genotype samples for the given phenotype.
     * Per combination, this indicates how many standard deviations (SDs) the actual phenotype minus the PGS,
     * deviates from the calculated residuals calculated between the original link of genotype and phenotype sample.
     *
     * @param phenotype The phenotype to calculate the standard (Z) score for.
     * @return The Z score matrix.
     */
    private DoubleMatrixDataset<String, String> calculateZScoreMatrix(String phenotype) throws PGSBasedMixupMapperException {
        DoubleMatrixDataset<String, String> zScoreMatrixOfPolyGenicScoreDeviations =
                initializeGenotypePhenotypeMatrix();

        // Initialize polygenic scores
        DoubleMatrixDataset<String, String> polyGenicScores = calculatePolyGenicScores(phenotype);
        polyGenicScoresMap.put(phenotype, polyGenicScores.duplicate());

        // Initialize a matrix of residuals
        DoubleMatrixDataset<String, String> residualsMatrix = new DoubleMatrixDataset<>(
                polyGenicScores.getRowObjects(), polyGenicScores.getColObjects());

        List<Integer> presentValuesIndices = new ArrayList<>();

        // Calculate the polygenic score for every genotype sample
        for (int i = 0; i < polyGenicScores.getColObjects().size(); i++) {
            String genotypeSample = polyGenicScores.getColObjects().get(i);
            // Calculate polygenic score
            // Obtain the actual phenotype score
            String phenotypeSampleIdentifier = genotypeSampleToPhenotypeSampleCoupling.get(genotypeSample);
            double actualPhenotypeScore = phenotypeMatrix.getElement(
                    phenotypeSampleIdentifier, phenotype);

            // Check if the Phenotype Score is missing
            // If this is missing, skip this phenotype score
            if (!missingPhenotypes.contains(new ImmutablePair<>(phenotypeSampleIdentifier, phenotype))) {
                presentValuesIndices.add(i);
            }

            // Calculate the residual
            residualsMatrix.viewCol(i).assign(actualPhenotypeScore);
        }

        // Get an array with an index for every sample for which the phenotype value is missing
        int[] indicesOfPresentPhenotypeValues = presentValuesIndices.stream().mapToInt(Integer::intValue).toArray();

        // Index corresponding to the index of the residuals array with the lowest sum,
        // indicating that the p-value is the best
        int minimumResidualsIndex = -1;
        // Initialize value to compare residuals sum with.
        double absoluteMinimumResiduals = Double.MAX_VALUE;
        // Calculate for every row the residuals and get the row with the lowest sum of residuals
        for (int rowIndex = 0; rowIndex < residualsMatrix.getRowObjects().size(); rowIndex++) {
            double[] values = rankArray(polyGenicScores.viewRow(rowIndex).toArray());
            polyGenicScores.viewRow(rowIndex).assign(
                    values);

            // Calculate the residuals
            residualsMatrix.viewRow(rowIndex).assign(polyGenicScores.viewRow(rowIndex),
                    DoubleFunctions.minus);

            // Determine if the current row is the lowest
            double absoluteSumOfResiduals = Math.abs(residualsMatrix.viewRow(rowIndex)
                    .viewSelection(indicesOfPresentPhenotypeValues)
                    .copy() // Copying since assigning the absolute values will
                    // otherwise be reflected in the residualsMatrix
                    .assign(DoubleFunctions.abs).zSum());
            if (absoluteSumOfResiduals < absoluteMinimumResiduals) {
                // If so, update the minimum
                minimumResidualsIndex = rowIndex;
                absoluteMinimumResiduals = absoluteSumOfResiduals;
            }
        }

        minimumAbsoluteResidualsIndices.put(phenotype, minimumResidualsIndex);

        // Get the residuals corresponding to the best P-value
        double[] residuals = residualsMatrix.viewRow(minimumResidualsIndex)
                .viewSelection(indicesOfPresentPhenotypeValues).toArray();

        // Calculate the standard deviation and the mean of the residuals.
        double variance = StatUtils.populationVariance(residuals);
        double sd = Math.sqrt(variance);
        double mean = StatUtils.mean(residuals);

        if (variance == 0) {
            throw new PGSBasedMixupMapperException("No variance in residuals");
        }

        // Calculate Z-score for every sample combination
        for (int sampleIndex = 0; sampleIndex < polyGenicScores.getColObjects().size(); sampleIndex++) {
            String genotypeSample = polyGenicScores.getColObjects().get(sampleIndex);
//            System.out.printf("%s (gen)/%s (phen): %s | %s%n",
//                    genotypeSample, genotypeSampleToPhenotypeSampleCoupling.get(genotypeSample),
//                    polyGenicScores.getElementQuick(minimumResidualsIndex, sampleIndex),
//                    phenotypeMatrix.getElement(genotypeSampleToPhenotypeSampleCoupling.get(genotypeSample), phenotype));
            for (String phenotypeSample : phenotypeSampleIdentifiers) {
                if (missingPhenotypes.contains(new ImmutablePair<>(phenotypeSample, phenotype))) {
                    // Set the deviation to zero
                    zScoreMatrixOfPolyGenicScoreDeviations.setElement(
                            genotypeSample,
                            phenotypeSample, 0);
                } else {

                    // For the current combination, calculate the residual.
                    double actualTraitScores = phenotypeMatrix.getElement(phenotypeSample, phenotype);
                    // Calculate the difference between the actual trait score and the PGS
                    double polygenicScoreDeviation = actualTraitScores - polyGenicScores
                            .getElementQuick(minimumResidualsIndex, sampleIndex);

                    // Divide every standard deviation with the standard deviations within the population of (...)
                    zScoreMatrixOfPolyGenicScoreDeviations.setElement(
                            genotypeSample, phenotypeSample, Math.abs(polygenicScoreDeviation - mean) / sd);
                }
            }
        }
        // Return the Z scores of polygenic score deviations.
        return zScoreMatrixOfPolyGenicScoreDeviations;
    }

    private List<String> getGenotypeSampleIdentifiers() {
        List<String> genotypeSampleIdentifiers = new ArrayList<>(this.genotypeSampleToPhenotypeSampleCoupling.keySet());
        if (!genotypeData.getSamples().stream()
                .map(Sample::getId).collect(Collectors.toList())
                .containsAll(genotypeSampleIdentifiers)) {
            throw new IllegalArgumentException("genotype data does not contain all samples from the coupling file");
        }
        return genotypeSampleIdentifiers;
    }

    private List<String> getPhenotypeSampleIdentifiers() {
        List<String> phenotypeSampleIdentifiers = new ArrayList<>(this.genotypeSampleToPhenotypeSampleCoupling.values());
        List<String> OrderedPhenotypeSampleIdentifiers = phenotypeMatrix.getRowObjects();
        if (!OrderedPhenotypeSampleIdentifiers
                .containsAll(phenotypeSampleIdentifiers)) {
            throw new IllegalArgumentException("phenotype data does not contain all samples from the coupling file");
        }
        return OrderedPhenotypeSampleIdentifiers;
    }

    private void assertUniquePhenotypeSamples(Map<String, String> genotypeSampleToTraitSample) {
        if (genotypeSampleToTraitSample.values().stream().distinct().count() != genotypeSampleToTraitSample.size()) {
            throw new IllegalArgumentException(String.format(
                    "Original sample coupling does not exclusively " +
                            "contain distinct phenotype samples." +
                            "%n%d duplicate labels found",
                    genotypeSampleToTraitSample.size() - genotypeSampleToTraitSample.values().stream().distinct().count()));
        }
    }

    private DoubleMatrixDataset<String, String> calculatePolyGenicScores(String phenotype) {
        List<MultiStudyGwasSummaryStatistics> summaryStatistics = gwasSummaryStatisticsMap.get(phenotype);


        List<Integer> windowSizeList = polyGenicScoreCalculator.getLdHandler().getWindowSize();
        int[] windowSize = windowSizeList.stream().mapToInt(i->i).toArray();
        boolean sumRisks = false;
        double[] pValThres = polyGenicScoreCalculator.getpValueThresholds()
                .stream().mapToDouble(Double::doubleValue).toArray();
        double rSquare = polyGenicScoreCalculator.getLdHandler().getRSquaredThreshold();
        boolean unweighted = false;
        String[] genomicRangesToExclude = polyGenicScoreCalculator.getGenomicRangesToExclude();
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = SimplePolyGenicScoreCalculator
                .summStatsToConvolutedDataStructure(genotypeData,
                        summaryStatistics, pValThres, genomicRangesToExclude,
                        unweighted, LOGGER.isDebugEnabled());

        DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = null;

        if (windowSize.length == 1) {
            geneticRiskScoreMatrix = SimplePolyGenicScoreCalculator.calculate(
                    genotypeData, risks, rSquare, windowSize[0], pValThres, sumRisks);
        } else if (windowSize.length == 2) {
            geneticRiskScoreMatrix = SimplePolyGenicScoreCalculator.calculateTwoStages(
                    genotypeData, risks, rSquare, windowSize, pValThres, sumRisks);
        } else {
            System.out.println("More than two window-sizes is not supported.");
            System.exit(0);
        }
        return geneticRiskScoreMatrix;
    }

    private void resolveMixUps(Map<String, String> bestMatchingPhenotypeSamplePerGenotypeSample) {
        Map<String, String> traitSampleToGenotypeSample =
                genotypeSampleToPhenotypeSampleCoupling.entrySet()
                        .stream()
                        .collect(Collectors.toMap(Map.Entry::getValue, Map.Entry::getKey));
        for (Map.Entry<String, String> newLink : bestMatchingPhenotypeSamplePerGenotypeSample.entrySet()) {
            // Based on the Z-score, is there enough evidence for a mix-up?
            String genotypeSample = newLink.getKey();
            String traitSample = newLink.getValue();
            if (!isBestMatchSufficient(genotypeSample, traitSample)) {
                bestMatchingPhenotypeSamplePerGenotypeSample.remove(genotypeSample);
            }
            // Check if the new link matches the original link
            if (genotypeSampleToPhenotypeSampleCoupling.get(genotypeSample).equals(traitSample)) {
                // No mix-up
                System.out.printf("%s Not mixed up%n", genotypeSample);
            } else {
                // Check if the originally matched genotype is correctly matched to the trait sample.
                if (traitSample.equals(bestMatchingPhenotypeSamplePerGenotypeSample.get(
                        traitSampleToGenotypeSample.get(traitSample)))) {
                    // Discard the genotype sample currently being assessed.
//                    System.out.println("should discard");
                } else {
                    // Treat as a mix-up
                    System.out.printf("%s Mixed up with %s%n", genotypeSample, traitSampleToGenotypeSample.get(traitSample));
                }
            }
            if (LOGGER.isDebugEnabled()) {
                LOGGER.debug(String.format("Assessed genotype sample '%s', matched with '%s' (z = %f)",
                        genotypeSample, traitSample, phenotypeMatrix.getElement(genotypeSample, traitSample)));
            }
        }
    }

    private boolean isBestMatchSufficient(String key, String value) {
        return true;
    }

    private DoubleMatrixDataset<String, String> initializeGenotypePhenotypeMatrix() {
        return new DoubleMatrixDataset<>(Arrays.asList(genotypeData.getSampleNames()), phenotypeSampleIdentifiers);
    }

    private void determineBestMatchingPhenotypeSample(String genotypeSample, Map<String, String> bestMatchingPhenotypeSamplePerGenotype) {
        DoubleMatrix1D scores = zScoreMatrix.getCol(genotypeSample);
        int minimumValueIndex = (int) scores.getMinLocation()[1];
        bestMatchingPhenotypeSamplePerGenotype.put(genotypeSample,
                phenotypeMatrix.getRowObjects().get(minimumValueIndex));
    }

    private void rankPhenotypes() {
        for (int columnIndex = 0; columnIndex < phenotypeMatrix.columns(); columnIndex++) {
            DoubleMatrix1D col = phenotypeMatrix.getCol(columnIndex);
            double[] rank = rankArray(col.toArray());
            phenotypeMatrix.viewCol(columnIndex).assign(rank);
        }
    }

    private double[] rankArray(double[] array) {
        double[] rank = new NaturalRanking(NaNStrategy.REMOVED, TiesStrategy.MINIMUM).rank(array);
        for (int j = 0; j < rank.length; j++) {
            rank[j] = rank[j] / rank.length;
        }
        return rank;
    }

    public Map<String, DoubleMatrixDataset<String, String>> getZScoresMap() {
        return zScoresMap;
    }

    public DoubleMatrixDataset<String, String> getZScoreMatrix() {
        return zScoreMatrix;
    }

    /**
     * @param args the command line arguments
     * @throws InterruptedException
     */
    public static void main(String[] args) throws InterruptedException {

        // Get the current date and time.
        Date currentDataTime = new Date();
        String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);

        // Print a header
        System.out.printf("%s%n%n" +
                "Version: %s%n%n" +
                "More information: http://molgenis.org/systemsgenetics%nCurrent date and time: %s%n",
                HEADER, VERSION, startDateTime);

        System.out.flush(); //flush to make sure header is before errors
        Thread.sleep(25); //Allows flush to complete

        // Initialize an instance for command line options.
        PGSBasedMixupMapperOptions options;

        // Print a help whenever the arguments have a length of zero.
        if (args.length == 0) {
            PGSBasedMixupMapperOptions.printHelp();
            return;
        }

        // Parse the arguments list
        options = getPgsBasedMixupMapperOptions(args);

        // Create a logger (set the correct file for output)
        createLogger(startDateTime, options);

        // Print the options
        options.printOptions();

        // Load the genotype to phenotype sample coupling map
        Map<String, String> genotypeToPhenotypeSampleCoupling = loadGenotypeToPhenotypeSampleCoupling(
                options.getGenotypeToPhenotypeSampleCouplingFile(), CSV_DELIMITER);

        // Load the gwas to phenotype coupling map
        Map<String, String> gwasPhenotypeCoupling = loadGwasSummaryStatisticsPhenotypeCouplings(
                options.getGwasSummaryStatisticsPhenotypeCouplingFile(), CSV_DELIMITER);

        // Load trait data, only including the samples specified in the coupling map.
        DoubleMatrixDataset<String, String> phenotypeData = loadPhenotypeData(options,
                new HashSet<>(genotypeToPhenotypeSampleCoupling.values()),
                new HashSet<>(gwasPhenotypeCoupling.values()));

        genotypeToPhenotypeSampleCoupling = genotypeToPhenotypeSampleCoupling.entrySet()
                .stream()
                .filter(map -> phenotypeData.getRowObjects().contains(map.getValue()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        // Load Genotype data, only including the samples specified in the coupling map.
        RandomAccessGenotypeData genotypeData = getGenotypeData(options,
                genotypeToPhenotypeSampleCoupling.keySet());

        // Get a variant filter that includes only variants als present in the given genotype data.
        VariantFilter variantFilter = getVariantFilter(genotypeData);

        // Load GWAS summary statistics applying the variant filter.
        Map<String, List<MultiStudyGwasSummaryStatistics>> gwasSummaryStatisticsMap = loadGwasSummaryStatisticsMap(
                options.getGwasSummaryStatisticsPath(), gwasPhenotypeCoupling, variantFilter);

        // Get the P-value thresholds to use in PGS calculation
        List<Double> pValueThresholds = options.getpValueThresholds();

        try {
            // Initialize an LD handler for clumping the effect entries into the top, unlinked, clumps
            LDHandler ldHandler = new Clumper(
                    genotypeData, // The genotype data to calculate an LD matrix in.
                    options.getWindowSize(), // Get the window size in number of base pairs
                    options.getrSquared(),  // Use R2 as a threshold
                    pValueThresholds.get(pValueThresholds.size() - 1)); // Use the last (least stringiest) p-value,
                    // effect entries with even higher p-values (lower log(p))should not for clumps

            // Initialize a PGS calculator.
            PolyGenicScoreCalculator pgsCalculator = new PolyGenicScoreCalculator(ShrinkageStrategy.THRESHOLDING,
                    ldHandler,
                    options.getGenomicRangesToExclude());
            pgsCalculator.setPValueThresholds(pValueThresholds);

//            DoubleMatrixDataset<String, Double> hdl_cholesterol = pgsCalculator.calculate(genotypeData,
//                    gwasSummaryStatisticsMap.get("HDL cholesterol").get(0));
//
//            hdl_cholesterol.save(String.format("%s_PGSs_test_%s.tsv",
//                    options.getOutputBasePath(), "hdl_cholesterol"));

            // Initialize the Mix-up mapper
            PGSBasedMixupMapper pgsBasedMixupMapper = new PGSBasedMixupMapper(
                    genotypeData, phenotypeData.duplicate(), genotypeToPhenotypeSampleCoupling,
                    gwasSummaryStatisticsMap, pgsCalculator);
            // Report results
            reportResults(options, pgsBasedMixupMapper, phenotypeData);

        } catch (PGSBasedMixupMapperException | PolyGenicScoreCalculatorException e) {
            System.err.println("Error running PGS Based Mixup Mapper: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error running PGS Based Mixup Mapper: " + e.getMessage(), e);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("Error saving output from PGS Based Mixup Mapper: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error saving output from PGS Based Mixup Mapper: " + e.getMessage(), e);
            System.exit(1);
        } catch (LDHandlerException | LdCalculatorException e) {
            System.err.println("Error calculating LD matrix: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error calculating LD matrix: " + e.getMessage(), e);
            System.exit(1);
        }
    }

    private static void reportResults(PGSBasedMixupMapperOptions options, PGSBasedMixupMapper pgsBasedMixupMapper, DoubleMatrixDataset<String, String> phenotypeData) throws IOException {
        for (Map.Entry<String, DoubleMatrixDataset<String, String>> zScores : pgsBasedMixupMapper.getZScoresMap().entrySet()) {
            zScores.getValue().save(String.format("%s_zScoreMatrix_%s.tsv",
                    options.getOutputBasePath(), zScores.getKey().replace(' ', '-')));
            DoubleMatrixDataset<String, String> polyGenicScores = pgsBasedMixupMapper
                    .polyGenicScoresMap.get(zScores.getKey());
            DoubleMatrix1D firstPolygenicScores = pgsBasedMixupMapper.getBestPolygenicScores(zScores.getKey());
            DoubleMatrix1D actualTraits = phenotypeData.viewRowSelection(polyGenicScores.getColObjects()).viewCol(zScores.getKey());
            DoubleMatrixDataset<String, String> comparison = new DoubleMatrixDataset<>(
                    polyGenicScores.getColObjects(),
                    new ArrayList<>(Arrays.asList("actual", "PGS")));
            comparison.viewCol(0).assign(actualTraits);
            comparison.viewCol(1).assign(firstPolygenicScores);
            comparison.save(String.format("%s_TraitsVsPGSs_%s.tsv",
                    options.getOutputBasePath(), zScores.getKey().replace(' ', '-')));
            polyGenicScores.save(String.format("%s_PGSs_%s.tsv",
                    options.getOutputBasePath(), zScores.getKey().replace(' ', '-')));
        }

        pgsBasedMixupMapper.getZScoreMatrix().save(
                String.format("%s_zScoreMatrix_%s.tsv",
                        options.getOutputBasePath(), "overall"));
    }

    private DoubleMatrix1D getBestPolygenicScores(String key) {
        return polyGenicScoresMap.get(key).viewRow(minimumAbsoluteResidualsIndices.get(key));
    }

    /**
     * Method that builds a variant filter from the genotype data for including only the variants present in
     * the given genotype data.
     *
     * @param genotypeData The genotype data that specifies which genetic variants the filter should include.
     * @return A variant filter that includes the variants present in the genotype data.
     */
    static VariantFilter getVariantFilter(RandomAccessGenotypeData genotypeData) {
        Set<String> variantIdentifiers = new HashSet<>();
        for (GeneticVariant variant : genotypeData) {
            variantIdentifiers.add(variant.getPrimaryVariantId());
        }
        return new VariantIdIncludeFilter(variantIdentifiers);
    }

    /**
     * Method loading the map that couples genotype samples to phenotype samples
     * The input file should contain two columns with the first one representing the genotype sample,
     * and the second representing the corresponding phenotype sample.
     *
     * @param couplingFilePath The path to the file that contains this information.
     * @param delimiter The delimiter that separate columns.
     * @return A map with GWAS summary statistic file prefix as keys and the corresponding trait as values.
     */
    static Map<String, String> loadGenotypeToPhenotypeSampleCoupling(String couplingFilePath,
                                                                     char delimiter) {
        Map<String, String> genotypeToPhenotypeSampleCoupling = new LinkedHashMap<>();

        try {
            // Load the coupling map.
            genotypeToPhenotypeSampleCoupling = loadCouplingMap(couplingFilePath, delimiter);

            // Log the number of samples that are loaded
            LOGGER.info(String.format("Loaded a genotype to phenotype sample coupling map for %d samples",
                    genotypeToPhenotypeSampleCoupling.size()));
        } catch (IOException | PGSBasedMixupMapperException e) {
            System.err.println("genotype to phenotype sample coupling file could not be loaded: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("genotype to phenotype sample coupling file could not be loaded: " + e.getMessage(), e);
            System.exit(1);
        }
        return genotypeToPhenotypeSampleCoupling;
    }

    /**
     * Method that loads a file of delimiter separated values, with two columns. The first column corresponds
     * to the genotype samples to include. The second column represents the phenotype samples linked to the
     * genotype samples from the first column
     *
     * @param couplingFilePath The path to the file that holds the values.
     * @param delimiter The character with which values are separated.
     * @return A linked hash map that holds the genotype samples as keys and the phenotype samples as values.
     * @throws IOException If an I/O exception occurs while reading the input file.
     * @throws PGSBasedMixupMapperException If the file was not formatted as described; either the file contains
     * a different number of columns or there are duplicate values within one of the columns.
     */
    private static Map<String, String> loadCouplingMap(String couplingFilePath, char delimiter) throws IOException, PGSBasedMixupMapperException {
        Map<String, String> loadedCouplingMap = new LinkedHashMap<>();
        // Create prepare csv parser
        final CSVParser parser = new CSVParserBuilder()
                .withSeparator(delimiter)
                .withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(
                new FileReader(couplingFilePath)))
                .withCSVParser(parser).build();

        String[] nextLine;
        int lineIndex = 0;
        while ((nextLine = reader.readNext()) != null) {
            // We want strictly 2 values every line:
            // one for the prefix of the gwas vcf file and
            // one for the phenotype label that this represents.
            if (nextLine.length != 2) {
                throw new PGSBasedMixupMapperException(String.format(
                        "Encountered %d values on line %d while 2 where expected",
                        nextLine.length, lineIndex));
            }
            String gwasSummaryStatisticsFilePrefix = nextLine[0];    // The first column (index zero (0)) on this
            // line respresents the prefix for a vcf(.gz(.tbi)) file(s)
            // Throw an exception if the prefix is already used.
            if (loadedCouplingMap.containsKey(gwasSummaryStatisticsFilePrefix)) {
                throw new PGSBasedMixupMapperException(String.format(
                        "Encountered a duplicate value '%s' on line %d, this is not supported.",
                        gwasSummaryStatisticsFilePrefix, lineIndex));
            }
            loadedCouplingMap.put(
                    gwasSummaryStatisticsFilePrefix,
                    nextLine[1]);    // The second column (index one (1)) on this line represents the phenotype
            // label that corresponds to the gwas summary stats in the vcf file.

            // Count the line
            lineIndex++;
        }
        return loadedCouplingMap;
    }

    /**
     * Method loading the map that couples GWAS summary statistics files to traits.
     * The input file should contain two columns with the first one representing the GWAS summary statistic file,
     * and the second representing the corresponding trait.
     *
     * @param couplingFilePath The path to the file that contains this information.
     * @param delimiter The delimiter that separate columns.
     * @return A map with GWAS summary statistic file prefix as keys and the corresponding trait as values.
     */
    static Map<String, String> loadGwasSummaryStatisticsPhenotypeCouplings(String couplingFilePath,
                                                                           char delimiter) {
        // Initialize coupling map
        Map<String, String> gwasSummaryStatisticsPhenotypeCoupling = null;

        try {
            // Load the map
            gwasSummaryStatisticsPhenotypeCoupling = loadCouplingMap(couplingFilePath, delimiter);

        } catch (IOException | PGSBasedMixupMapperException e) {
            System.err.println("GWAS summary statistics phenotype coupling file could not be loaded: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("GWAS summary statistics phenotype coupling file could not be loaded: " + e.getMessage(), e);
            System.exit(1);
        }

        // Return the stuff...
        return gwasSummaryStatisticsPhenotypeCoupling;
    }

    /**
     * Method that loads gwas summary statistics data from GWAS VCF files.
     * These files correspond to the prefixes from the given gwas phenotype coupling map, + .vcf.gz.
     * An additional .tbi index file per VCF file is also required to load the file.
     *
     * @param gwasSummaryStatisticsPath The path to the directory where GWAS VCF files can be found.
     * @param gwasPhenotypeCoupling A map with GWAS VCF file prefixes as keys and corresponding trait identifiers as
     *                              values.
     * @param variantFilter A filter to apply to the GWAS VCF file. When writing this documentation this is used to
     *                      exclude all variants not present in the genotype data that is to be assigned to phenotype
     *                      data.
     * @return A map with for every phenotype (key), a list of gwas summary statistics.
     */
    static Map<String, List<MultiStudyGwasSummaryStatistics>> loadGwasSummaryStatisticsMap(
            String gwasSummaryStatisticsPath, Map<String, String> gwasPhenotypeCoupling, VariantFilter variantFilter) {

        // Initialize map with, for every phenotype, a lists of vcf summary statistics as the value.
        Map<String, List<MultiStudyGwasSummaryStatistics>> summaryStatisticsMap = new LinkedHashMap<>();

        // Create the path for the directory where the gwas summary statistic files are stored.
        Path traitSpecificGwasSummaryStatisticsVcfPath =
                Paths.get(gwasSummaryStatisticsPath);
        try {
            // Try to fill the summary statistics map.
            for (String gwasSummaryStatisticsFilePrefix : gwasPhenotypeCoupling.keySet()) {

                // Get the phenotype identifier corresponding to the prefix.
                String phenotype = gwasPhenotypeCoupling.get(gwasSummaryStatisticsFilePrefix);
                // Load the vcf gwas summary statistics data with.
                File bzipVcfFile = traitSpecificGwasSummaryStatisticsVcfPath
                        .resolve(gwasSummaryStatisticsFilePrefix + ".vcf.gz").toFile();

                MultiStudyGwasSummaryStatistics summaryStatistics = new VcfGwasSummaryStatistics(
                        bzipVcfFile,
                        50000, // Represents the cache size in number of variants,
                        // value copied from the GeneticRiskScoreCalculator module
                        1); // Represents the minimum posterior probability to call,
                        // 0.4 is generally the default value

                // Report loaded status and the variant count
                int variantCount = summaryStatistics.getVariantIdMap().size();
                LOGGER.info(String.format("Loaded GWAS summary statistics from '%s'",
                        bzipVcfFile.toString()));
                LOGGER.info(String.format("%d variants present prior to filtering",
                        variantCount));

                // A VCF gwas file can contain more than one study.
                // We currently do not know how to deal with this, so only the first study is used
                if (summaryStatistics.getSamples().size() > 1) {
                    LOGGER.warn(String.format("Encountered %d studies (VCF samples) while only 1 is expected",
                            summaryStatistics.getSamples().size()));
                    LOGGER.warn("Only the first study will be used.");
                }

                // Remove the variants that are not according to the variant filter.
                // This should correspond to including only the variants that are in the genotype data as well.
                MultiStudyGwasSummaryStatistics filteredSummaryStatistics =
                        new VariantFilterableGwasSummaryStatisticsDecorator(summaryStatistics, variantFilter);
                int filteredVariantCount = filteredSummaryStatistics.getVariantIdMap().size();
                LOGGER.info(String.format("Removing %d variants not present in the genotype data, keeping %d",
                        variantCount - filteredVariantCount,
                        filteredVariantCount));

                // Check if the map already contains the key,
                // if so add the new summary statistics data in the list.
                if (summaryStatisticsMap.containsKey(phenotype)) {
                    summaryStatisticsMap.get(phenotype).add(filteredSummaryStatistics);
                } else {
                    // If the map does not exist yet, initialize this.
                    summaryStatisticsMap.put(phenotype, new ArrayList<>(Collections.singletonList(filteredSummaryStatistics)));
                }
            }

        } catch (IOException | VcfGwasSummaryStatisticsException e) {
            System.err.println("GWAS summary statistics data could not be loaded: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("GWAS summary statistics data could not be loaded: " + e.getMessage(), e);
            System.exit(1);
        }
        return summaryStatisticsMap;
    }

    /**
     * Load the phenotype data.
     *
     * @param options The user specified PGS options.
     * @param phenotypeSampleIdentifiersToInclude A set of strings that contains all the sample identifiers to
     *                                            load phenotypes for.
     * @return A List of samples with the phenotypes annotated within these samples.
     */
    static DoubleMatrixDataset<String, String> loadPhenotypeData(PGSBasedMixupMapperOptions options,
                                                                 Set<String> phenotypeSampleIdentifiersToInclude,
                                                                 Set<String> traitsToInclude) {
        DoubleMatrixDataset<String, String> phenotypeMatrix = null;
        try {
            phenotypeMatrix = loadPhenotypeMatrix(options.getInputPhenotypePath(),
                    phenotypeSampleIdentifiersToInclude, traitsToInclude,'\t');
            phenotypeMatrix.save(options.getOutputBasePath() + "prior_phen_mat.tsv");
            LOGGER.info(String.format("Loaded phenotype data for %d samples",
                    phenotypeMatrix.getRowObjects().size()));
        } catch (IOException e) {
            System.err.println("Error accessing input phenotype data: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error accessing input phenotype data: " + e.getMessage(), e);
            System.exit(1);
        }
        return phenotypeMatrix;
    }

    /**
     * Method for loading genotype data given the PGS based mixup mapper options.
     *
     * @param options The options from which the genotype data path and type can be obtained.
     * @param sampleIdentifiersToInclude The sample identifiers to use for
     *                                   filtering the samples in the genotype data
     * @return Random access genotype data.
     */
    static RandomAccessGenotypeData getGenotypeData(PGSBasedMixupMapperOptions options,
                                                    Set<String> sampleIdentifiersToInclude) {
        // Initialize the genotype data
        RandomAccessGenotypeData randomAccessGenotypeData = null;
        // Get a new variant filter. This is used to only keep biallelec variants.
        VariantFilter variantFilter = new VariantCombinedFilter(
                new VariantFilterBiAllelic(),
                new VariantFilterAmbigousSnp());

        // If a MAF has been set, extend the current variant filter
        if (options.getMafFilter() != 0) {
            VariantFilter mafFilter = new VariantFilterMaf(options.getMafFilter());
            variantFilter = new VariantCombinedFilter(variantFilter, mafFilter);
        }

        // Filter the genotype data on the given set of sample identifiers.
        SampleFilter sampleFilter = new SampleIdIncludeFilter(sampleIdentifiersToInclude);

//        String[][] inputGenotypePaths = options.getInputGenotypePaths();
        String[] inputGenotypePathSet = options.getInputGenotypePath();

        try {
//            Set<RandomAccessGenotypeData> genotypeDataSet = new HashSet<>();
//
//            for (String[] inputGenotypePathSet : inputGenotypePaths) {
//                RandomAccessGenotypeData genotypeData = options.getInputGenotypeType().createFilteredGenotypeData(
//                        inputGenotypePathSet,
//                        750000,
//                        variantFilter,
//                        sampleFilter,
//                        options.getForceSeqName(),
//                        0.5f);
//
//                LOGGER.info(String.format("Loaded genotype data (%s) from '%s'.",
//                        options.getInputGenotypeType().getName(),
//                        Arrays.toString(inputGenotypePathSet)));
//
//                int originalVariantCount = genotypeData.getVariantIdMap().size();
//                LOGGER.info(String.format("%d samples and %d variants present after filtering",
//                        genotypeData.getSamples().size(),
//                        originalVariantCount));
//
//                genotypeDataSet.add(genotypeData);
//            }
//            randomAccessGenotypeData = new MultiPartGenotypeData(genotypeDataSet);
            randomAccessGenotypeData = options.getInputGenotypeType().createFilteredGenotypeData(
                    inputGenotypePathSet,
                    50000,
                    variantFilter,
                    sampleFilter,
                    options.getForceSeqName(),
                    0.5f);

            LOGGER.info(String.format("Loaded genotype data (%s) from '%s'.",
                    options.getInputGenotypeType().getName(),
                    Arrays.toString(inputGenotypePathSet)));

            int originalVariantCount = randomAccessGenotypeData.getVariantIdMap().size();
            LOGGER.info(String.format("%d samples and %d variants present after filtering",
                    randomAccessGenotypeData.getSamples().size(),
                    originalVariantCount));


        } catch (IOException e) {
            System.err.println("Error accessing input genotype data: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error accessing input genotype data: " + e.getMessage(), e);
            System.exit(1);
        }
        return randomAccessGenotypeData;
    }

    /**
     * Method that initializes a Logger instance with a log file.
     *
     * @param startDateTime The date and time at which the program is started.
     * @param options The specified command line options.
     */
    static void createLogger(String startDateTime, PGSBasedMixupMapperOptions options) {
        // Try to make the parent directory of the log file if it is not present.
        if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
            if (!options.getLogFile().getParentFile().mkdirs()) {
                LOGGER.fatal("Failed to create output folder: " + options.getLogFile().getParent());
                System.exit(1);
            }
        }

        try {
            // Apply a new log file.
            FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath());
            Logger.getRootLogger().addAppender(logFileAppender);

            // Log some first info.
            LOGGER.info("PGSBasedMixupMapper" + VERSION);
            LOGGER.info("Current date and time: " + startDateTime);

            if (options.isDebugMode()) {
                Logger.getRootLogger().setLevel(Level.DEBUG);
                if (!options.getDebugFolder().exists() && !options.getDebugFolder().mkdir()) {
                    throw new IOException("Could not create debug logger.");
                }
            } else {
                Logger.getRootLogger().setLevel(Level.INFO);
            }

        } catch (IOException e) {
            System.err.println("Failed to create logger: " + e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Load a CSV of phenotype values as a list of samples.
     *
     * @param inputPhenotypePath The path that points towards the file with the phenotype data.
     * @param phenotypeSampleIdentifiersToInclude A set of strings that contains all the sample identifiers to
     *                                            load phenotypes for.
     * @param delimiter The delimiter with which columns are separated.
     * @return A List of samples with the phenotypes annotated within these samples.
     * @throws IOException If an I/O error occurs while reading the phenotype data.
     */
    private static DoubleMatrixDataset<String, String> loadPhenotypeMatrix(File inputPhenotypePath,
                                                         Set<String> phenotypeSampleIdentifiersToInclude,
                                                         Set<String> traitsToInclude,
                                                         char delimiter) throws IOException {

        // Create the CSV reader object with which the phenotypes can be read through
        final CSVParser parser = new CSVParserBuilder()
                .withSeparator(delimiter)
                .withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(
                new FileReader(inputPhenotypePath)))
                .withCSVParser(parser).build();

        // Create a phenotype matrix with the phenotype sample identifiers as
        DoubleMatrixDataset<String, String> phenotypeMatrix = new DoubleMatrixDataset<>(
                phenotypeSampleIdentifiersToInclude,
                traitsToInclude);

        // Get the header containing column names.
        String[] header = reader.readNext();
        int numberOfColumns = header.length;

        // Create a list from the set of traits to include to be able to iterate through these traits
        // by using an indices
        List<String> traitsToIncludeList = new ArrayList<>(traitsToInclude);
        // Create a list from the first row (header) in the file to look up which column id corresponds to
        // which trait using .indexOf()
        List<String> headerAsList = Arrays.asList(header);

        // Get the index for every trait to include
        int[] traitIndices = new int[traitsToInclude.size()];
        for (int i = 0; i < traitsToInclude.size(); i++) {
            traitIndices[i] = headerAsList
                    .indexOf(traitsToIncludeList.get(i));
        }

        LinkedHashSet<String> ids = new LinkedHashSet<>();
        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {

            // The first columns should contain the individual id
            String individual_id = nextLine[0];

            // If the individual id is in the set indicating which samples to include,
            // include it.
            if (!phenotypeSampleIdentifiersToInclude.contains(individual_id)) {
                // Otherwise, continue on the next iteration with a new next line.
                continue;
            }

            // Check if the number of columns remains consistent
            if (numberOfColumns != nextLine.length) {
                throw new IllegalArgumentException("Different number of ids");
            }

            // Fill the phenotype matrix
            try {
                for (int traitIndex : traitIndices) {
                    if (individual_id.equals("LLDeep_0001")) {
                        System.out.println(header[traitIndex]);
                        System.out.println(nextLine[traitIndex]);
                    }
                    String value = nextLine[traitIndex];
                    phenotypeMatrix.setElement(individual_id, header[traitIndex], Double.parseDouble(value));
                }
                // The individual id should not have been added already
                if (!ids.add(individual_id)) {
                    throw new IllegalArgumentException("Duplicate individual id name: " + individual_id);
                }
            } catch (NumberFormatException e) {
                System.out.println("Encountered error: " + e.getMessage());
                System.out.println("Skipping: " + individual_id);
            }
        }
//        System.out.println(phenotypeMatrix.getRowObjects());
//        System.out.println(phenotypeMatrix.getCol(0));
        return phenotypeMatrix.viewRowSelection(ids);
    }

    /**
     * Method that processes the command line options that are specified by the user.
     *
     * @param args The command line arguments that are specified.
     * @return A PGSBasedMixupMapperOptions instance that holds the processed command line arguments.
     */
    static PGSBasedMixupMapperOptions getPgsBasedMixupMapperOptions(String[] args) {
        try {
            return new PGSBasedMixupMapperOptions(args);
        } catch (ParseException ex) {
            System.err.println("Error parsing commandline: " + ex.getMessage());
            PGSBasedMixupMapperOptions.printHelp();
            System.exit(1);
            return null;
        }
    }
}