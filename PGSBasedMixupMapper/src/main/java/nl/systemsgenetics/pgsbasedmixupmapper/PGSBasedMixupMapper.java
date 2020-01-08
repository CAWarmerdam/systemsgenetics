package nl.systemsgenetics.pgsbasedmixupmapper;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.jet.math.tdouble.DoubleFunctions;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import gnu.trove.map.hash.THashMap;
import nl.systemsgenetics.gwassummarystatistics.VcfGwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.VcfGwasSummaryStatisticsException;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import nl.systemsgenetics.polygenicscorecalculator.SimplePolyGenicScoreCalculator;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Patrick Deelen
 */
public class PGSBasedMixupMapper {

    private static final char CSV_DELIMITER = ',';
    private static final Set<String> TRAITS_TO_MATCH = new HashSet<>(
            Arrays.asList("trait1", "trait2"));
    public static final DecimalFormat LARGE_INT_FORMAT = new DecimalFormat("###,###");
    public static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static final Logger LOGGER = Logger.getLogger(PGSBasedMixupMapper.class);
    private static final String HEADER
            = "  /---------------------------------------\\\n"
            + "  |          PGSBasedMixupMapper          |\n"
            + "  |                                       |\n"
            + "  |  University Medical Center Groningen  |\n"
            + "  \\---------------------------------------/";
    private final RandomAccessGenotypeData genotypeData;
    private Map<String, List<VcfGwasSummaryStatistics>> gwasSummaryStatisticsMap;
    private final Map<String, String> genotypeSampleToPhenotypeSampleCoupling;
    private final DoubleMatrixDataset<String, String> phenotypeMatrix;
    private final DoubleMatrixDataset<String, String> zScoreMatrix;
    private List<String> genotypeSampleIdentifiers;
    private List<String> phenotypeSampleIdentifiers;

    public PGSBasedMixupMapper(RandomAccessGenotypeData genotypeData,
							   ArrayList<Sample> phenotypeSamples,
							   Map<String, String> genotypeSampleToPhenotypeSampleCoupling, Map<String, List<VcfGwasSummaryStatistics>> gwasSummaryStatisticsMap) throws PGSBasedMixupMapperException {
        this.genotypeData = genotypeData;
        this.gwasSummaryStatisticsMap = gwasSummaryStatisticsMap;
        this.genotypeSampleToPhenotypeSampleCoupling = genotypeSampleToPhenotypeSampleCoupling;
        this.setGenotypeSampleIdentifiers();
        this.setPhenotypeSampleIdentifiers(phenotypeSamples);

        // Check if all phenotype samples are unique.
        assertUniquePhenotypeSamples(genotypeSampleToPhenotypeSampleCoupling);

        // Check if there are more than 0 phenotype samples.
        if (phenotypeSamples.size() < 1) {
            throw new IllegalArgumentException("The number of samples cannot be zero (0)");
        }

        // Scale phenotypes
        try {
            phenotypeMatrix = phenotypesToDoubleMatrixDataset(phenotypeSamples);
        } catch (Exception e) {
            throw new PGSBasedMixupMapperException("An error occured while processing phenotype data", e);
        }

        zScoreMatrix = getzScoreMatrix();

        for (String phenotype : phenotypeMatrix.getColObjects()) {
            // Calculate the Z scores for every phenotype
            DoubleMatrixDataset<String, String> phenotypeSpecificZscoreMatrix = calculateZscoreMatrix(phenotype);
            // Sum the Z score.
            zScoreMatrix.getMatrix()
                    .assign(phenotypeSpecificZscoreMatrix.getMatrix(), DoubleFunctions.plus);
        }
        zScoreMatrix.getMatrix().assign(DoubleFunctions.div(phenotypeMatrix.getColObjects().size()));

        // Get the samples
        Map<String, String> bestMatchingPhenotypeSamplePerGenotype = new HashMap<>();
        for (String genotypeSample : genotypeSampleIdentifiers) {
            determineBestMatchingPhenotypeSample(genotypeSample, bestMatchingPhenotypeSamplePerGenotype);
        }

        resolveMixUps(bestMatchingPhenotypeSamplePerGenotype);
    }

    private DoubleMatrixDataset<String, String> calculateZscoreMatrix(String phenotype) {
        DoubleMatrixDataset<String, String> zScoreMatrixOfPolyGenicScoreDeviations =
                new DoubleMatrixDataset<>(genotypeSampleIdentifiers, phenotypeSampleIdentifiers);

        // Initialize residuals of trait prediction
        double[] residuals = new double[genotypeSampleIdentifiers.size()];
        // Initialize polygenic scores
        double[] polyGenicScores = new double[genotypeSampleIdentifiers.size()];

        // Calculate the polygenic score for every genotype sample
        for (int i = 0; i < genotypeSampleIdentifiers.size(); i++) {
            String genotypeSample = genotypeSampleIdentifiers.get(i);
            // Calculate polygenic score
            polyGenicScores[i] = calculatePolyGenicScore(phenotype, genotypeSample);
            // Obtain the actual phenotype score
            double actualPhenotypeScore = phenotypeMatrix.getElement(
                    genotypeSampleToPhenotypeSampleCoupling.get(genotypeSample), phenotype);
            // Calculate the residual
            residuals[i] = actualPhenotypeScore - polyGenicScores[i];
        }

        // Calculate the standard deviation and the mean of the residuals.
        double variance = StatUtils.populationVariance(residuals);
        double sd = Math.sqrt(variance);
        double mean = StatUtils.mean(residuals);

        // Calculate Z-score for every sample combination
        for (int i = 0; i < genotypeSampleIdentifiers.size(); i++) {
            String genotypeSample = genotypeSampleIdentifiers.get(i);
            for (String phenotypeSample : phenotypeSampleIdentifiers) {
                double actualTraitScores = phenotypeMatrix.getElement(phenotypeSample, phenotype);
                double polyGenicScoreDeviation = Math.abs(actualTraitScores - polyGenicScores[i]);

                // Divide every standard deviation with the standard deviations within the population of (...)
                zScoreMatrixOfPolyGenicScoreDeviations.setElement(
                        genotypeSample, phenotypeSample, (polyGenicScoreDeviation - mean) / sd);
            }
        }
        return zScoreMatrixOfPolyGenicScoreDeviations;
    }

    private void setGenotypeSampleIdentifiers() {
        genotypeSampleIdentifiers = new ArrayList<>(this.genotypeSampleToPhenotypeSampleCoupling.keySet());
        if (!genotypeData.getSamples().stream()
                .map(Sample::getId).collect(Collectors.toList())
                .containsAll(genotypeSampleIdentifiers)) {
            throw new IllegalArgumentException("genotype data does not contain all samples from the link file");
        }
    }

    private void setPhenotypeSampleIdentifiers(ArrayList<Sample> phenotypeSamples) {
        phenotypeSampleIdentifiers = new ArrayList<>(this.genotypeSampleToPhenotypeSampleCoupling.values());
        if (!genotypeData.getSamples().stream()
                .map(Sample::getId).collect(Collectors.toList())
                .containsAll(genotypeSampleIdentifiers)) {
            throw new IllegalArgumentException("phenotype data does not contain all samples from the link file");
        }
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

    private double calculatePolyGenicScore(String phenotype, String genotypeSampleIdentifier) {
        List<VcfGwasSummaryStatistics> summaryStatistics = gwasSummaryStatisticsMap.get(phenotype);

        int[] windowSize = new int[]{50000};

        File outputFolder = new File("");
        boolean sumRisks = false;
        double[] pValThres = new double[]{Math.pow(5 * 10, -8), 0.05};
        double rSquare = 1.0d;
        boolean unweighted = false;
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = SimplePolyGenicScoreCalculator
                .summStatsToConvolutedDataStructure(
                        summaryStatistics, pValThres, null,
                        unweighted, LOGGER.isDebugEnabled());

        if (windowSize.length == 1) {
            DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = SimplePolyGenicScoreCalculator.calculate(genotypeData, risks, outputFolder, rSquare, windowSize[0], LOGGER.isDebugEnabled(), pValThres, sumRisks);
//			} else if (windowSize.length == 2) {
//				DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = SimplePolyGenicScoreCalculator.calculateTwoStages(genotypeData, risks, outputFolder, rSquare, windowSize, debugMode, pValThres, sumRisks);
        } else {
            System.out.println("More than two window-sizes is not supported.");
            System.exit(0);
        }

        return 1;
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
                System.out.println("no mix up");
            } else {
                // Check if the originally matched genotype is correctly matched to the trait sample.
                if (traitSample.equals(bestMatchingPhenotypeSamplePerGenotypeSample.get(
                        traitSampleToGenotypeSample.get(traitSample)))) {
                    // Discard the genotype sample currently being assessed.
                    System.out.println("should discard");
                } else {
                    // Treat as a mix-up
                    System.out.println("mix up");
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

    private DoubleMatrixDataset<String, String> getzScoreMatrix() {
        return new DoubleMatrixDataset<>(genotypeSampleIdentifiers, phenotypeSampleIdentifiers);
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
            double[] rank = new NaturalRanking(NaNStrategy.REMOVED, TiesStrategy.MAXIMUM).rank(col.toArray());
            for (int j = 0; j < rank.length; j++) {
                rank[j] = rank[j] / rank.length;
            }
            phenotypeMatrix.viewCol(columnIndex).assign(rank);
        }
    }

    private DoubleMatrixDataset<String, String> phenotypesToDoubleMatrixDataset(ArrayList<Sample> phenotypeSamples) throws Exception {
        List<String> phenotypeLabels = getPhenotypeLabels(phenotypeSamples);
        // Initialize phenotype matrix (samples as rows; phenotypes as columns)
        DoubleMatrixDataset<String, String> phenotypeMatrix = new DoubleMatrixDataset<>(
                phenotypeSampleIdentifiers,
                phenotypeLabels);

        for (int i = 0; i < phenotypeLabels.size(); i++) {
            String key = phenotypeLabels.get(i);
            // Filter to make sure that the same phenotype samples are used.
            double[] values = phenotypeSamples.stream()
                    .filter(sample -> phenotypeSampleIdentifiers.contains(sample.getId()))
                    .map(sample -> Double.valueOf(String.valueOf(sample.getAnnotationValues().get(key))))
                    .mapToDouble(Double::doubleValue)
                    .toArray();
            phenotypeMatrix.viewCol(i).assign(values);
        }
        return phenotypeMatrix;
    }

    private ArrayList<String> getPhenotypeLabels(ArrayList<Sample> phenotypeSamples) {
        Set<String> traitLabels = new HashSet<>(phenotypeSamples.get(0).getAnnotationValues().keySet());
        traitLabels.retainAll(TRAITS_TO_MATCH);
        return new ArrayList<>(traitLabels);
    }

    /**
     * @param args the command line arguments
     * @throws InterruptedException
     */
    public static void main(String[] args) throws InterruptedException {

        System.out.println(HEADER);
        System.out.println();
        System.out.println("          --- Version: " + VERSION + " ---");
        System.out.println();
        System.out.println("More information: http://molgenis.org/systemsgenetics");
        System.out.println();

        Date currentDataTime = new Date();
        String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);
        System.out.println("Current date and time: " + startDateTime);
        System.out.println();

        System.out.flush(); //flush to make sure header is before errors
        Thread.sleep(25); //Allows flush to complete

        PGSBasedMixupMapperOptions options;

        // Print a help whenever the arguments have a length of zero.
        if (args.length == 0) {
            PGSBasedMixupMapperOptions.printHelp();
            return;
        }

        // Parse the arguments list
        options = getPgsBasedMixupMapperOptions(args);

        // Check if the output directory is a prefix for output files. May not be a directory
        if (new File(options.getOutputBasePath()).isDirectory()) {
            System.err.println("Specified output path is a directory. Please include a prefix for the output files.");
            return;
        }

        // Create a logger (set the correct file for output)
        createLogger(startDateTime, options);

        // Print the options
        options.printOptions();

        // Load Genotype and trait data
        RandomAccessGenotypeData genotypeData = getGenotypeData(options);
        ArrayList<Sample> traitSamples = loadTraitSamples(options);

        // Load GWAS summary statistics
        Map<String, String> gwasPhenotypeCoupling = loadGwasSummaryStatisticsPhenotypeCouplings(
                options.getGwasSummaryStatisticsPhenotypeCouplingFile(), CSV_DELIMITER);
        Map<String, List<VcfGwasSummaryStatistics>> gwasSummaryStatisticsMap = loadGwasSummaryStatisticsMap(
                options.getGwasSummaryStatisticsPath(), gwasPhenotypeCoupling);

        Map<String, String> genotypeToPhenotypeSampleCoupling = loadGenotypeToPhenotypeSampleCoupling(
                options.getGenotypeToPhenotypeSampleCouplingFile(), CSV_DELIMITER);

        // Run the mixup mapper stuff
        try {
            PGSBasedMixupMapper pgsBasedMixupMapper = new PGSBasedMixupMapper(
                    genotypeData, traitSamples, genotypeToPhenotypeSampleCoupling, gwasSummaryStatisticsMap);
            pgsBasedMixupMapper.getResults();

        } catch (PGSBasedMixupMapperException e) {
            System.err.println("Error running PGS Based Mixup Mapper: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error accessing input data: " + e.getMessage(), e);
            System.exit(1);
        }
    }

    private static Map<String, String> loadGenotypeToPhenotypeSampleCoupling(String couplingFilePath,
                                                                             char delimiter) {
        Map<String, String> genotypeToPhenotypeSampleCoupling = new LinkedHashMap<>();

        try {
            genotypeToPhenotypeSampleCoupling = loadCouplingMap(couplingFilePath, delimiter);

        } catch (IOException | PGSBasedMixupMapperException e) {
            System.err.println("genotype to phenotype sample coupling file could not be loaded: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("genotype to phenotype sample coupling file could not be loaded: " + e.getMessage(), e);
            System.exit(1);
        }
        return genotypeToPhenotypeSampleCoupling;
    }

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

    private static Map<String, String> loadGwasSummaryStatisticsPhenotypeCouplings(String couplingFilePath,
                                                                                   char delimiter) {
        Map<String, String> gwasSummaryStatisticsPhenotypeCoupling = null;

        try {
            // Create prepare csv parser
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

    private static Map<String, List<VcfGwasSummaryStatistics>> loadGwasSummaryStatisticsMap(
    		String gwasSummaryStatisticsPath, Map<String, String> gwasPhenotypeCoupling) {

        // Initialize map with, for every phenotype, a lists of vcf summary statistics as the value.
    	Map<String, List<VcfGwasSummaryStatistics>> summaryStatisticsMap = new LinkedHashMap<>();

    	// Create the path for the directory where the gwas summary statistic files are stored.
        Path traitSpecificGwasSummaryStatisticsVcfPath =
                Paths.get(gwasSummaryStatisticsPath);
        try {
        	// Try to fill the summary statistics map.
            for (String gwasSummaryStatisticsFilePrefix : gwasPhenotypeCoupling.keySet()) {

            	// Get the phenotype identifier corresponding to the prefix.
				String phenotype = gwasPhenotypeCoupling.get(gwasSummaryStatisticsFilePrefix);
				// Load the vcf gwas summary statistics data with.
				VcfGwasSummaryStatistics summaryStatistics = new VcfGwasSummaryStatistics(
						traitSpecificGwasSummaryStatisticsVcfPath
								.resolve(gwasSummaryStatisticsFilePrefix + ".vcf.gz").toFile(),
						750000, // Represents the cache size in number of variants,
						// value copied from the GeneticRiskScoreCalculator module
						0.4); // Represents the minimum posterior probability to call,
						// 0.4 is generally the default value

				// Check if the map already contains the key,
				// if so add the new summary statistics data in the list.
				if (summaryStatisticsMap.containsKey(phenotype)) {
					summaryStatisticsMap.get(phenotype).add(summaryStatistics);
				} else {
					// If the map does not exist yet, initialize this.
					summaryStatisticsMap.put(phenotype, new ArrayList<>(Collections.singletonList(summaryStatistics)));
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

    private void getResults() {

    }

    private static ArrayList<Sample> loadTraitSamples(PGSBasedMixupMapperOptions options) {
        ArrayList<Sample> phenotypeSamples = null;
        try {
            phenotypeSamples = loadPhenotypeMatrix(options.getInputPhenotypePath(), CSV_DELIMITER);
        } catch (IOException e) {
            System.err.println("Problem running newMixupMapper");
            System.err.println("Error accessing input phenotype data: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error accessing input data: " + e.getMessage(), e);
            System.exit(1);
        }
        System.out.println("phenotypeSamples = " + phenotypeSamples);
        return phenotypeSamples;
    }

    private static RandomAccessGenotypeData getGenotypeData(PGSBasedMixupMapperOptions options) {
        RandomAccessGenotypeData genotypeData = null;

        try {
            genotypeData = loadGenotypes(options);
        } catch (IOException e) {
            System.err.println("Problem running newMixupMapper");
            System.err.println("Error accessing input genotype data: " + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal("Error accessing input data: " + e.getMessage(), e);
            System.exit(1);
        }

        LOGGER.info("Done loading genotype data");
        return genotypeData;
    }

    private static void createLoggerParentDirectory(PGSBasedMixupMapperOptions options) {
        if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
            if (!options.getLogFile().getParentFile().mkdirs()) {
                System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
                System.exit(1);
            }
        }
    }

    private static void createLogger(String startDateTime, PGSBasedMixupMapperOptions options) {
        createLoggerParentDirectory(options);

        try {
            FileAppender logFileAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath());
            Logger.getRootLogger().removeAllAppenders();
            Logger.getRootLogger().addAppender(logFileAppender);

            LOGGER.info("DEPICT" + VERSION);
            LOGGER.info("Current date and time: " + startDateTime);

            if (options.isDebugMode()) {
                Logger.getRootLogger().setLevel(Level.DEBUG);
                if (options.getDebugFolder().mkdir()) {
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

    private static PGSBasedMixupMapperOptions getPgsBasedMixupMapperOptions(String[] args) {
        try {
            return new PGSBasedMixupMapperOptions(args);
        } catch (ParseException ex) {
            System.err.println("Error parsing commandline: " + ex.getMessage());
            PGSBasedMixupMapperOptions.printHelp();
            System.exit(1);
            return null;
        }
    }

    private static RandomAccessGenotypeData loadGenotypes(PGSBasedMixupMapperOptions options) throws IOException {
        final RandomAccessGenotypeData referenceGenotypeData;

        // Check if the forceSeqName option is incorrectly used.
        if (
                options.getForceSeqName() != null &&
                        options.getInputGenotypeType() != RandomAccessGenotypeDataReaderFormats.SHAPEIT2 &&
                        options.getInputGenotypeType() != RandomAccessGenotypeDataReaderFormats.GEN
        ) {
            LOGGER.fatal("Error cannot force sequence name of: " + options.getInputGenotypeType().getName());
            System.err.println("Error cannot force sequence name of: " + options.getInputGenotypeType().getName());
            System.exit(1);
        }

        final SampleFilter sampleFilter;
        if (options.getGenotypeSamplesFile() != null) {
            sampleFilter = readSampleFile(options.getGenotypeSamplesFile());
        } else {
            sampleFilter = null;
        }

        VariantFilter variantFilter = null;

        if (options.getMafFilter() != 0) {
            VariantFilter mafFilter = new VariantFilterMaf(options.getMafFilter());
            if (variantFilter == null) {
                variantFilter = mafFilter;
            } else {
                variantFilter = new VariantCombinedFilter(variantFilter, mafFilter);
            }
        }

        referenceGenotypeData = options.getInputGenotypeType().createFilteredGenotypeData(
                options.getInputGenotypePath(),
                10000,
                variantFilter,
                sampleFilter,
                options.getForceSeqName(),
                0.34f);

        return referenceGenotypeData;
    }

    private static SampleIdIncludeFilter readSampleFile(File sampleFile) throws IOException {

        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withCSVParser(parser).withSkipLines(0).build();

        final HashSet<String> samples = new HashSet<>();

        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {

            samples.add(nextLine[0]);

        }

        return new SampleIdIncludeFilter(samples);

    }

    private static ArrayList<Sample> loadPhenotypeMatrix(String inputPhenotypePath, char delimiter) throws IOException {

        final CSVParser parser = new CSVParserBuilder()
                .withSeparator(delimiter)
                .withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(
                new FileReader(inputPhenotypePath)))
                .withCSVParser(parser).build();

        // Get the header containing column names.
        String[] header = reader.readNext();
        int numberOfColumns = header.length;

        LinkedHashSet<String> ids = new LinkedHashSet<>();
        ArrayList<Sample> samples = new ArrayList<>();
        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {

            // The first columns should contain the individual id
            String individual_id = nextLine[0];
            // The individual id should not have been added already
            if (!ids.add(individual_id)) {
                throw new IllegalArgumentException("Duplicate individual id name: " + individual_id);
            }
            // Initialize annotation map
            Map<String, Object> annotationValues = new HashMap<>();

            // Check if the number of columns remains consistent
            if (numberOfColumns != nextLine.length) {
                throw new IllegalArgumentException("Different number of ids");
            }
            // Fill the annotation map
            for (int i = 1; i < header.length; i++) {
                annotationValues.put(header[i], nextLine[i]);
            }
            // Add the sample to the samples list
            samples.add(new Sample(individual_id, null, annotationValues));
        }
        return samples;
    }
}