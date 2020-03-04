package nl.systemsgenetics.pgsbasedmixupmapper.robustnessanalyzer;

import com.opencsv.CSVWriter;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatistics;
import nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapper;
import nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapperException;
import nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapperOptions;
import nl.systemsgenetics.polygenicscorecalculator.SimplePolygenicScoreCalculator;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.RandomAccessGenotypeData;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

import static nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapper.*;

/**
 * Tests and reports the robustness of the PGSBasedMixupMapper
 * @author Robert Warmerdam
 */
public class RobustnessAnalyzer {

    private static final String VERSION = ResourceBundle.getBundle("version").getString("application.version");
    private static final Logger LOGGER = Logger.getLogger(RobustnessAnalyzer.class);
    private static final char CSV_DELIMITER = ',';
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    private Random random = new Random();
    private final Map<String, String> correctCouplingMap;
    private Map<String, String> permutedMap;

    public RobustnessAnalyzer(Map<String, String> correctCouplingMap) {
        this.correctCouplingMap = correctCouplingMap;
    }

    public int permute(int nMixUpsToIntroduce) {
        Map<String, String> permutedSampleCouplingMap = new HashMap<>(correctCouplingMap);
        List<String> keyList = new ArrayList<>(permutedSampleCouplingMap.keySet());


        Collections.shuffle(keyList, random);

        boolean[] hasBeenSwapped = new boolean[nMixUpsToIntroduce];

        System.out.println("Mixing up samples: " + keyList.subList(0, nMixUpsToIntroduce));

        for (int i = nMixUpsToIntroduce - 1; 0 < i; i--) {

            // Check if the current value must be swapped.
            boolean mustSwapCurrent = (!hasBeenSwapped[i] || (i == 1 && !hasBeenSwapped[0]));

            int j = i;

            if (mustSwapCurrent) {
                // sample random integer j so that 0 <= j < i
                j = random.nextInt(i);
            } else {
                // sample random integer j so that 0 <= j <= i
                j = random.nextInt(i + 1);
            }

            if (j != i) {
                // swapping j with i guarantees that the every key is swapped:
                // key i can never be swapped with itself, and can neither be swapped in the future.
                String a = keyList.get(i);
                String b = keyList.get(j);
                permutedSampleCouplingMap.put(a, permutedSampleCouplingMap.put(b, permutedSampleCouplingMap.get(a)));
                Collections.swap(keyList, i, j);
                hasBeenSwapped[i] = true;
                hasBeenSwapped[j] = true;
            }
        }
        this.permutedMap = permutedSampleCouplingMap;

        int numberOfShuffledItems = 0;
        for (boolean b : hasBeenSwapped) {
            numberOfShuffledItems += b ? 1 : 0;
        }
        return numberOfShuffledItems;
    }

    /**
     * Method that writes a table with the genotype identifiers in the first column,
     * the correct phenotype coupling in the second column,
     * the used phenotype coupling in the third column and
     * the output coupling in the last column.
     *
     * @param resolvedSampleCouplings A map with genotype identifiers as keys and
     *                                mapped phenotype identifiers as values.
     * @param outputFile The file to write the results to.
     * @throws IOException If an I/O error occurred.
     */
    public void save(Map<String, String> resolvedSampleCouplings, File outputFile) throws IOException {

        // Initialize the CSV writer for outputting the differences.
        CSVWriter writer = new CSVWriter(new FileWriter(
                String.format("%s_inducedMixUps.tsv", outputFile)),
                CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        writer.writeNext(new String[]{"gen", "correct", "permuted", "corrected"});

        for (String key : correctCouplingMap.keySet()) {
            String correctValue = correctCouplingMap.get(key);
            String usedValue = permutedMap.get(key);
            String potentiallyResolvedValue = resolvedSampleCouplings.get(key);

            writer.writeNext(new String[]{key, correctValue, usedValue, potentiallyResolvedValue});
        }
        writer.close();
    }

    /**
     * Method that analyses results and checks whether induced sample mix-ups are
     * identified and resolved, if sample mix-ups are incorrectly identified.
     *
     * @param resolvedSampleCouplings A map with genotype identifiers as keys and
     *                                mapped phenotype identifiers as values.
     * @return The
     */
    public Map<String, Integer> analyzeMapperResults(Map<String, String> resolvedSampleCouplings) {
        // Check what the new sample assignments look like
        // Measure:
        // True positives (TP); the number of mixed-up samples that are correctly identified.
        // False positives (FP); the number samples that were falsely deemed a mix-up
        // True negatives (TN); the number of samples correctly identified as not being mixed-up
        Map<String, Integer> statistics = new HashMap<>();

        statistics.put("TP", 0);
        statistics.put("FP", 0);
        statistics.put("TN", 0);
        statistics.put("FN", 0);

        for (String key : correctCouplingMap.keySet()) {
            String correctValue = correctCouplingMap.get(key);
            String usedValue = permutedMap.get(key);
            String potentiallyResolvedValue = resolvedSampleCouplings.get(key);
//            System.out.println("------------------------------");
//
//            System.out.println("correctValue = " + correctValue);
//            System.out.println("usedValue = " + usedValue);
//            System.out.println("potentiallyResolvedValue = " + potentiallyResolvedValue);

            // A value has been induced if the used value is not equal to the correct value
            boolean isMixedUp = !correctValue.equals(usedValue);
            boolean hasChanged = !usedValue.equals(potentiallyResolvedValue);
            boolean isMixUpMissed = isMixedUp && !hasChanged;
            boolean isMixUpIncorrectlyIdentified = !isMixedUp && hasChanged;
            boolean isMixUpCorrectlyIdentified = isMixedUp && hasChanged;
            boolean isCorrectlyResolved = isMixUpCorrectlyIdentified && correctValue.equals(potentiallyResolvedValue);

            if (isCorrectlyResolved) statistics.put("TP", statistics.get("TP") + 1);
            if (isMixUpIncorrectlyIdentified) statistics.put("FP", statistics.get("FP") + 1);
            if (!isMixedUp && !hasChanged) statistics.put("TN", statistics.get("TN") + 1);
            if (isMixUpMissed) statistics.put("FN", statistics.get("FN") + 1);

//            System.out.println(statistics);
        }
        return statistics;
    }

    public static void main(String[] args) {
        // Parse the options in the PGSBasedMixupMapper, with additionally
        // options with which to get an idea of the robustness of the mixup mapper
        // strategy. This might include:
        // - induced mix-up percentages

        // Get the current date and time.
        Date currentDataTime = new Date();
        String startDateTime = DATE_TIME_FORMAT.format(currentDataTime);

        // Initialize an instance for command line options.
        RobustnessAnalyzerOptions options;

        // Parse the arguments list
        options = getOptions(args);

        // Create a logger (set the correct file for output)
        createLogger(startDateTime, options);

        // Print the options
        options.printOptions();

        // Read the data like in the PGSBasedMixupMapper

        // Load the genotype to phenotype sample coupling map
        Map<String, String> correctSampleCouplingMap = loadGenotypeToPhenotypeSampleCoupling(
                options.getGenotypeToPhenotypeSampleCouplingFile(), CSV_DELIMITER);

        // Load the gwas to phenotype coupling map
        Map<String, String> gwasPhenotypeCoupling = loadGwasSummaryStatisticsPhenotypeCouplings(
                options.getGwasSummaryStatisticsPhenotypeCouplingFile(), CSV_DELIMITER);

        // Load trait data, only including the samples specified in the coupling map.
        DoubleMatrixDataset<String, String> phenotypeData = loadPhenotypeData(
                new HashSet<>(correctSampleCouplingMap.values()),
                new HashSet<>(gwasPhenotypeCoupling.values()), options.getInputPhenotypePath());

        // Get the filter out the samples from the coupling file that could not be found.
        correctSampleCouplingMap = correctSampleCouplingMap.entrySet()
                .stream()
                .filter(map -> phenotypeData.getRowObjects().contains(map.getValue()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        // Load Genotype data, only including the samples specified in the coupling map.
        RandomAccessGenotypeData genotypeData = loadGenotypeData(
                options.getInputGenotypePaths(),
                options.getInputGenotypeType(),
                correctSampleCouplingMap.keySet(),
                options.getGenomicRangesToExclude(),
                options.getMinorAlleleFrequencyThreshold(),
                options.shouldForceSeqName());

        // Initialize an Simple polygenic score calculator
        SimplePolygenicScoreCalculator polygenicScoreCalculator = new SimplePolygenicScoreCalculator(
                genotypeData,
                options.getWindowSize(),
                options.getpValueThresholds(),
                options.getrSquared(),
                false,
                options.getGenomicRangesToExclude());

        // Get the mix-up percentages.
        List<Double> mixUpPercentages = options.getMixUpPercentages();

        // Define a matrix to use for filling in the different statistics from the robustness analyzer.
        ArrayList<String> columnNames = new ArrayList<>(
                Arrays.asList("N", "Introduced", "TP", "FP", "FN", "TN"));
        DoubleMatrixDataset<Double, String> inducedMixUpResults = new DoubleMatrixDataset<>(
                mixUpPercentages,
                columnNames);
        
        // Check the difference between the g
        Map<String, Integer> statistics = new HashMap<>();

        // Get a new robustness analyzer.
        RobustnessAnalyzer robustnessAnalyzer = new RobustnessAnalyzer(correctSampleCouplingMap);

        for (Double mixUpPercentage : mixUpPercentages) {
            // For every number / percentage of induced sample mix-ups
            // Shuffle the number in the coupling map.
            int nMixUpsToIntroduce = (int) Math.round(
                    ((float) correctSampleCouplingMap.size()) / 100 * mixUpPercentage);

            // Permute the samples.
            int numberOfShuffledSamples = robustnessAnalyzer.permute(nMixUpsToIntroduce);

            // Log the number of samples that are mixed up and the total number of samples
            LOGGER.info(String.format("Total number of samples: %d", correctSampleCouplingMap.size()));
            inducedMixUpResults.setElement(mixUpPercentage, "N", correctSampleCouplingMap.size());
            LOGGER.info(String.format("Number of permuted samples: %d", numberOfShuffledSamples));
            inducedMixUpResults.setElement(mixUpPercentage, "Introduced", numberOfShuffledSamples);

            // Calculate the actual percentage of mixed up samples.
            float percentageOfMixedUpSamples =
                    (numberOfShuffledSamples / (float) correctSampleCouplingMap.size()) * 100;

            LOGGER.info(String.format("Actual percentage of permuted samples: %f", percentageOfMixedUpSamples));

            // Get the permuted map.
            Map<String, String> permutedSampleCouplingMap = robustnessAnalyzer.getPermutedMap();

            Path path = Paths.get(String.format(Locale.ROOT, "%s_robustness%.1f", options.getOutputBasePath(),
                    mixUpPercentage));
            boolean mkdir = path.toFile().mkdir();
            assert mkdir;
            File mixUpPercentageOutputPath = path.resolve("output").toFile();
            
            try {
                // Get the gwas summary statistics map
                Map<String, GwasSummaryStatistics> gwasSummaryStatisticsMap = getGwasSummaryStatisticsMap(options,
                        permutedSampleCouplingMap, gwasPhenotypeCoupling, phenotypeData, genotypeData);

                // Initialize the Mix-up mapper
                PGSBasedMixupMapper pgsBasedMixupMapper = new PGSBasedMixupMapper(
                        genotypeData, phenotypeData.duplicate(), permutedSampleCouplingMap, polygenicScoreCalculator);

                pgsBasedMixupMapper.calculatePolygenicScores(gwasSummaryStatisticsMap);

                // Run the mix-up mapper
                pgsBasedMixupMapper.calculateZScoreMatrix();

                // Report results
                Map<String, String> resolvedSampleCouplingMap = pgsBasedMixupMapper.getSampleAssignments();

                reportResults(mixUpPercentageOutputPath, pgsBasedMixupMapper, phenotypeData, resolvedSampleCouplingMap);

                statistics = robustnessAnalyzer.analyzeMapperResults(resolvedSampleCouplingMap);
                robustnessAnalyzer.save(
                        resolvedSampleCouplingMap,
                        mixUpPercentageOutputPath);

                System.out.println("statistics = " + statistics);

            } catch (PGSBasedMixupMapperException e) {
                String errorMessage = String.format(
                        "Error running PGS based Mix-up mapper with %d / %d (%.2f%%) mixed-up samples: ",
                        nMixUpsToIntroduce, permutedSampleCouplingMap.size(), percentageOfMixedUpSamples);

                System.err.println(errorMessage + e.getMessage());
                System.err.println("See log file for stack trace");
                LOGGER.fatal(errorMessage + e.getMessage(), e);
                System.exit(1);
            } catch (IOException e) {
                String errorMessage = String.format(
                        "Error saving robustness results for %d / %d (%.2f%%) mixed-up samples: ",
                        nMixUpsToIntroduce, permutedSampleCouplingMap.size(), percentageOfMixedUpSamples);

                System.err.println(errorMessage + e.getMessage());
                System.err.println("See log file for stack trace");
                LOGGER.warn(errorMessage + e.getMessage(), e);
            }

            // Write the true positives, false positives, false negatives, true negatives.
            for (String key : statistics.keySet()) {
                inducedMixUpResults.setElement(mixUpPercentage, key, statistics.get(key));
            }
        }
        
        try {
            inducedMixUpResults.save(options.getOutputBasePath() + "_robustnessMatrix");
        } catch (IOException e) {
            String errorMessage = "Error saving robustness results: ";

            System.err.println(errorMessage + e.getMessage());
            System.err.println("See log file for stack trace");
            LOGGER.fatal(errorMessage + e.getMessage(), e);
            System.exit(1);
        }
    }

    /**
     * Method that initializes a Logger instance with a log file.
     *
     * @param startDateTime The date and time at which the program is started.
     * @param options The specified command line options.
     */
    static void createLogger(String startDateTime, RobustnessAnalyzerOptions options) {
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
     * Method that processes the command line options that are specified by the user.
     *
     * @param args The command line arguments that are specified.
     * @return A PGSBasedMixupMapperOptions instance that holds the processed command line arguments.
     */
    static RobustnessAnalyzerOptions getOptions(String[] args) {
        try {
            return new RobustnessAnalyzerOptions(args);
        } catch (ParseException ex) {
            System.err.println("Error parsing commandline: " + ex.getMessage());
            PGSBasedMixupMapperOptions.printHelp();
            System.exit(1);
            return null;
        }
    }

    public Map<String, String> getPermutedMap() {
        return permutedMap;
    }
}
