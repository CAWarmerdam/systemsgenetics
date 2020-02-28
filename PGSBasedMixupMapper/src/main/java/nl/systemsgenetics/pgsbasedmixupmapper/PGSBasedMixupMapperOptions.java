/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.pgsbasedmixupmapper;

import edu.emory.mathcs.utils.ConcurrencyUtils;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author patri
 */
public class PGSBasedMixupMapperOptions {

    private static final Options OPTIONS;
    private static final double DEFAULT_CLUMPING_R_SQUARED = 0.2;
    private static final double DEFAULT_MAF_THRESHOLD = 0;
    private static final List<Double> DEFAULT_P_VALUE_THRESHOLDS = new ArrayList<>(Collections.singletonList(1e-5));
    private static final int DEFAULT_NUMBER_OF_FOLDS = 10;
    private static final String[] HSA_DEFAULT_SEQUENCES =
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                    "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

    private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();//Might be changed
    private static final Logger LOGGER = Logger.getLogger(PGSBasedMixupMapperOptions.class);

    private final boolean debugMode;
    private final boolean calculateNewGenomeWideAssociationsEnabled;
    private final boolean writeNewGenomeWideAssociationsEnabled;
    private final double minorAlleleFrequencyThreshold;
    private final double rSquared;
    private final int folds;
    private final String forceSeqName;
    private final String gwasSummaryStatisticsPhenotypeCouplingFile;
    private final String genotypeToPhenotypeSampleCouplingFile;
    private final List<Integer> windowSize;
    private final List<Double> pValueThresholds;
    private final List<String> genomicRangesToExclude;
    private final RandomAccessGenotypeDataReaderFormats inputGenotypeType;
    private final String[] inputGenotypePath;
    private final String[] sequences;
    private final Path gwasSummaryStatisticsPath;
    private final File outputBasePath;
    private final File inputPhenotypePath;
    private final File logFile;
    private final File debugFolder;
    private final CommandLine commandLine;



    static {

        OPTIONS = new Options();

        OptionBuilder.withArgName("string");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The PGS calculation mode");
        OptionBuilder.withLongOpt("PGSmode");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("pgs"));

        OptionBuilder.withArgName("basePath");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("The input genotype file(s)");
        OptionBuilder.withLongOpt("input");
        OPTIONS.addOption(OptionBuilder.create("i"));

        OptionBuilder.withArgName("type");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The input genotype data type. " +
                "If not defined will attempt to automatically select the first matching dataset on the specified path\n"
                + "* PED_MAP - plink PED MAP files.\n"
                + "* PLINK_BED - plink BED BIM FAM files.\n"
                + "* VCF - bgziped vcf with tabix index file\n"
                + "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
                + "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
                + "* GEN - Oxford .gen & .sample\n"
                + "* BGEN - Oxford .bgen & optionally .sample\n"
                + "* TRITYPER - TriTyper format folder");
        OptionBuilder.withLongOpt("inputType");
        OPTIONS.addOption(OptionBuilder.create("I"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The output path");
        OptionBuilder.withLongOpt("output");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("o"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("The input phenotype path");
        OptionBuilder.withLongOpt("inputPhenotype");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("p"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("GWAS summary statistics path");
        OptionBuilder.withLongOpt("gwasStatsPath");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create('s'));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("CSV file representing the coupling between " +
                "GWAS summary statistics files and phenotypes.");
        OptionBuilder.withLongOpt("gwasStatsToPhenCoupling");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("wpc"));

        OptionBuilder.withArgName("path");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("CSV file representing the current coupling between genotype and phenotype samples.");
        OptionBuilder.withLongOpt("genToPhenSampleCoupling");
        OptionBuilder.isRequired();
        OPTIONS.addOption(OptionBuilder.create("gpc"));

        OptionBuilder.withArgName("string");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Shapeit2 does not output the sequence name in the first column of " +
                "the haplotype file. Use this option to force the chromosome for all variants. " +
                "This option is only valid in combination with --inputType SHAPEIT2");
        OptionBuilder.withLongOpt("forceChr");
        OPTIONS.addOption(OptionBuilder.create("f"));

        OptionBuilder.withArgName("int");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Maximum number of calculation threads");
        OptionBuilder.withLongOpt("threads");
        OPTIONS.addOption(OptionBuilder.create("t"));

        OptionBuilder.withArgName("String");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Exclude genomic range(s) from the risk score calculation. " +
                "Range needs to be specified as: \"6:101-110 6:250000-350000. Warning: " +
                "Chr name must be specified as expected in the genotype dataset.");
        OptionBuilder.withLongOpt("excludeRange");
        OPTIONS.addOption(OptionBuilder.create("er"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Activate debug mode. " +
                "This will result in a more verbose log file and will save many intermediate results to files. " +
                "Not recommended for large analysis.");
        OptionBuilder.withLongOpt("debug");
        OPTIONS.addOption(OptionBuilder.create("d"));

        OptionBuilder.withArgName("double");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Minimum MAF");
        OptionBuilder.withLongOpt("maf");
        OPTIONS.addOption(OptionBuilder.create("maf"));

        OptionBuilder.withArgName("double");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("P-value thresholds for genetic risk score inclusion, " +
                "colon separated should be ordered from most stringent to least stringent.");
        OptionBuilder.withLongOpt("pValueThresholds");
        OPTIONS.addOption(OptionBuilder.create("pv"));

        OptionBuilder.withArgName("integer[:integer]");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Window sizes for clumping, if given two window-sizes (colon separated), " +
                "a two step window approach is used.");
        OptionBuilder.withLongOpt("windowSizes");
        OPTIONS.addOption(OptionBuilder.create("bp"));

        OptionBuilder.withArgName("double");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("R2 for clumping.");
        OptionBuilder.withLongOpt("rSquared");
        OPTIONS.addOption(OptionBuilder.create("r2"));

        OptionBuilder.withArgName("strings");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Chromosome identifiers to expand genotype input with");
        OptionBuilder.withLongOpt("chromosomes");
        OPTIONS.addOption(OptionBuilder.create('c'));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Calculates genome-wide associations for the phenotypes in the " +
                "CSV coupling these phenotypes to filenames for use in calculating polygenic scores. " +
                "Files will be overwritten!");
        OptionBuilder.withLongOpt("calculateNewAssociations");
        OPTIONS.addOption(OptionBuilder.create("n"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Write summary statistics of newly calculated genome wide associations. " +
                "Option is ignored without -n / --calculateNewAssociations");
        OptionBuilder.withLongOpt("writeSummaryStatistics");
        OPTIONS.addOption(OptionBuilder.create("w"));

        OptionBuilder.withArgName("int");
        OptionBuilder.withDescription("K folds to iterate over and split samples.");
        OptionBuilder.withLongOpt("folds");
        OPTIONS.addOption(OptionBuilder.create('K'));
    }

    PGSBasedMixupMapperOptions(String... args) throws ParseException {

        // Parse the raw command line input
        final CommandLineParser parser = new PosixParser();
        commandLine = parser.parse(OPTIONS, args, false);
        setNumberOfThreadsToUse(commandLine);

        outputBasePath = getOutputBasePath(commandLine);
        logFile = new File(outputBasePath + ".log");
        debugMode = commandLine.hasOption('d');
        debugFolder = new File(outputBasePath + "_debugFiles");
        gwasSummaryStatisticsPath = Paths.get(commandLine.getOptionValue('s'));
        gwasSummaryStatisticsPhenotypeCouplingFile = commandLine.getOptionValue("wpc");
        genotypeToPhenotypeSampleCouplingFile = commandLine.getOptionValue("gpc");

        // Get the path and type for the genotype data
        inputGenotypePath = getInputGenotypePath(commandLine);
        inputGenotypeType = getInputGenotypeType(commandLine);

        // Get the force seq name
        forceSeqName = getForceSeqName(commandLine);

        // Is the custom GWA enabled?
        calculateNewGenomeWideAssociationsEnabled = commandLine.hasOption("n");
        writeNewGenomeWideAssociationsEnabled = commandLine.hasOption("w");

        // Parse the command line option representing sequence identifiers to expand the genotype input with.
        sequences = parseSequences(commandLine);

        // Get the R squared from the command line arguments.
        rSquared = getRSquared(commandLine);

        // Get the window sizes from the command line arguments.
        windowSize = parseWindowSizes(commandLine);

        // Get the P-value thresholds from the command line arguments.
        pValueThresholds = parsePValueThresholds(commandLine);

        // Get the input phenotype path.
        inputPhenotypePath = parseInputPhenotypePath(commandLine);

        // Parse the genomic ranges to exclude in PGS calculations
        genomicRangesToExclude = parseRangesToExclude(commandLine);

        // Parse the minor allele frequency threshold to use for the input genotype data
        minorAlleleFrequencyThreshold = parseMinorAlleleFrequencyThreshold(commandLine);

        folds = parseKFolds(commandLine);
    }

    private int parseKFolds(CommandLine commandLine) throws ParseException {
        int folds = DEFAULT_NUMBER_OF_FOLDS;

        if (commandLine.hasOption("K")) {
            try {
                folds = Integer.parseInt(commandLine.getOptionValue("K"));
            } catch (NumberFormatException e) {
                throw new ParseException(String.format(
                        "Error parsing --folds / -K: \"%s\" is not an integer", commandLine.getOptionValue("K")));
            }
        }
        return folds;
    }

    /**
     * Parses the minor allele frequency thresholds from the provided command line arguments.
     *
     * @param commandLine the command line that could contain the minor allele frequency threshold in option "maf".
     * @return a double with the MAF threshold provided or its the default value.
     * @throws ParseException the MAF value was not able to be parsed to a double.
     */
    private double parseMinorAlleleFrequencyThreshold(CommandLine commandLine) throws ParseException {
        // Initialize minor allele frequency to return
        double mafThreshold = DEFAULT_MAF_THRESHOLD;

        if (commandLine.hasOption("maf")) {
            try {
                mafThreshold =  Double.parseDouble(commandLine.getOptionValue("maf"));
            } catch (NumberFormatException e) {
                throw new ParseException(String.format(
                        "Error parsing --maf: \"%s\" is not a double", commandLine.getOptionValue("maf")));
            }
        }
        return mafThreshold;
    }

    /**
     * Gets the ranges to exclude from the command line.
     *
     * @param commandLine The commandline that could contain the ranges to exclude in option "excludeRange".
     * @return A list of ranges to exclude. Could be empty if the option was not present.
     */
    private List<String> parseRangesToExclude(CommandLine commandLine) {
        List<String> ranges = new ArrayList<>();

        if (commandLine.hasOption("er")) {
            // initialise the member variable
            Collections.addAll(ranges, commandLine.getOptionValues("excludeRange"));
        }
        return ranges;
    }

    /**
     * Gets the P-values thresholds from the provided command line arguments.
     *
     * @param commandLine the command line that could contain the P-value thresholds in option 'pv'.
     * @return a list with every P-value threshold provided or the default value of 1e-5.
     * @throws ParseException if any of the P-values were not able to be parsed to a double.
     */
    private List<Double> parsePValueThresholds(CommandLine commandLine) throws ParseException {

        // Initialize an array list with
        List<Double> pValueThresholds = new ArrayList<>();

        if (!commandLine.hasOption("pv")) {
            return DEFAULT_P_VALUE_THRESHOLDS;
        }

        // Loop through the P-value thresholds
        for (String pValue : commandLine.getOptionValue("pv").split(":")) {
            // Parse the P-value to a double
            try {
                double e = Double.parseDouble(pValue);
                pValueThresholds.add(e);
            // If the supposed P-value cannot be parsed to a double, a NumberFormatException is thrown,
            // catch this and throw a parse exception.
            } catch (NumberFormatException e) {
                throw new ParseException(String.format(
                        "Error parsing -pv / --pValueThresholds: \"%s\" could not be parsed to a double",
                        pValue));
            }
        }
        return pValueThresholds;
    }

    /**
     * Gets the window sizes from the provided command line arguments.
     *
     * @param commandLine the command line that could contain the window sizes in option 'bp'.
     * @return a list with every window size provided or the default value of 50000.
     * @throws ParseException if any of the window size values were not able to be parsed to an integer.
     */
    private List<Integer> parseWindowSizes(CommandLine commandLine) throws ParseException {

        // Get the window size parameter or the default value of 500.000
        String windowSizeArgument = commandLine.getOptionValue("bp", String.valueOf(500000));
        // Split the given option
        String[] split = windowSizeArgument.split(":");
        // Initialize the window sizes list
        List<Integer> windowSizes = new ArrayList<>();

        // Loop through the window sizes
        for (String splitWindowSize : split) {
            // Try to parse the window size to a integer
            try {
                windowSizes.add(Integer.parseInt(splitWindowSize));
            // If the supposed window size cannot be parsed to an integer, a NumberFormatException is thrown,
            // catch this and throw a parse exception.
            } catch (NumberFormatException e) {
                throw new ParseException(String.format("Error parsing -bp / --windowSizes: \"%s\" is not a valid integer",
                        splitWindowSize));
            }
        }
        return windowSizes;
    }

    /**
     * Gets the R-squared threshold for clumping from the provided command line arguments.
     *
     * @param commandLine the command line that could contain the R-squared option 'r2'.
     * @return the provided R-squared value or the default.
     * @throws ParseException If the value was not able to be parsed to a double.
     */
    private double getRSquared(CommandLine commandLine) throws ParseException {
    	if (commandLine.hasOption("r2")) {
			try {
			    // Try to parse the r2
				return Double.parseDouble(commandLine.getOptionValue("r2"));
            // Throw a parse exception if the r2 cannot be parsed to a double
			} catch (NumberFormatException e) {
				throw new ParseException(String.format("Error parsing -r2 / --rSquared \"%s\" is not a valid double",
						commandLine.getOptionValue("r2")));
			}
		}
    	return DEFAULT_CLUMPING_R_SQUARED;
	}

    /**
     * Method parsing a the chromosomes option from the given options.
     *
     * @param commandLine the command line that could contain the chromosomes option.
     * @return the provided chromosomes or the default human sequences.
     */
	private String[] parseSequences(CommandLine commandLine) {
        // If the command line has a chromosomes option, return this.
        if (commandLine.hasOption("chromosomes")) {
            return commandLine.getOptionValues("chromosomes");
        }
        // If the command line does not have a chromosomes option return the default sequences
        return HSA_DEFAULT_SEQUENCES;
    }

    /**
     * Method that gets the input phenotype filepath from the command line.
     *
     * @param commandLine A command line that contains the input phenotype filepath in option "inputPhenotype"
     * @return the provided input phenotype filepath
     * @throws ParseException if either the option "inputPhenotype" was not present in the command line,
     * or if the provided path does not point to an existing file.
     */
    private File parseInputPhenotypePath(CommandLine commandLine) throws ParseException {
        if (!commandLine.hasOption("inputPhenotype")) {
            throw new ParseException("--inputPhenotype not specified");
        } else {
            File inputPhenotypePath = new File(commandLine.getOptionValue("inputPhenotype"));
            if (inputPhenotypePath.isFile()) {
                return inputPhenotypePath;
            }
            throw new ParseException(String.format("Input phenotype file \"%s\"does not exist",
                    commandLine.getOptionValue("inputPhenotype")));
        }
    }

    /**
     * Method that gets the output base path from the command line.
     *
     * @param commandLine The command line that contains the output base path in '-o'
     * @return the output base path as a file.
     * @throws ParseException if either the option '-o' was not present in the command line,
     * or if the provided path does not point to an existing directory.
     */
    private File getOutputBasePath(CommandLine commandLine) throws ParseException {
        if (!commandLine.hasOption('o')) {
            throw new ParseException("-o / --output not specified");
        }

        File outputBasePath = new File(commandLine.getOptionValue('o'));
        if (outputBasePath.isDirectory()) {
            throw new ParseException(String.format("Specified output path '%s' is a directory. " +
                            "Please include a prefix for the output files.",
                    outputBasePath.toString()));
        }
        return outputBasePath;
    }

    private String getForceSeqName(CommandLine commandLine) throws ParseException {
        boolean isForceSeqPresent = commandLine.hasOption('f');
        if (
                isForceSeqPresent && this.getInputGenotypeType() != null &&
                        this.getInputGenotypeType() != RandomAccessGenotypeDataReaderFormats.SHAPEIT2 &&
                        this.getInputGenotypeType() != RandomAccessGenotypeDataReaderFormats.GEN
        ) {
            throw new ParseException("Error cannot force sequence name of: " + inputGenotypeType.getName());
        }
        return isForceSeqPresent ? commandLine.getOptionValue('f') : null;
    }

    private void setNumberOfThreadsToUse(CommandLine commandLine) throws ParseException {
        if (commandLine.hasOption('t')) {
            try {
                numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
                System.setProperty("Djava.util.concurrent.ForkJoinPool.common.parallelism", commandLine.getOptionValue('t'));
                ConcurrencyUtils.setNumberOfThreads(numberOfThreadsToUse);
            } catch (NumberFormatException e) {
                throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
            }
        }
    }

    private String[] getInputGenotypePath(CommandLine commandLine) throws ParseException {
        if (!commandLine.hasOption('i')) {
            throw new ParseException("--input not specified");
        }
        return commandLine.getOptionValues('i');
    }

    private RandomAccessGenotypeDataReaderFormats getInputGenotypeType(CommandLine commandLine) throws ParseException {
        RandomAccessGenotypeDataReaderFormats inputGenotypeType;

        try {
            if (commandLine.hasOption('I')) {
                inputGenotypeType = RandomAccessGenotypeDataReaderFormats.valueOfSmart(commandLine.getOptionValue('I').toUpperCase());
            } else {
                if (inputGenotypePath[0].endsWith(".vcf")) {
                    throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
                }
                try {
                    inputGenotypeType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(inputGenotypePath);
                } catch (GenotypeDataException e) {
                    throw new ParseException("Unable to determine reference type based on specified path. Please specify --inputType");
                }
            }

        } catch (IllegalArgumentException e) {
            throw new ParseException("Error parsing --inputType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
        }

        return inputGenotypeType;
    }

    public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }

    public void printOptions() {

        LOGGER.info("Supplied options:");

        LOGGER.info(" * Ouput path: " + outputBasePath.getAbsolutePath());

        LOGGER.info(" * Debug mode: " + (debugMode ? "on (this will result in many intermediate output files)" : "off"));

    }

    /**
     * Method that returns the input genotype paths expanded for every sequence
     * if one of these paths contain a '#'. The input genotype paths are copied
     * for every sequence with which the '#' is subsequently replaced by the
     * sequence name.
     *
     * @return a map with sequence names as the keys and a corresponding array
     * of input paths as values. the sequence name is null if a '#' is not present
     * and force seq name was not set.
     */
    public Map<String, String[]> getInputGenotypePaths() {

        // Initialize the output 2d array
        Map<String, String[]> expandedInputGenotypePaths = new LinkedHashMap<>();

        // Check if one of the input paths for the genotype data contains a '#',
        // If this is the case expand the input genotype paths.
        if (Arrays.stream(inputGenotypePath).anyMatch(path -> path.contains("#"))) {

            // Loop through the sequences and the genotype paths, replacing every
            // '#' with the sequence name.
            for (String sequence : sequences) {
                String[] strings = new String[inputGenotypePath.length];
                for (int i = 0; i < inputGenotypePath.length; i++) {
                    String filePath = inputGenotypePath[i];
                    strings[i] = filePath.replace("#", sequence);
                }
                expandedInputGenotypePaths.put(sequence, strings);
            }
        } else {
            // If the paths do not contain a '#', return the map with size 1.
            expandedInputGenotypePaths.put(forceSeqName, inputGenotypePath);
        }

        return expandedInputGenotypePaths;
    }

    public RandomAccessGenotypeDataReaderFormats getInputGenotypeType() {
        return inputGenotypeType;
    }

    public File getOutputBasePath() {
        return outputBasePath;
    }

    public File getLogFile() {
        return logFile;
    }

    public File getInputPhenotypePath() {
        return inputPhenotypePath;
    }

    public File getDebugFolder() {
        return debugFolder;
    }

    public double getMinorAlleleFrequencyThreshold() {
        return minorAlleleFrequencyThreshold;
    }

    public String getForceSeqName() {
        return forceSeqName;
    }

    public Path getGwasSummaryStatisticsPath() {
        return gwasSummaryStatisticsPath;
    }

    public String getGwasSummaryStatisticsPhenotypeCouplingFile() {
        return gwasSummaryStatisticsPhenotypeCouplingFile;
    }

    public String getGenotypeToPhenotypeSampleCouplingFile() {
        return genotypeToPhenotypeSampleCouplingFile;
    }

    public double getrSquared() {
        return rSquared;
    }

    public List<Integer> getWindowSize() {
        return windowSize;
    }

    public List<Double> getpValueThresholds() {
        return pValueThresholds;
    }

    public String[] getGenomicRangesToExclude() {
        return genomicRangesToExclude.stream().toArray(String[] ::new);
    }

    public boolean isDebugMode() {
        return debugMode;
    }

    public String[] getInputGenotypePath() {
        return inputGenotypePath;
    }

    public boolean isCalculateNewGenomeWideAssociationsEnabled() {
        return calculateNewGenomeWideAssociationsEnabled;
    }

    public boolean shouldForceSeqName() {
        return forceSeqName != null;
    }

    public boolean isWriteNewGenomeWideAssociationsEnabled() {
        return writeNewGenomeWideAssociationsEnabled;
    }

    protected static Options getOptions() {
        return OPTIONS;
    }

    public CommandLine getCommandLine() {
        return commandLine;
    }

    public int getFolds() {
        return folds;
    }
}
