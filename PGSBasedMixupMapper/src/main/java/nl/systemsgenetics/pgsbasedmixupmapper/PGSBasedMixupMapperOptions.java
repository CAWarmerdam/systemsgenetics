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
import java.util.*;

/**
 * @author patri
 */
public class PGSBasedMixupMapperOptions {

    private static final double DEFAULT_CLUMPING_R_SQUARED = 0.2;
    private static final String[] HSA_DEFAULT_SEQUENCES =
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                    "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

    private static final Options OPTIONS;
    private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();//Might be changed
    private static final Logger LOGGER = Logger.getLogger(PGSBasedMixupMapperOptions.class);

    private final String[] inputGenotypePath;
    private final RandomAccessGenotypeDataReaderFormats inputGenotypeType;
    private final File outputBasePath;
    private final String gwasSummaryStatisticsPath;
    private final boolean debugMode;
    private final boolean calculateNewGenomeWideAssociationsEnabled;
    private final boolean writeNewGenomeWideAssociationsEnabled;
    private final double mafFilter;
    private final double rSquared;
    private final String forceSeqName;
    private final String gwasSummaryStatisticsPhenotypeCouplingFile;
    private final List<Integer> windowSize;
    private final List<Double> pValueThresholds;
    private final String[] genomicRangesToExclude;

    private final String[] sequences;

    private final File inputPhenotypePath;
    private final File logFile;
    private final File debugFolder;
    private final String genotypeToPhenotypeSampleCouplingFile;


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
        OptionBuilder.withDescription("The input genotype data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
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
        OptionBuilder.withDescription("CSV file representing the coupling between GWAS summary statistics files and phenotypes.");
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
                "Range needs to be specified as: \"6:101-110;6:250000-350000. Warning: " +
                "Chr name must be specified as expected in the genotype dataset.");
        OptionBuilder.withLongOpt("excludeRange");
        OPTIONS.addOption(OptionBuilder.create("er"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Activate debug mode. This will result in a more verbose log file and will save many intermediate results to files. Not recommended for large analysis.");
        OptionBuilder.withLongOpt("debug");
        OPTIONS.addOption(OptionBuilder.create("d"));

        OptionBuilder.withArgName("double");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Minimum MAF");
        OptionBuilder.withLongOpt("maf");
        OPTIONS.addOption(OptionBuilder.create("maf"));

        OptionBuilder.withArgName("boolean");
        OptionBuilder.withDescription("Quantile normalize the permutations gene p-values before pathway enrichments to fit the real gene p-value distribution");
        OptionBuilder.withLongOpt("qnorm");
        OPTIONS.addOption(OptionBuilder.create("qn"));

        OptionBuilder.withArgName("double");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("P-value thresholds for genetic risk score inclusion, " +
                "colon separated should be ordered from most stringent to least stringent.");
        OptionBuilder.withLongOpt("pValueThresholds");
        OPTIONS.addOption(OptionBuilder.create("pv"));

        OptionBuilder.withArgName("integer[:integer]");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Window size for clumping, if given two window-sizes (colon separated), " +
                "a two step window approach is used.");
        OptionBuilder.withLongOpt("windowSize");
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
    }

    PGSBasedMixupMapperOptions(String... args) throws ParseException {

        // Parse the raw command line input
        final CommandLineParser parser = new PosixParser();
        CommandLine commandLine = parser.parse(OPTIONS, args, false);
        setNumberOfThreadsToUse(commandLine);

        outputBasePath = getOutputBasePath(commandLine);
        logFile = new File(outputBasePath + ".log");
        debugMode = commandLine.hasOption('d');
        debugFolder = new File(outputBasePath + "_debugFiles");
        gwasSummaryStatisticsPath = commandLine.getOptionValue('s');
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
        sequences = getSequences(commandLine);

        // Get the R squared from the command line arguments.
        rSquared = getRSquared(commandLine);

        // Get the window sizes from the command line arguments
        windowSize = getWindowSizes(commandLine);

        try {
            pValueThresholds = new ArrayList<>();
            for (String pvalue : commandLine.getOptionValue("pv", String.valueOf(1e-5)).split(":")) {
                double e = Double.parseDouble(pvalue);
                pValueThresholds.add(e);
            }
        } catch (NumberFormatException e) {
            throw new ParseException(String.format("Error parsing --pv \"%s\" could not be parsed to valid doubles",
                    commandLine.getOptionValue("pv")));
        }

        if (commandLine.hasOption("excludeRange") || commandLine.hasOption("er")) {
            // initialise the member variable
            genomicRangesToExclude = commandLine.getOptionValue("excludeRange").split(";");
        } else {
            genomicRangesToExclude = new String[]{};
        }

        this.inputPhenotypePath = getInputPhenotypePath(commandLine);

        if (commandLine.hasOption("maf")) {
            try {
                mafFilter = Double.parseDouble(commandLine.getOptionValue("maf"));
            } catch (NumberFormatException e) {
                throw new ParseException("Error parsing --maf \"" + commandLine.getOptionValue("maf") + "\" is not an double");
            }
        } else {
            mafFilter = 0;
        }
    }

    /**
     * Gets the window sizes from the provided command line arguments.
     *
     * @param commandLine the command line that could contain the window sizes in option 'bp'.
     * @return a list with every window size provided or the default value of 50000.
     * @throws ParseException If any of the window size values were not able to be parsed to an integer.
     */
    private List<Integer> getWindowSizes(CommandLine commandLine) throws ParseException {
        String windowSizeArgument = commandLine.getOptionValue("bp", String.valueOf(500000));
        String[] split = windowSizeArgument.split(":");
        List<Integer> windowSizes = new ArrayList<>();
        for (String splitWindowSize :
                split) {
            try {
                windowSizes.add(Integer.parseInt(splitWindowSize));
            } catch (NumberFormatException e) {
                throw new ParseException(String.format("Error parsing --bp \"%s\" is not a valid integer",
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
				return Double.parseDouble(commandLine.getOptionValue("r2"));
			} catch (NumberFormatException e) {
				throw new ParseException(String.format("Error parsing --r2 \"%s\" is not a valid double",
						commandLine.getOptionValue("r2")));
			}
		}
    	return DEFAULT_CLUMPING_R_SQUARED;
	}

	private String[] getSequences(CommandLine commandLine) {
        if (commandLine.hasOption("chromosomes")) {
            return commandLine.getOptionValues("chromosomes");
        }
        return HSA_DEFAULT_SEQUENCES;
    }

    private File getInputPhenotypePath(CommandLine commandLine) throws ParseException {
        if (!commandLine.hasOption("inputPhenotype")) {
            throw new ParseException("--inputPhenotype not specified");
        } else {
            File inputPhenotypePath = new File(commandLine.getOptionValue("inputPhenotype"));
            if (inputPhenotypePath.exists()) {
                return inputPhenotypePath;
            }
            throw new ParseException(String.format("Input phenotype file \"%s\"does not exist",
                    commandLine.getOptionValue("inputPhenotype")));
        }
    }

    private File getOutputBasePath(CommandLine commandLine) throws ParseException {
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
                isForceSeqPresent &&
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

    public double getMafFilter() {
        return mafFilter;
    }

    public String getForceSeqName() {
        return forceSeqName;
    }

    public String getGwasSummaryStatisticsPath() {
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
        return genomicRangesToExclude;
    }

    public boolean isDebugMode() {
        return debugMode;
    }

    public String[] getInputGenotypePath() {
        return inputGenotypePath;
    }

    public Map<String, String[]> getInputGenotypePaths() {

        // Initialize the output 2d array
        Map<String, String[]> expandedInputGenotypePaths = new LinkedHashMap<>();

        if (Arrays.stream(inputGenotypePath).anyMatch(path -> path.contains("#"))) {
            for (String sequence : sequences) {
                String[] strings = new String[inputGenotypePath.length];
                for (int i = 0; i < inputGenotypePath.length; i++) {
                    String filePath = inputGenotypePath[i];
                    strings[i] = filePath.replace("#", sequence);
                }
                expandedInputGenotypePaths.put(sequence, strings);
            }
        } else {
            expandedInputGenotypePaths.put(forceSeqName, inputGenotypePath);
        }

        return expandedInputGenotypePaths;
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
}
