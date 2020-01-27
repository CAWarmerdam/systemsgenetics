/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.pgsbasedmixupmapper;

import edu.emory.mathcs.utils.ConcurrencyUtils;
import nl.systemsgenetics.polygenicscorecalculator.PolyGenicScoreCalculatorMode;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author patri
 */
public class PGSBasedMixupMapperOptions {

	private static final Options OPTIONS;
	private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();//Might be changed
	private static final Logger LOGGER = Logger.getLogger(PGSBasedMixupMapperOptions.class);

	private final PolyGenicScoreCalculatorMode pgsMode;

	private final String[] inputGenotypePath;
	private final RandomAccessGenotypeDataReaderFormats inputGenotypeType;
	private final File outputBasePath;
	private final String inputPhenotypePath;
	private final String gwasSummaryStatisticsPath;
	private final File logFile;
	private final boolean debugMode;
	private final File debugFolder;
	private final double mafFilter;
	private final String forceSeqName;
	private final String gwasSummaryStatisticsPhenotypeCouplingFile;
	private final double rSquared;
	private final int windowSize;
	private final List<Double> pValueThresholds;
	private final String[] genomicRangesToExclude;
	private String genotypeToPhenotypeSampleCouplingFile;

	public boolean isDebugMode() {
		return debugMode;
	}

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
		OptionBuilder.withDescription("The input genotype file");
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

		OptionBuilder.withArgName("integer");
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
	}

	public PGSBasedMixupMapperOptions(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		CommandLine commandLine = parser.parse(OPTIONS, args, false);
		setNumberOfThreadsToUse(commandLine);

		pgsMode = getPolyGenicScoreCalculatorMode(commandLine);
		outputBasePath = getOutputBasePath(commandLine);
		logFile = new File(outputBasePath + ".log");
		debugMode = commandLine.hasOption('d');
		gwasSummaryStatisticsPath = commandLine.getOptionValue('s');
		gwasSummaryStatisticsPhenotypeCouplingFile = commandLine.getOptionValue("wpc");
		genotypeToPhenotypeSampleCouplingFile = commandLine.getOptionValue("gpc");
		debugFolder = new File(outputBasePath + "_debugFiles");

		inputGenotypePath = getInputGenotypePath(commandLine);
		inputGenotypeType = getInputGenotypeType(commandLine);

		forceSeqName = getForceSeqName(commandLine);

		try {
			rSquared = Double.parseDouble(commandLine.getOptionValue("r2", String.valueOf(0.2)));
		} catch (NumberFormatException e) {
			throw new ParseException(String.format("Error parsing --r2 \"%s\" is not a valid double",
					commandLine.getOptionValue("r2")));
		}

		try {
			windowSize = Integer.parseInt(commandLine.getOptionValue("bp", String.valueOf(500000)));
		} catch (NumberFormatException e) {
			throw new ParseException(String.format("Error parsing --bp \"%s\" is not a valid integer",
					commandLine.getOptionValue("bp")));
		}

		try {
			pValueThresholds = new ArrayList<>();
			for (String pvalue : commandLine.getOptionValue("pv", String.valueOf(1e-5)).split(":")) {
				double e = -Math.log10(Double.parseDouble(pvalue));
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

		if (!commandLine.hasOption("inputPhenotype")) {
			throw new ParseException("--inputPhenotype not specified");
		} else {
			inputPhenotypePath = commandLine.getOptionValue("inputPhenotype");
		}

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

	private PolyGenicScoreCalculatorMode getPolyGenicScoreCalculatorMode(CommandLine commandLine) throws ParseException {
		try {
			return PolyGenicScoreCalculatorMode.valueOf(commandLine.getOptionValue("pgs").toUpperCase());
		} catch (IllegalArgumentException e) {
			throw new ParseException(String.format("Error parsing --PGSmode \"%s\" is not a valid mode",
					commandLine.getOptionValue("pgs")));
		}
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

	public String[] getInputGenotypePath() {
		return inputGenotypePath;
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

	public String getInputPhenotypePath() {
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

	public PolyGenicScoreCalculatorMode getPgsMode() {
		return pgsMode;
	}

	public double getrSquared() {
		return rSquared;
	}

	public int getWindowSize() {
		return windowSize;
	}

	public List<Double> getpValueThresholds() {
		return pValueThresholds;
	}

	public String[] getGenomicRangesToExclude() {
		return genomicRangesToExclude;
	}

}
