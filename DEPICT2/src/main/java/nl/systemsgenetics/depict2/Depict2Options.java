/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

/**
 *
 * @author patri
 */
public class Depict2Options {

	private static final Options OPTIONS;
	private static int numberOfThreadsToUse = Runtime.getRuntime().availableProcessors();

	private final Depict2Mode mode;

	private final String[] genotypeBasePath;
	private final RandomAccessGenotypeDataReaderFormats genotypeType;
	private final File outputFile;
	private final File geneInfoFile;
	private final String gwasZscoreMatrixPath;
	private final int numberOfPermutations;
	private final int windowExtend;
	private final double maxRBetweenVariants;

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("On of the following modes:\n"
				+ "* RUN - Run the DEPICT2 prioritization.\n"
				+ "* CONVERT_TXT - Convert a txt z-score matrix to binary. Use --gwas and --output\n"
				+ "* CONVERT_EQTL - Convert binary matrix with eQTL z-scores from our pipeline. Use --gwas and --output");
		OptionBuilder.withLongOpt("mode");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("m"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("GWAS Z-sccore binary matrix. Rows variants, Cols phenotypes. Without .dat");
		OptionBuilder.withLongOpt("gwas");
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The reference genotypes");
		OptionBuilder.withLongOpt("referenceGenotypes");
		OPTIONS.addOption(OptionBuilder.create("r"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The reference genotype data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("referenceGenotypeFormat");
		OPTIONS.addOption(OptionBuilder.create("R"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The output path");
		OptionBuilder.withLongOpt("output");
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Maximum number of calculation threads");
		OptionBuilder.withLongOpt("threads");
		OPTIONS.addOption(OptionBuilder.create("t"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Number of permutations");
		OptionBuilder.withLongOpt("permutations");
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Number of bases to add left and right of gene window");
		OptionBuilder.withLongOpt("window");
		OPTIONS.addOption(OptionBuilder.create("w"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Max correlation between variants to use (recommend = 0.95)");
		OptionBuilder.withLongOpt("variantCorrelation");
		OPTIONS.addOption(OptionBuilder.create("v"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("File with gene info. col1: geneName (ensg) col2: chr col3: startPos col4: stopPos");
		OptionBuilder.withLongOpt("genes");
		OPTIONS.addOption(OptionBuilder.create("ge"));

	}

	public Depict2Options(String... args) throws ParseException {

		final CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);

		if (commandLine.hasOption('t')) {
			try {
				numberOfThreadsToUse = Integer.parseInt(commandLine.getOptionValue('t'));
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --threads \"" + commandLine.getOptionValue('t') + "\" is not an int");
			}
		}

		outputFile = new File(commandLine.getOptionValue('o'));

		try {
			mode = Depict2Mode.valueOf(commandLine.getOptionValue("m"));
		} catch (IllegalArgumentException e) {
			throw new ParseException("Error parsing --mode \"" + commandLine.getOptionValue("m") + "\" is not a valid mode");
		}

		if (mode == Depict2Mode.CONVERT_TXT || mode == Depict2Mode.RUN || mode == Depict2Mode.CONVERT_EQTL) {

			if (!commandLine.hasOption("g")) {
				throw new ParseException("Please provide --gwas for mode: " + mode.name());
			}

			gwasZscoreMatrixPath = commandLine.getOptionValue('g');
		} else {
			gwasZscoreMatrixPath = null;
		}

		if (mode == Depict2Mode.RUN) {

			if (!commandLine.hasOption("p")) {
				throw new ParseException("--permutations not specified");
			} else {

				try {
					numberOfPermutations = Integer.parseInt(commandLine.getOptionValue('p'));
				} catch (NumberFormatException e) {
					throw new ParseException("Error parsing --permutations \"" + commandLine.getOptionValue('p') + "\" is not an int");
				}
			}

			if (!commandLine.hasOption("w")) {
				throw new ParseException("--window not specified");
			} else {
				try {
					windowExtend = Integer.parseInt(commandLine.getOptionValue('w'));
				} catch (NumberFormatException e) {
					throw new ParseException("Error parsing --window \"" + commandLine.getOptionValue('w') + "\" is not an int");
				}
			}

			if (!commandLine.hasOption("v")) {
				throw new ParseException("--variantCorrelation not specified");
			} else {
				try {
					maxRBetweenVariants = Double.parseDouble(commandLine.getOptionValue('v'));
				} catch (NumberFormatException e) {
					throw new ParseException("Error parsing --variantCorrelation \"" + commandLine.getOptionValue('v') + "\" is not an double");
				}
			}

			if (!commandLine.hasOption("ge")) {
				throw new ParseException("--genes not specified");
			} else {
				geneInfoFile = new File(commandLine.getOptionValue("ge"));
			}

			if (!commandLine.hasOption('r')) {

				throw new ParseException("--referenceGenotypes not specified");

			} else {

				genotypeBasePath = commandLine.getOptionValues('r');

				try {
					if (commandLine.hasOption('R')) {
						genotypeType = RandomAccessGenotypeDataReaderFormats.valueOfSmart(commandLine.getOptionValue('R').toUpperCase());
					} else {
						if (genotypeBasePath[0].endsWith(".vcf")) {
							throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
						}
						try {
							genotypeType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(genotypeBasePath);
						} catch (GenotypeDataException e) {
							throw new ParseException("Unable to determine reference type based on specified path. Please specify --refType");
						}
					}

				} catch (IllegalArgumentException e) {
					throw new ParseException("Error parsing --refType \"" + commandLine.getOptionValue('R') + "\" is not a valid reference data format");
				}

			}
		} else {
			genotypeBasePath = null;
			genotypeType = null;
			geneInfoFile = null;
			maxRBetweenVariants = 0d;
			numberOfPermutations = 0;
			windowExtend = 0;
		}

	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
	}

	public void printOptions() {

		System.out.println(" - Mode: " + mode.name());
		System.out.println(" - Ouput path: " + outputFile.getPath());

		switch (mode) {
			case CONVERT_EQTL:
				System.out.println(" - eQTL Z-score matrix: " + gwasZscoreMatrixPath);
				break;
			case CONVERT_TXT:
				System.out.println(" - Gwas Z-score matrix: " + gwasZscoreMatrixPath);
				break;
			case RUN:
				System.out.println(" - Gwas Z-score matrix: " + gwasZscoreMatrixPath);

				if (genotypeBasePath != null) {
					StringBuilder genotypeBasePaths = new StringBuilder();
					for (String path : genotypeBasePath) {
						genotypeBasePaths.append(path);
						genotypeBasePaths.append(' ');
					}
					System.out.println(" - Reference genotype data: " + genotypeBasePaths);
					System.out.println(" - Reference genotype data type: " + genotypeType.getName());
				}
				System.out.println(" - Gene window extend in bases: " + windowExtend);
				System.out.println(" - Number of permutations: " + numberOfPermutations);
				System.out.println(" - Max correlation between variants: " + maxRBetweenVariants);
				System.out.println(" - Number of threads to use: " + numberOfThreadsToUse);
				System.out.println(" - Gene info file: " + geneInfoFile.getAbsolutePath());
				break;
		}

	}

	public String[] getGenotypeBasePath() {
		return genotypeBasePath;
	}

	public RandomAccessGenotypeDataReaderFormats getGenotypeType() {
		return genotypeType;
	}

	public File getOutputFile() {
		return outputFile;
	}

	public String getGwasZscoreMatrixPath() {
		return gwasZscoreMatrixPath;
	}

	public int getNumberOfPermutations() {
		return numberOfPermutations;
	}

	public int getWindowExtend() {
		return windowExtend;
	}

	public double getMaxRBetweenVariants() {
		return maxRBetweenVariants;
	}

	public static int getNumberOfThreadsToUse() {
		return numberOfThreadsToUse;
	}

	public File getGeneInfoFile() {
		return geneInfoFile;
	}

	public Depict2Mode getMode() {
		return mode;
	}

}
