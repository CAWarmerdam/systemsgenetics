package nl.systemsgenetics.depict2;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.ResourceBundle;
import org.apache.commons.cli.ParseException;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;

/**
 *
 * @author Patrick Deelen
 */
public class Depict2 {

	private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |                DEPICT2                |\n"
			+ "  |                                       |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	/**
	 * @param args the command line arguments
	 * @throws java.lang.InterruptedException
	 */
	public static void main(String[] args) throws InterruptedException {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		Date currentDataTime = new Date();
		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		Thread.sleep(25); //Allows flush to complete

		Depict2Options options;

		if (args.length == 0) {
			Depict2Options.printHelp();
			return;
		}

		try {
			options = new Depict2Options(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			Depict2Options.printHelp();
			return;
		}

		options.printOptions();
		
		readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath()))

		if (options.getOutputFile().getParentFile() != null && !options.getOutputFile().getParentFile().isDirectory()) {
			if (!options.getOutputFile().getParentFile().mkdirs()) {
				System.err.println("Failed to create output folder: " + options.getOutputFile().getParent());
				System.exit(1);
			}
		}

		RandomAccessGenotypeData referenceGenotypeData = loadGenotypes(options);

		//DoubleMatrixDataset<String, String>genePvalues = CalculateGenePvalues.calculatorGenePvalues(variantPhenotypeZscoreMatrix, new GenotypeCovarianceGenotypes(referenceGenotypeData), genes, options.getWindowExtend(), options.getMaxRBetweenVariants(), options.getNumberOfPermutations());
	}

	private static RandomAccessGenotypeData loadGenotypes(Depict2Options options) {
		final RandomAccessGenotypeData referenceGenotypeData;
		try {
			referenceGenotypeData = options.getGenotypeType().createFilteredGenotypeData(options.getGenotypeBasePath(), 10000, null, null, null, 0.34f);
		} catch (TabixFileNotFoundException e) {
			System.err.println("Tabix file not found for input data at: " + e.getPath() + "\n"
					+ "Please see README on how to create a tabix file");
			System.exit(1);
			return null;
		} catch (IOException e) {
			System.err.println("Error accessing input data: " + e.getMessage());
			System.exit(1);
			return null;
		} catch (IncompatibleMultiPartGenotypeDataException e) {
			System.err.println("Error combining the impute genotype data files: " + e.getMessage());
			System.exit(1);
			return null;
		} catch (GenotypeDataException e) {
			System.err.println("Error reading input data: " + e.getMessage());
			System.exit(1);
			return null;
		}
		return referenceGenotypeData;
	}

	protected static final List<String> readMatrixAnnotations(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).build();

		ArrayList<String> identifiers = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			identifiers.add(nextLine[0]);

		}

		return identifiers;

	}

}
