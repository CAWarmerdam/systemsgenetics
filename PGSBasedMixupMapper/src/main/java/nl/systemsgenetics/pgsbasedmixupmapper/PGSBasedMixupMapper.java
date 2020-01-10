package nl.systemsgenetics.pgsbasedmixupmapper;

import org.apache.commons.cli.ParseException;
import org.apache.log4j.*;

import java.io.*;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 *
 * @author Patrick Deelen
 */
public class PGSBasedMixupMapper {

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

		if (args.length == 0) {
			PGSBasedMixupMapperOptions.printHelp();
			return;
		}

		try {
			options = new PGSBasedMixupMapperOptions(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			PGSBasedMixupMapperOptions.printHelp();
			return;
		}

		if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
			if (!options.getLogFile().getParentFile().mkdirs()) {
				System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
				System.exit(1);
			}
		}

		if (new File(options.getOutputBasePath()).isDirectory()) {
			System.err.println("Specified output path is a directory. Please include a prefix for the output files.");
			return;
		}
	}
}