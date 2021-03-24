package nl.systemsgenetics.gwassummarystatistics;

import com.google.common.collect.ImmutableMap;
import nl.systemsgenetics.gwassummarystatistics.writer.GwasSummaryStatisticsWriter;
import nl.systemsgenetics.gwassummarystatistics.writer.OutputFormat;
import org.apache.commons.cli.*;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class VcfSummaryStatisticsReformatter {

    private static final Logger LOGGER;
    private static final Options OPTIONS;

    static {
        LOGGER = Logger.getLogger(VcfSummaryStatisticsReformatter.class);

        OPTIONS = new Options();

        OPTIONS.addOption(Option.builder("i")
                .longOpt("inputVcf")
                .hasArg()
                .desc("The path of the VCF summary statistics file to convert.")
                .required().build());

        OPTIONS.addOption(Option.builder("o")
                .longOpt("outputPath")
                .hasArg()
                .desc("The output bash path")
                .required().build());

        OPTIONS.addOption(Option.builder("O")
                .longOpt("OutputFormat")
                .hasArgs()
                .desc("The output format").build());
    }

    public static void main(String[] args) {

        try {
            CommandLineParser parser = new DefaultParser();
            final CommandLine commandLine = parser.parse(OPTIONS, args, true);

            final String inputVcfPath = commandLine.getOptionValue('i');
            File vcfFile;

            if (inputVcfPath.endsWith(".vcf.gz")) {
                vcfFile = new File(inputVcfPath);
            } else {
                vcfFile = new File(inputVcfPath + ".vcf.gz");
            }

            if (!vcfFile.exists()) {
                throw new ParseException(String.format(
                        "Vcf file '%s' does not exist", vcfFile
                ));
            }

            File vcfTabixFile = new File(vcfFile.getPath() + ".tbi");

            if (!vcfTabixFile.exists()) {
                throw new ParseException(String.format(
                        "Tabix file '%s' does not exist", vcfFile
                ));
            }

            final File outputPath = new File(commandLine.getOptionValue('o'));

            final OutputFormat outputFormat = parseOutputFormatOption(commandLine.getOptionValues('O'));

            GwasSummaryStatistics summaryStatistics = new ReadOnlyGwasSummaryStatistics(
                    new VcfGwasSummaryStatistics(vcfFile, vcfTabixFile)
            );

            GwasSummaryStatisticsWriter gwasSummaryStatisticsWriter = new GwasSummaryStatisticsWriter(
                    summaryStatistics, outputFormat);

            if (outputPath.mkdirs()) {
                System.out.println("Created base directory");
            }

            gwasSummaryStatisticsWriter.save(outputPath);

        } catch (ParseException e) {
            LOGGER.fatal("Invalid command line arguments: ");
            LOGGER.fatal(e.getMessage());
            System.err.println();
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", OPTIONS);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static OutputFormat parseOutputFormatOption(String[] format) throws ParseException {
        OutputFormat outputFormat = null;
        if (format.length == 1) {
            String givenFormatAsString = format[0];
            if (!OutputFormat.getPredefinedOutputFormats().contains(givenFormatAsString)) {
                throw new ParseException(String.format(
                        "The given format '%s' is not recognized.%n" +
                                "Use '-O/--outputFormat CUSTOM <custom format>' for a custom output format.",
                        givenFormatAsString));
            }
            outputFormat = OutputFormat.getPredefinedOutputFormat(givenFormatAsString);
        } else if (format.length == 2) {
            if ("CUSTOM".equals(format[0])) {
                outputFormat = parseCustomOutputFormat(format[1]);
            }
        } else {
            throw new ParseException(String.format("Found %d arguments while 1-2 were expected.", format.length));
        }

        return outputFormat;
    }

    private static OutputFormat parseCustomOutputFormat(String customOutputFormatAsString) {
        // Find all substrings matching % plus one of the columns
        List<OutputFormat.Column> enumValues = Stream.of(OutputFormat.Column.values())
                .sorted(Comparator.comparing(e -> e.name().length()))
                .collect(Collectors.toList());
        List<String> columns = new ArrayList<>();

        String[] split = customOutputFormatAsString.split("%");
        for (String formatUnitAsString : split) {
            for (OutputFormat.Column column : enumValues) {
                if (formatUnitAsString.startsWith(column.name())) {
                    columns.add(column.name());
                }
            }
        }

        return new OutputFormat(columns);
    }
}