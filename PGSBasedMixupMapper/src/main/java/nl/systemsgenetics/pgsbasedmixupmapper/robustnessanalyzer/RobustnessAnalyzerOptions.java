package nl.systemsgenetics.pgsbasedmixupmapper.robustnessanalyzer;

import nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapperOptions;
import org.apache.commons.cli.*;

import java.util.ArrayList;

public class RobustnessAnalyzerOptions extends PGSBasedMixupMapperOptions {

    private static Options OPTIONS = getOptions();
    private final ArrayList<Double> mixUpPercentages;

    static {
        OPTIONS.addOption(getMixUpPercentagesOption());
    }

    public RobustnessAnalyzerOptions(String... args) throws ParseException {
        super(args);
        mixUpPercentages = parseMixUpPercentages(getCommandLine());
    }

    private ArrayList<Double> parseMixUpPercentages(CommandLine commandLine) throws ParseException {
        ArrayList<Double> parsedMixUpPercentages = new ArrayList<>();

        String[] mixUpPercentages = commandLine.getOptionValues("mixUpPercentages");
        for (String mixUpPercentage : mixUpPercentages) {
            // Parse the percentage to a double
            try {
                double percentage = Double.parseDouble(mixUpPercentage);
                if (percentage > 100) {
                    throw new ParseException(String.format(
                            "Error parsing -pc / --mixUpPercentages: \"%f\" is greater than 100.",
                            percentage));
                }
                parsedMixUpPercentages.add(percentage);
                // If the supposed percentage cannot be parsed to a double, a NumberFormatException is thrown,
                // catch this and throw a parse exception.
            } catch (NumberFormatException e) {
                throw new ParseException(String.format(
                        "Error parsing -pc / --mixUpPercentages: \"%s\" could not be parsed to a double",
                        mixUpPercentage));
            }
        }
        return parsedMixUpPercentages;
    }

    private static Option getMixUpPercentagesOption() {
        OptionBuilder.withArgName("float");
        OptionBuilder.hasArgs();
        OptionBuilder.withDescription("Percentages of mix-ups to introduce");
        OptionBuilder.withLongOpt("mixUpPercentages");
        OptionBuilder.isRequired();
        return OptionBuilder.create("pc");
    }

    public ArrayList<Double> getMixUpPercentages() {
        return mixUpPercentages;
    }
}
