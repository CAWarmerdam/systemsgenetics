package nl.systemsgenetics.gwassummarystatistics.writer;

import com.opencsv.CSVWriter;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.ReadOnlyGwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;

public class GwasSummaryStatisticsWriter {

    private GwasSummaryStatistics summaryStatistics;
    private final OutputFormat outputFormat;

    public GwasSummaryStatisticsWriter(
            ReadOnlyGwasSummaryStatistics summaryStatistics,
            OutputFormat outputFormat) {

        this.summaryStatistics = summaryStatistics;
        this.outputFormat = outputFormat;
    }

    public void readConfig() {}

    public void save(String pathPrefix) throws IOException {
        // Format the output path.
        String outputPath = String.format("%s_%s.tsv",
                pathPrefix, summaryStatistics.getGwasId().replace(" ", "_"));

        // Create a CSV writer without quotation characters and with tabs as separators to mimic
        // the legacy gwas summary statistics files.
        CSVWriter writer = new CSVWriter(new FileWriter(outputPath),
                '\t', CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        // Get the columns to write
        Map<String, Function<EffectAllele, String>> columns = outputFormat.getColumns();
        String[] rowArray = new String[columns.size()];

        // Write a header
        writer.writeNext(columns.keySet().toArray(rowArray));

        // Write risk entries
        for (EffectAllele effectAllele : summaryStatistics) {
            List<String> row = new ArrayList<>();

            for (String column : columns.keySet()) {
                row.add(columns.get(column).apply(effectAllele));
            }

            writer.writeNext(row.toArray(rowArray));
        }

        writer.close();
    }
}
