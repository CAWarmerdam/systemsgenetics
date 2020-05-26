package nl.systemsgenetics.gwassummarystatistics;

import com.opencsv.CSVWriter;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import org.apache.log4j.Logger;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public abstract class WritableGwasSummaryStatistics implements GwasSummaryStatistics {

    private HashMap<String, Function<EffectAllele, String>> columnMap = createColumnMap();
    private String[] requestedColumns;
    private CharSequence alleleDelimiter = ",";

    private HashMap<String, Function<EffectAllele, String>> createColumnMap() {

        HashMap<String, Function<EffectAllele, String>> columnMap = new HashMap<>();

        columnMap.put("ID", EffectAllele::getPrimaryVariantId);
        columnMap.put("CHR", EffectAllele::getSequenceName);
        columnMap.put("POS", (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getStartPos()));

        columnMap.put("EA", (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getAllele()));

        columnMap.put("NEA", (EffectAllele effectAllele) ->
                effectAllele.getNonEffectAlleles().getAlleles().stream()
                .map(object -> Objects.toString(object, null))
                .collect(Collectors.joining(alleleDelimiter)));

        columnMap.put("EAF", (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getAlleleFrequency()));

        columnMap.put("PVAL", (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getPValue()));

        columnMap.put("LOG10PVAL", (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getLogTransformedPValue()));

        columnMap.put("ES", (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getEffectSize()));

        return columnMap;
    }

    public void readConfig() {}

    public void save(String pathPrefix) throws IOException {
        // Format the output path.
        String outputPath = String.format("%s_%s.tsv",
                pathPrefix, getGwasId().replace(" ", "_"));

        // Create a CSV writer without quotation characters and with tabs as separators to mimic
        // the legacy gwas summary statistics files.
        CSVWriter writer = new CSVWriter(new FileWriter(outputPath),
                '\t', CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        // Write a header
        writer.writeNext(requestedColumns);

        // Write risk entries
        for (EffectAllele effectAllele : this) {
            String[] rowArray = new String[requestedColumns.length];
            List<String> row = new ArrayList<>();

            for (String column : requestedColumns) {
                row.add(columnMap.get(column).apply(effectAllele));
            }

            writer.writeNext(row.toArray(rowArray));
        }

        writer.close();
    }
}
