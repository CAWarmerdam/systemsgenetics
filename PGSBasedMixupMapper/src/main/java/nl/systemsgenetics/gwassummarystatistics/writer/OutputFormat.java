package nl.systemsgenetics.gwassummarystatistics.writer;

import com.google.common.collect.ImmutableMap;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;

import java.util.*;
import java.util.function.Function;

public class OutputFormat {

    private static final Map<Column, Function<EffectAllele, String>> DEFAULT_FUNCTIONS = createColumnMap();
    private static final Map<String, OutputFormat> predefinedOutputFormats = createPredefinedOutputFormats();

    public static OutputFormat getPredefinedOutputFormat(String key) {
        return predefinedOutputFormats.get(key);
    }

    Map<String, Function<EffectAllele, String>> getColumns() {
        return columns;
    }

    private Map<String, Function<EffectAllele, String>> columns;

    private static Map<String, Function<EffectAllele, String>> mapColumns(Map<Column, String> columnNames) {
        Map<String, Function<EffectAllele, String>> columns = new LinkedHashMap<>();

        for (Column column : columnNames.keySet()) {
            columns.put(columnNames.get(column), DEFAULT_FUNCTIONS.get(column));
        }
        return columns;
    }

    private OutputFormat(Map<String, Function<EffectAllele, String>> columns) {
        this.columns = columns;
    }

    private static Map<String, OutputFormat> createPredefinedOutputFormats() {
        Map<String, OutputFormat> outputFormats = new LinkedHashMap<>();

        outputFormats.put("PRScs", new OutputFormat(mapColumns(ImmutableMap.<Column, String>builder()
                .put(Column.ID, "SNP")
                .put(Column.EA, "A1")
                .put(Column.NEA, "A2")
                .put(Column.ES, "BETA")
                .put(Column.P, "P").build())));
        return outputFormats;
    }

    public OutputFormat(List<String> columnNames) {
        LinkedHashMap<Column, String> columns = new LinkedHashMap<>();

        for (String columnName : columnNames) {
            columns.put(Column.valueOf(columnName), columnName);
        }
        this.columns = mapColumns(columns);
    }

    public static Set<String> getPredefinedOutputFormats() {
        return predefinedOutputFormats.keySet();
    }

    private static HashMap<Column, Function<EffectAllele, String>> createColumnMap() {

        HashMap<Column, Function<EffectAllele, String>> columnMap = new HashMap<>();

        columnMap.put(Column.ID, EffectAllele::getPrimaryVariantId);
        columnMap.put(Column.CHROM, EffectAllele::getSequenceName);
        columnMap.put(Column.POS, (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getStartPos()));

        columnMap.put(Column.EA, (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getAllele()));

        columnMap.put(Column.NEA, (EffectAllele effectAllele) ->
                effectAllele.getNonEffectAllele().toString());

        columnMap.put(Column.EAF, (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getAlleleFrequency()));

        columnMap.put(Column.P, (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getPValue()));

        columnMap.put(Column.LP, (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getLogTransformedPValue()));

        columnMap.put(Column.ES, (EffectAllele effectAllele) ->
                String.valueOf(effectAllele.getEffectSize()));

        return columnMap;
    }

    public enum Column {
        ID, EA, NEA, ES, P, CHROM, POS, LP, EAF
    }
}
