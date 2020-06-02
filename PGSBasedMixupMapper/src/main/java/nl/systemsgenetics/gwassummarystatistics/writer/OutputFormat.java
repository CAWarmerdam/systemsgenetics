package nl.systemsgenetics.gwassummarystatistics.writer;

import com.google.common.collect.ImmutableMap;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatisticsException;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import org.molgenis.genotype.Allele;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

public enum OutputFormat {
    PRS_CS(mapColumns(ImmutableMap.<Column, String>builder()
            .put(Column.ID, "SNP")
            .put(Column.EA, "A1")
            .put(Column.NEA, "A2")
            .put(Column.ES, "BETA")
            .put(Column.P, "P").build()));

    private static final Map<Column, Function<EffectAllele, String>> DEFAULT_FUNCTIONS = createColumnMap();

    public Map<String, Function<EffectAllele, String>> getColumns() {
        return columns;
    }

    private Map<String, Function<EffectAllele, String>> columns;

    static Map<String, Function<EffectAllele, String>> mapColumns(Map<Column, String> columnNames) {
        Map<String, Function<EffectAllele, String>> columns = new HashMap<>();

        for (Column column : columnNames.keySet()) {
            columns.put(columnNames.get(column), DEFAULT_FUNCTIONS.get(column));
        }
        return columns;
    }

    OutputFormat(Map<String, Function<EffectAllele, String>> columns) {
        this.columns = columns;
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
