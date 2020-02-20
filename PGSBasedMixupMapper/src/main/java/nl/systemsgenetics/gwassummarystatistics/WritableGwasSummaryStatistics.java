package nl.systemsgenetics.gwassummarystatistics;

import com.opencsv.CSVWriter;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.apache.log4j.Logger;

import java.io.FileWriter;
import java.io.IOException;

public abstract class WritableGwasSummaryStatistics implements GwasSummaryStatistics {

    private static final Logger LOGGER = Logger.getLogger(WritableGwasSummaryStatistics.class);
    private String gwasId;

    public WritableGwasSummaryStatistics(String gwasId) {
        this.gwasId = gwasId;
    }

    @Override
    public String getGwasId() {
        return gwasId;
    }

    public void save(String pathPrefix) throws IOException {
        // Format the output path.
        String outputPath = String.format("%s_%s.tsv",
                pathPrefix, gwasId.replace(" ", "_"));

        // Create a CSV writer without quotation characters and with tabs as separators to mimic
        // the legacy gwas summary statistics files.
        CSVWriter writer = new CSVWriter(new FileWriter(outputPath),
                '\t', CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        // Write a header
        writer.writeNext(new String[]{"variant", "AssessedAllele", "ES", "pvalue"});

        // Write risk entries
        for (EffectAllele effectAllele : this) {
            writer.writeNext(new String[]{
                    effectAllele.getPrimaryVariantId(),
                    String.valueOf(effectAllele.getAllele()),
                    String.valueOf(effectAllele.getEffectSize()),
                    String.valueOf(effectAllele.getPValue())});
        }

        writer.close();
    }
}
