package nl.systemsgenetics.gwassummarystatistics;

import nl.systemsgenetics.gwassummarystatistics.writer.GwasSummaryStatisticsWriter;
import nl.systemsgenetics.gwassummarystatistics.writer.OutputFormat;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

public class VcfSummaryStatisticsReformatter {

    private static final Logger LOGGER;
    private static final Options OPTIONS;

    static {
        LOGGER = Logger.getLogger(VcfSummaryStatisticsReformatter.class);

        OPTIONS = new Options();

        OPTIONS.addOption(Option.builder("i")
                .longOpt("inputVcf")
                .hasArgs()
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
            ReadOnlyGwasSummaryStatistics summaryStatistics =
                    new ReadOnlyGwasSummaryStatistics(new VcfGwasSummaryStatistics(
                            new File("/Users/cawarmerdam/Documents/projects/pgs_based_mixup_correction/data/summarystatistics/vcf/ukb-d-1747_2.vcf.gz"),
                            new File("/Users/cawarmerdam/Documents/projects/pgs_based_mixup_correction/data/summarystatistics/vcf/ukb-d-1747_2.vcf.gz.tbi")
                    ), "ukb-d-1747_2");

            GwasSummaryStatisticsWriter gwasSummaryStatisticsWriter = new GwasSummaryStatisticsWriter(
                    summaryStatistics, OutputFormat.PRS_CS);
            gwasSummaryStatisticsWriter.save("/Users/cawarmerdam/Documents/projects/pgs_based_mixup_correction/data/summarystatistics/PRScs/ukb-d-1747_2");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
