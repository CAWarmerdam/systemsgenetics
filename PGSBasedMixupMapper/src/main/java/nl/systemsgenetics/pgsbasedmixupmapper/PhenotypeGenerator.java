package nl.systemsgenetics.pgsbasedmixupmapper;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.jet.math.tdouble.DoubleFunctions;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatistics;
import nl.systemsgenetics.polygenicscorecalculator.*;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variantFilter.VariantFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import static nl.systemsgenetics.pgsbasedmixupmapper.PGSBasedMixupMapper.*;

public class PhenotypeGenerator {
    private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static final char CSV_DELIMITER = ',';
    private static final org.apache.log4j.Logger LOGGER = Logger.getLogger(PhenotypeGenerator.class);


    public static void main(String[] args) throws InterruptedException, IOException {

        System.out.println();
        System.out.println("          --- Version: " + " ---");
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

        // Print a help whenever the arguments have a length of zero.
        if (args.length == 0) {
            PGSBasedMixupMapperOptions.printHelp();
            return;
        }

        // Parse the arguments list
        options = getPgsBasedMixupMapperOptions(args);

        // Check if the output directory is a prefix for output files. May not be a directory
        if (options.getOutputBasePath().isDirectory()) {
            System.err.println("Specified output path is a directory. Please include a prefix for the output files.");
            return;
        }

        // Create a logger (set the correct file for output)
        createLogger(startDateTime, options);

        // Print the options
        options.printOptions();

        // Load genotype to phenotype sample coupling map
        Map<String, String> genotypeToPhenotypeSampleCoupling = loadGenotypeToPhenotypeSampleCoupling(
                options.getGenotypeToPhenotypeSampleCouplingFile(), CSV_DELIMITER);

        // Load Genotype and trait data
        RandomAccessGenotypeData genotypeData = loadGenotypeData(options, genotypeToPhenotypeSampleCoupling.keySet());
        VariantFilter variantFilter = getVariantFilter(genotypeData);

        // Load GWAS summary statistics
        Map<String, String> gwasPhenotypeCoupling = loadGwasSummaryStatisticsPhenotypeCouplings(
                options.getGwasSummaryStatisticsPhenotypeCouplingFile(), CSV_DELIMITER);
        Map<String, GwasSummaryStatistics> gwasSummaryStatisticsMap = loadGwasSummaryStatisticsMap(
                options.getGwasSummaryStatisticsPath(), gwasPhenotypeCoupling, variantFilter);
//        System.out.println("Height");
//        DoubleMatrix1D height = generatePhenotype(genotypeData, gwasSummaryStatisticsMap.get("Height"));
        System.out.println("HDL");
        DoubleMatrix1D hdl = generatePhenotype(genotypeData, gwasSummaryStatisticsMap.get("HDL cholesterol"), options);
        System.out.println("LDL");
        DoubleMatrix1D ldl = generatePhenotype(genotypeData, gwasSummaryStatisticsMap.get("LDL cholesterol"), options);
//        DoubleMatrix1D bmi = generatePhenotype(genotypeData, gwasSummaryStatisticsMap.get("BMI"));

        DoubleMatrixDataset<String, String> phenotypeMatrix = new DoubleMatrixDataset<>(Arrays.asList(genotypeData.getSampleNames()),
                new ArrayList<>(Arrays.asList("LDL cholesterol", "HDL cholesterol")));

//        phenotypeMatrix.viewCol("Height").assign(height);
        phenotypeMatrix.viewCol("LDL cholesterol").assign(ldl);
        phenotypeMatrix.viewCol("HDL cholesterol").assign(hdl);
//        phenotypeMatrix.viewCol("LDL cholesterol").assign(ldl);
//        phenotypeMatrix.viewCol("BMI").assign(bmi);
        phenotypeMatrix.save(options.getOutputBasePath() + "_simphenotypes.tsv");
    }

    private static DoubleMatrix1D generatePhenotype(RandomAccessGenotypeData genotypeData,
                                                    GwasSummaryStatistics vcfGwasSummaryStatistics,
                                                    PGSBasedMixupMapperOptions options) {
        DoubleMatrixDataset<String, String> scores = calculatePolyGenicScores(
                vcfGwasSummaryStatistics, genotypeData, options);
        DoubleMatrix1D row = scores.getRow(0).assign(DoubleFunctions.mult(4));
//        double[] errorArray = getErrorArray(genotypeData.getSamples().size());
//        row.assign(new DenseDoubleMatrix1D(errorArray), DoubleFunctions.plus);
        return row;
    }

    private static double[] getErrorArray(int size) {
        Random fRandom = new Random();
        double aMean = 0.0f;
        double aVariance = 0.1f;
        double[] error = new double[size];
        for (int i = 0; i < error.length; i++) {
            error[i] = aMean + fRandom.nextGaussian() * aVariance;
        }
        return error;
    }

    private static DoubleMatrixDataset<String, String> calculatePolyGenicScores(
            GwasSummaryStatistics summaryStatistics,
            RandomAccessGenotypeData genotypeData, PGSBasedMixupMapperOptions options) {

        int[] windowSize = new int[]{options.getWindowSize().get(0)};
        LOGGER.info(String.format("Using window size of %s for phenotype generation", Arrays.toString(windowSize)));

        boolean sumRisks = false;
        double v = options.getpValueThresholds().get(options.getpValueThresholds().size() - 1);
        LOGGER.info(String.format("Using a pvalue threshold of %f for phenotype generation", v));
        List<Double> pValThres = Collections.singletonList(v);
        double rSquare = options.getrSquared();
        LOGGER.info(String.format("Using r2 size of %f for phenotype generation", rSquare));
        boolean unweighted = false;
        LOGGER.info(String.format("Using gr to exclude '%s' for phenotype generation",
                Arrays.toString(options.getGenomicRangesToExclude())));

        System.out.println(genotypeData.hashCode());

        DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = null;

        // Initialize an Simple polygenic score calculator
        SimplePolyGenicScoreCalculator polyGenicScoreCalculator = new SimplePolyGenicScoreCalculator(
                genotypeData, // The genotype data to calculate an LD matrix in.
                options.getWindowSize(), // Get the window size in number of base pairs,
                pValThres,
                options.getrSquared(),
                false,
                options.getGenomicRangesToExclude());

        geneticRiskScoreMatrix = polyGenicScoreCalculator.calculate(summaryStatistics);
//        if (windowSize.length == 1) {
//            geneticRiskScoreMatrix = SimplePolyGenicScoreCalculator.calculate(genotypeData, rSquare, windowSize[0], pValThres, sumRisks);
////			} else if (windowSize.length == 2) {
////				DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = SimplePolyGenicScoreCalculator.calculateTwoStages(genotypeData, risks, outputFolder, rSquare, windowSize, debugMode, pValThres, sumRisks);
//        } else {
//            System.out.println("More than two window-sizes is not supported.");
//            System.exit(0);
//        }

        return geneticRiskScoreMatrix;
    }
}
