/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.polygenicscorecalculator;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatisticsException;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.effectAllele.RiskEntry;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleFilterableGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilterableGenotypeDataDecorator;
import org.molgenis.genotype.sampleFilter.SampleFilteredReadOnlyGeneticVariant;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.Map.Entry;

import static org.molgenis.genotype.util.LdCalculator.calculateRsquare;

/**
 * @author MarcJan
 */
public class SimplePolygenicScoreCalculator {

    private static final Logger LOGGER = Logger.getLogger(SimplePolygenicScoreCalculator.class);
    private static final String[] chrOrder = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};
    private final RandomAccessGenotypeData genotypeData;
    private final List<Integer> windowSizeList;
    private final List<Double> pValueThresholds;
    private final double rSquared;
    private final boolean sumRisk;
    private final String[] genomicRangesToExclude;

    public SimplePolygenicScoreCalculator(RandomAccessGenotypeData genotypeData,
                                          List<Integer> windowSizeList,
                                          List<Double> pValueThresholds,
                                          double rSquared, boolean sumRisk,
                                          String[] genomicRangesToExclude) {

        this.genotypeData = genotypeData;
        this.windowSizeList = windowSizeList;
        this.pValueThresholds = pValueThresholds;
        this.rSquared = rSquared;
        this.sumRisk = sumRisk;
        this.genomicRangesToExclude = genomicRangesToExclude;
    }

    public DoubleMatrixDataset<String, String> calculate(
            GwasSummaryStatistics summaryStatistics) {

        double[] pValThres = this.getpValueThresholds()
                .stream().mapToDouble(Double::doubleValue).toArray();
        boolean unweighted = false;
        String[] genomicRangesToExclude = this.getGenomicRangesToExclude();
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks =
                riskEntries(summaryStatistics, genotypeData,
                        pValThres, genomicRangesToExclude,
                        unweighted);

        DoubleMatrixDataset<String, String> scores = initializePolygenicScoreMatrix(genotypeData, risks);

        if (windowSizeList.size() == 1) {
            return calculate(scores, risks, windowSizeList.get(0));
        } else if (windowSizeList.size() == 2) {
            return calculateTwoStages(scores, risks, windowSizeList.stream().mapToInt(i->i).toArray());
        }
        throw new UnsupportedOperationException("More than two window sizes not supported");
    }

    public DoubleMatrixDataset<String, String> calculate(
            GwasSummaryStatistics summaryStatistics,
            SampleFilter referenceSampleFilter,
            SampleFilter responseSampleFilter) {

        double[] pValThres = this.getpValueThresholds()
                .stream().mapToDouble(Double::doubleValue).toArray();
        boolean unweighted = false;
        String[] genomicRangesToExclude = this.getGenomicRangesToExclude();
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks =
                riskEntries(summaryStatistics, genotypeData,
                        pValThres, genomicRangesToExclude,
                        unweighted);

        SampleFilterableGenotypeDataDecorator referenceGenotypeData =
                new SampleFilterableGenotypeDataDecorator(genotypeData, referenceSampleFilter);
        SampleFilterableGenotypeDataDecorator responsePhenotypeData =
                new SampleFilterableGenotypeDataDecorator(genotypeData, responseSampleFilter);

        DoubleMatrixDataset<String, String> scores = initializePolygenicScoreMatrix(responsePhenotypeData, risks);

        if (windowSizeList.size() == 1) {
            throw new UnsupportedOperationException();
        } else if (windowSizeList.size() == 2) {
            return calculateTwoStages(scores, risks, windowSizeList.stream().mapToInt(i->i).toArray(),
                    referenceGenotypeData, responsePhenotypeData);
        }
        throw new UnsupportedOperationException("More than two window sizes not supported");
    }

    public DoubleMatrixDataset<String, String> calculate(
            DoubleMatrixDataset<String, String> scores,
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks,
            double windowSize, SampleFilterableGenotypeData referenceGenotypeData,
            SampleFilterableGenotypeData responseGenotypeData) {

//        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        // Loop through the chromosomes
        for (int counter = 0; counter < chrOrder.length; counter++) {
            // Loop through the risk entries for every file
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                // Loop through the p-values
                HashSet<String> excludeList = new HashSet<String>();
                for (double pVal : pValueThresholds) {
                    // Generate a _P<pvalue> key
                    String key = "_P" + pVal;
                    THashMap<String, ArrayList<RiskEntry>> riskScorePheno2 = riskScorePheno.getValue().get(key);
                    // Generate a <filename>_P<pvalue> string
                    String NameOfEntry = key;
                    int rowNr = scores.getHashRows().get(NameOfEntry);
                    if (LOGGER.isDebugEnabled()) {
                        System.out.println(NameOfEntry);

                        LOGGER.debug("SNPs used for GRS calculation:\n");
                    }
                    int nrSNPs = 0;

//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                    if (riskScorePheno2.containsKey(chrOrder[counter])) {

                        // Get the riskentries for this sequence / chromosome
                        ArrayList<RiskEntry> valueE2 = riskScorePheno2.get(chrOrder[counter]);

                        int nrSNPsThisChr = valueE2.size();
                        boolean[] excludeSNPs = new boolean[nrSNPsThisChr];

                        // Get the original entries back, so we are sure we dont need to do to many look ups.
                        if (excludeList.size() > 0) {
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (excludeList.contains(valueE2.get(snp).getPrimaryVariantId())) {
                                    excludeSNPs[snp] = true;
                                }
                            }
                        }

                        // Actual scoring.

                        // Loop through the individual risk scores (these should be sorted)
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
                                //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant originalVariant = genotypeData.getSnpVariantByPos(riskE.getSequenceName(), riskE.getStartPos());
                                GeneticVariant referenceVariant = new SampleFilteredReadOnlyGeneticVariant(originalVariant, referenceGenotypeData);

                                //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                if (referenceVariant.getCallRate() < 0.75) {
                                    excludeSNPs[snp] = true;
                                    excludeList.add(riskE.getPrimaryVariantId());
                                    continue;
                                }

                                if (LOGGER.isDebugEnabled()) {
                                    LOGGER.debug(riskE.toString());
                                }

                                // Get the reference allele from this variant 1
                                Allele var1RefAllele = referenceVariant.getRefAllele();
                                Allele alternativeAllele;
                                if (var1RefAllele == null) {
                                    var1RefAllele = referenceVariant.getAlternativeAlleles().get(0);
                                    alternativeAllele = referenceVariant.getAlternativeAlleles().get(1);
                                } else {
                                    alternativeAllele = referenceVariant.getAlternativeAlleles().get(0);
                                }

                                if (!String.valueOf(riskE.getAllele()).equals(String.valueOf(alternativeAllele))) {
                                    System.out.printf("Allele match '%s' | '%s' (effect allele | genotype data)%n",
                                            riskE.getAllele(), alternativeAllele);
                                }

                                double or = riskE.getEffectSize();
//                                System.out.println("or = " + or);
                                boolean riskCodedAsTwo;
                                if (sumRisk && or < 0) {
                                    or = or * -1;
                                    riskCodedAsTwo = false;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = true;
                                    }
                                } else {
                                    riskCodedAsTwo = true;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = false;
                                    }
                                }
//                                System.out.println("riskCodedAsTwo = " + riskCodedAsTwo);
//                                    StringBuilder Genos = new StringBuilder();
//                                    StringBuilder Genos1 = new StringBuilder();
//                                    StringBuilder Genos2 = new StringBuilder();
//                                    StringBuilder Genos3 = new StringBuilder();
//                                System.out.println("var1.getSampleCalledDosages() = " + Arrays.toString(var1.getSampleCalledDosages()));
//                                System.out.println("or = " + or);
                                GeneticVariant responseVariant =
                                        new SampleFilteredReadOnlyGeneticVariant(originalVariant, responseGenotypeData);

                                byte[] sampleCalledDosages = responseVariant.getSampleCalledDosages();
                                for (int sample = 0; sample < sampleCalledDosages.length; sample++) {
                                    if (sampleCalledDosages[sample] != -1) {
                                        if (riskCodedAsTwo) {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * sampleCalledDosages[sample])));
//                                                Genos2.append(or * var1.getSampleCalledDosages()[sample]).append(",");
                                        } else {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * Math.abs(sampleCalledDosages[sample] - 2))));
//                                                Genos2.append(or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                        }
//                                            Genos.append(var1.getSampleCalledDosages()[sample]).append(",");
//                                            Genos1.append(var1.getSampleVariants().get(sample).toString());Genos1.append(",");
//                                            Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                    }
                                }
                                //I have now removed the conversion from 0 to 2 and vice versa from the calculations. Why is this needed?
//                                            System.out.println("");
//                                            System.out.println("SNP: "+riskE.getRsName());
//                                            System.out.println("Allele in data: "+var1.getRefAllele().toString());
//                                            System.out.println("Allele: "+riskE.getAllele());
//                                            System.out.println("Diction: "+direction);
//                                            System.out.println("Or: "+or);
//                                            System.out.println("Genotypes: "+Genos1.toString());
//                                            System.out.println("Dosages: "+Genos.toString());
//                                            System.out.println("Scores: "+Genos2.toString());
//                                            System.out.println("ScoresT: "+Genos3.toString());
//                                            System.out.println("");

                                nrSNPs++;

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getStartPos() - riskE.getStartPos()) <= windowSize) {
                                            GeneticVariant var2 = referenceGenotypeData.getSnpVariantByPos(
                                                    riskE2.getSequenceName(), riskE2.getStartPos());
                                            if (var2.getCallRate() < 0.75) {
                                                excludeSNPs[t] = true;
                                                excludeList.add(riskE2.getPrimaryVariantId());
                                                continue;
                                            }
                                            try {
                                                if (calculateRsquare(referenceVariant, var2) >= rSquared) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getPrimaryVariantId());
                                                }
                                            } catch (LdCalculatorException ex) {
                                                LOGGER.error(ex);
                                            }

                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (LOGGER.isDebugEnabled()) {
                        LOGGER.debug("Total SNPs used: " + nrSNPs);
                    }
                }
//                p.iterate();
            }
        }
//        p.close();

        return scores;
    }


    public DoubleMatrixDataset<String, String> calculate(
            DoubleMatrixDataset<String, String> scores,
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks,
            double windowSize) {

        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        // Loop through the chromosomes
        for (int counter = 0; counter < chrOrder.length; counter++) {
            // Loop through the risk entries for every file
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                // Loop through the p-values
                HashSet<String> excludeList = new HashSet<String>();
                for (double pVal : pValueThresholds) {
                    // Generate a _P<pvalue> key
                    String key = "_P" + pVal;
                    THashMap<String, ArrayList<RiskEntry>> riskScorePheno2 = riskScorePheno.getValue().get(key);
                    // Generate a <filename>_P<pvalue> string
                    String NameOfEntry = key;
                    int rowNr = scores.getHashRows().get(NameOfEntry);
                    if (LOGGER.isDebugEnabled()) {
                        System.out.println(NameOfEntry);

                        LOGGER.debug("SNPs used for GRS calculation:\n");
                    }
                    int nrSNPs = 0;

//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                    if (riskScorePheno2.containsKey(chrOrder[counter])) {

                        // Get the riskentries for this sequence / chromosome
                        ArrayList<RiskEntry> valueE2 = riskScorePheno2.get(chrOrder[counter]);

                        int nrSNPsThisChr = valueE2.size();
                        boolean[] excludeSNPs = new boolean[nrSNPsThisChr];

                        // Get the original entries back, so we are sure we dont need to do to many look ups.
                        if (excludeList.size() > 0) {
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (excludeList.contains(valueE2.get(snp).getPrimaryVariantId())) {
                                    excludeSNPs[snp] = true;
                                }
                            }
                        }

                        // Actual scoring.

                        // Loop through the individual risk scores (these should be sorted)
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
                                //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getSequenceName(), riskE.getStartPos());
                                //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                if (var1.getCallRate() < 0.75) {
                                    excludeSNPs[snp] = true;
                                    excludeList.add(riskE.getPrimaryVariantId());
                                    continue;
                                }

                                if (LOGGER.isDebugEnabled()) {
                                    LOGGER.debug(riskE.toString());
                                }

                                // Get the reference allele from this variant 1
                                Allele var1RefAllele = var1.getRefAllele();
                                Allele alternativeAllele;
                                if (var1RefAllele == null) {
                                    var1RefAllele = var1.getAlternativeAlleles().get(0);
                                    alternativeAllele = var1.getAlternativeAlleles().get(1);
                                } else {
                                    alternativeAllele = var1.getAlternativeAlleles().get(0);
                                }

                                if (!String.valueOf(riskE.getAllele()).equals(String.valueOf(alternativeAllele))) {
                                    System.out.printf("Allele match '%s' | '%s' (effect allele | genotype data)%n",
                                            riskE.getAllele(), alternativeAllele);
                                }

                                double or = riskE.getEffectSize();
//                                System.out.println("or = " + or);
                                boolean riskCodedAsTwo;
                                if (sumRisk && or < 0) {
                                    or = or * -1;
                                    riskCodedAsTwo = false;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = true;
                                    }
                                } else {
                                    riskCodedAsTwo = true;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = false;
                                    }
                                }
//                                System.out.println("riskCodedAsTwo = " + riskCodedAsTwo);
//                                    StringBuilder Genos = new StringBuilder();
//                                    StringBuilder Genos1 = new StringBuilder();
//                                    StringBuilder Genos2 = new StringBuilder();
//                                    StringBuilder Genos3 = new StringBuilder();
//                                System.out.println("var1.getSampleCalledDosages() = " + Arrays.toString(var1.getSampleCalledDosages()));
//                                System.out.println("or = " + or);
                                for (int sample = 0; sample < var1.getSampleCalledDosages().length; sample++) {
                                    if (var1.getSampleCalledDosages()[sample] != -1) {
                                        if (riskCodedAsTwo) {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * var1.getSampleCalledDosages()[sample])));
//                                                Genos2.append(or * var1.getSampleCalledDosages()[sample]).append(",");
                                        } else {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * Math.abs(var1.getSampleCalledDosages()[sample] - 2))));
//                                                Genos2.append(or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                        }
//                                            Genos.append(var1.getSampleCalledDosages()[sample]).append(",");
//                                            Genos1.append(var1.getSampleVariants().get(sample).toString());Genos1.append(",");
//                                            Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                    }
                                }
                                //I have now removed the conversion from 0 to 2 and vice versa from the calculations. Why is this needed?
//                                            System.out.println("");
//                                            System.out.println("SNP: "+riskE.getRsName());
//                                            System.out.println("Allele in data: "+var1.getRefAllele().toString());
//                                            System.out.println("Allele: "+riskE.getAllele());
//                                            System.out.println("Diction: "+direction);
//                                            System.out.println("Or: "+or);
//                                            System.out.println("Genotypes: "+Genos1.toString());
//                                            System.out.println("Dosages: "+Genos.toString());
//                                            System.out.println("Scores: "+Genos2.toString());
//                                            System.out.println("ScoresT: "+Genos3.toString());
//                                            System.out.println("");

                                nrSNPs++;

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getStartPos() - riskE.getStartPos()) <= windowSize) {
                                            GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getSequenceName(), riskE2.getStartPos());
                                            if (var2.getCallRate() < 0.75) {
                                                excludeSNPs[t] = true;
                                                excludeList.add(riskE2.getPrimaryVariantId());
                                                continue;
                                            }
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquared) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getPrimaryVariantId());
                                                }
                                            } catch (LdCalculatorException ex) {
                                                LOGGER.error(ex);
                                            }

                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (LOGGER.isDebugEnabled()) {
                        LOGGER.debug("Total SNPs used: " + nrSNPs);
                    }
                }
                p.iterate();
            }
        }
        p.close();

        return scores;
    }

    public DoubleMatrixDataset<String, String> calculateTwoStages(
            DoubleMatrixDataset<String, String> scores,
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks,
            int[] windowSize) {

        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        for (int counter = 0; counter < chrOrder.length; counter++) {
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                HashSet<String> excludeList = new HashSet<String>();

                for (double pVal : pValueThresholds) {
                    String key = "_P" + pVal;
                    THashMap<String, ArrayList<RiskEntry>> riskScorePheno2 = riskScorePheno.getValue().get(key);
                    String NameOfEntry = riskScorePheno.getKey() + key;
                    int rowNr = scores.getHashRows().get(NameOfEntry);

                    if (LOGGER.isDebugEnabled()) {
                        System.out.println(NameOfEntry);

                        LOGGER.debug("SNPs used for GRS calculation:\n");
                    }

                    int nrSNPs = 0;

//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                    if (riskScorePheno2.containsKey(chrOrder[counter])) {
                        ArrayList<RiskEntry> valueE2 = riskScorePheno2.get(chrOrder[counter]);

                        int nrSNPsThisChr = valueE2.size();
                        boolean[] excludeSNPs = new boolean[nrSNPsThisChr];

                        //Get the original entries back, so we are sure we dont need to do to many look ups.
                        if (excludeList.size() > 0) {
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (excludeList.contains(valueE2.get(snp).getPrimaryVariantId())) {
                                    excludeSNPs[snp] = true;
                                }
                            }
                        }
                        //Loop 1, pre-filtering.
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
                                //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getSequenceName(), riskE.getStartPos());
                                //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                if (var1 == null || var1.getCallRate() < 0.75) {
                                    excludeSNPs[snp] = true;
                                    excludeList.add(riskE.getPrimaryVariantId());
                                    continue;
                                }

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getStartPos() - riskE.getStartPos()) <= windowSize[0]) {
                                            GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getSequenceName(), riskE2.getStartPos());
                                            if (var2 == null || var2.getCallRate() < 0.75) {
                                                excludeSNPs[t] = true;
                                                excludeList.add(riskE2.getPrimaryVariantId());
                                                continue;
                                            }
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquared) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getPrimaryVariantId());
                                                }
                                            } catch (LdCalculatorException ex) {
                                                LOGGER.error(ex);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        //Loop 2, Actual scoring.
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
//                                    System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant var1 = genotypeData.getSnpVariantByPos(
                                        riskE.getSequenceName(), riskE.getStartPos());

                                if (LOGGER.isDebugEnabled()) {
                                    LOGGER.debug(riskE.toString());
                                }

                                // Get the reference allele from this variant 1
                                Allele var1RefAllele = var1.getRefAllele();
                                Allele alternativeAllele;
                                if (var1RefAllele == null) {
                                    var1RefAllele = var1.getAlternativeAlleles().get(0);
                                    alternativeAllele = var1.getAlternativeAlleles().get(1);
                                } else {
                                    alternativeAllele = var1.getAlternativeAlleles().get(0);
                                }

//                                System.out.println("p-value = " + riskE.getpValue() + " | beta = " + riskE.getOr());

                                double or = riskE.getEffectSize();
                                boolean riskCodedAsTwo;
                                if (sumRisk && or < 0) {
                                    or = or * -1; // NOTE: please make sure we're using betas here, and not ORS
                                    riskCodedAsTwo = false;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = true;
                                    }
                                } else {
                                    riskCodedAsTwo = true;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = false;
                                    }
                                }

//                                    StringBuilder Genos = new StringBuilder();
//                                    StringBuilder Genos1 = new StringBuilder();
//                                    StringBuilder Genos2 = new StringBuilder();
//                                    StringBuilder Genos3 = new StringBuilder();

                                // Sample called dosages
                                byte[] sampleCalledDosages = var1.getSampleCalledDosages();
                                for (int sample = 0; sample < sampleCalledDosages.length; sample++) {
                                    if (sampleCalledDosages[sample] != -1) {

                                        if (riskCodedAsTwo) {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * sampleCalledDosages[sample])));
//                                                Genos2.append(or * var1.getSampleCalledDosages()[sample]).append(",");
                                        } else {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * Math.abs(sampleCalledDosages[sample] - 2))));
//                                                Genos2.append(or * Math.abs(var1.getSampleCalledDosages()[sample]-2)).append(",");
                                        }

                                        // The above stuff does the following
                                        // If the risk allele matches the reference allele, the first clause is entered
                                        // If the risk allele does not match the reference allele, the alternative should match (biallelic)
                                        // In this case the second clause is entered

                                        // Imagine the following scenario
                                        // -- Genotype:
                                        // ID       REF ALT Dosage_1    Dosage_2    Dosage_3
                                        // rs123    A   C   2           1           0

                                        // With the following risk entry
                                        // -- Risk entry:
                                        // ID       ALT ES
                                        // rs123    C   0.05

                                        // Resulting scores
                                        // Sample_1: |2 - 2| = 0 -> 0 * 0.05 = 0.00
                                        // Sample_2: |2 - 1| = 1 -> 1 * 0.05 = 0.05
                                        // Sample_3: |2 - 0| = 2 -> 2 * 0.05 = 0.10

                                        // OR
                                        // -- Risk entry:
                                        // ID       ALT ES
                                        // rs123    A   -0.05

                                        // -- Resulting scores
                                        // Sample_1: 2 -> 2 * -0.05 = -0.10
                                        // Sample_2: 1 -> 1 * -0.05 = -0.05
                                        // Sample_3: 0 -> 0 * -0.05 = 0.00
//                                            Genos.append(var1.getSampleCalledDosages()[sample]).append(",");
//                                            Genos1.append(var1.getSampleVariants().get(sample).toString());Genos1.append(",");
//                                            Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                    }
                                }

                                //I have now removed the conversion from 0 to 2 and vice versa from the calculations. Why is this needed?
//                                    System.out.println("");
//                                    System.out.println("SNP: "+riskE.getRsName());
//                                    System.out.println("Allele in data: "+var1.getRefAllele().toString());
//                                    System.out.println("Allele: "+riskE.getAllele());
//                                    System.out.println("Risk coded as two: "+riskCodedAsTwo);
//                                    System.out.println("Or: "+or);
//                                    System.out.println("Dosages: "+Genos.toString());
//                                    System.out.println("Genotypes: "+Genos1.toString());
//                                    System.out.println("Scores: "+Genos2.toString());
//                                    System.out.println("ScoresT: "+Genos3.toString());
//                                    System.out.println("");

                                nrSNPs++;

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getStartPos() - riskE.getStartPos()) <= windowSize[1]) {
                                            GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getSequenceName(), riskE2.getStartPos());
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquared) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getPrimaryVariantId());
                                                }
                                            } catch (LdCalculatorException ex) {
                                                LOGGER.error(ex);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (nrSNPs > 0) {
                        LOGGER.info(String.format("Risk entries used for chromosome %s, P-value %f: %d",
                                chrOrder[counter], pVal, nrSNPs));
                    }
                }
                p.iterate();
            }
        }
        p.close();
        return scores;
    }


    public DoubleMatrixDataset<String, String> calculateTwoStages(
            DoubleMatrixDataset<String, String> scores,
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks,
            int[] windowSize, SampleFilterableGenotypeData referenceGenotypeData,
            SampleFilterableGenotypeData responseGenotypeData) {

        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        for (int counter = 0; counter < chrOrder.length; counter++) {
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                HashSet<String> excludeList = new HashSet<String>();

                for (double pVal : pValueThresholds) {
                    String key = "_P" + pVal;
                    THashMap<String, ArrayList<RiskEntry>> riskScorePheno2 = riskScorePheno.getValue().get(key);
                    String NameOfEntry = riskScorePheno.getKey() + key;
                    int rowNr = scores.getHashRows().get(NameOfEntry);

                    if (LOGGER.isDebugEnabled()) {
                        System.out.println(NameOfEntry);

                        LOGGER.debug("SNPs used for GRS calculation:\n");
                    }

                    int nrSNPs = 0;

//                        System.out.println("Processing chromosome:\t" + chrOrder[counter]);
                    if (riskScorePheno2.containsKey(chrOrder[counter])) {
                        ArrayList<RiskEntry> valueE2 = riskScorePheno2.get(chrOrder[counter]);

                        int nrSNPsThisChr = valueE2.size();
                        boolean[] excludeSNPs = new boolean[nrSNPsThisChr];

                        //Get the original entries back, so we are sure we dont need to do to many look ups.
                        if (excludeList.size() > 0) {
                            for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                                if (excludeList.contains(valueE2.get(snp).getPrimaryVariantId())) {
                                    excludeSNPs[snp] = true;
                                }
                            }
                        }
                        //Loop 1, pre-filtering.
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
                                //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant var1 = referenceGenotypeData.getSnpVariantByPos(riskE.getSequenceName(), riskE.getStartPos());
                                //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                if (var1 == null || var1.getCallRate() < 0.75) {
                                    excludeSNPs[snp] = true;
                                    excludeList.add(riskE.getPrimaryVariantId());
                                    continue;
                                }

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getStartPos() - riskE.getStartPos()) <= windowSize[0]) {
                                            GeneticVariant var2 = referenceGenotypeData.getSnpVariantByPos(riskE2.getSequenceName(), riskE2.getStartPos());
                                            if (var2 == null || var2.getCallRate() < 0.75) {
                                                excludeSNPs[t] = true;
                                                excludeList.add(riskE2.getPrimaryVariantId());
                                                continue;
                                            }
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquared) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getPrimaryVariantId());
                                                }
                                            } catch (LdCalculatorException ex) {
                                                LOGGER.error(ex);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        //Loop 2, Actual scoring.
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
//                                    System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant originalVariant = genotypeData.getSnpVariantByPos(
                                        riskE.getSequenceName(), riskE.getStartPos());

                                SampleFilteredReadOnlyGeneticVariant referenceVariantOne =
                                        new SampleFilteredReadOnlyGeneticVariant(originalVariant, referenceGenotypeData);

                                if (LOGGER.isDebugEnabled()) {
                                    LOGGER.debug(riskE.toString());
                                }

                                // Get the reference allele from this variant 1
                                Allele var1RefAllele = referenceVariantOne.getRefAllele();
                                if (var1RefAllele == null) {
                                    var1RefAllele = referenceVariantOne.getAlternativeAlleles().get(0);
                                }

//                                System.out.println("p-value = " + riskE.getpValue() + " | beta = " + riskE.getOr());

                                double or = riskE.getEffectSize();
                                boolean riskCodedAsTwo;
                                if (sumRisk && or < 0) {
                                    or = or * -1; // NOTE: please make sure we're using betas here, and not ORS
                                    riskCodedAsTwo = false;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = true;
                                    }
                                } else {
                                    riskCodedAsTwo = true;
                                    if (!(riskE.getAlleleAsSnp() == (var1RefAllele.getAlleleAsSnp())
                                            || riskE.getAlleleAsSnp() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = false;
                                    }
                                }

                                SampleFilteredReadOnlyGeneticVariant sampleFilteredVariantOne =
                                        new SampleFilteredReadOnlyGeneticVariant(originalVariant, responseGenotypeData);

//                                    StringBuilder Genos = new StringBuilder();
//                                    StringBuilder Genos1 = new StringBuilder();
//                                    StringBuilder Genos2 = new StringBuilder();
//                                    StringBuilder Genos3 = new StringBuilder();

                                // Sample called dosages
                                byte[] sampleCalledDosages = sampleFilteredVariantOne.getSampleCalledDosages();
                                for (int sample = 0; sample < sampleCalledDosages.length; sample++) {
                                    if (sampleCalledDosages[sample] != -1) {

                                        if (riskCodedAsTwo) {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * sampleCalledDosages[sample])));
//                                                Genos2.append(or * referenceVariantOne.getSampleCalledDosages()[sample]).append(",");
                                        } else {
                                            scores.getMatrix().setQuick(rowNr, sample, (scores.getMatrix().getQuick(rowNr, sample) + (or * Math.abs(sampleCalledDosages[sample] - 2))));
//                                                Genos2.append(or * Math.abs(referenceVariantOne.getSampleCalledDosages()[sample]-2)).append(",");
                                        }

                                        // The above stuff does the following
                                        // If the risk allele matches the reference allele, the first clause is entered
                                        // If the risk allele does not match the reference allele, the alternative should match (biallelic)
                                        // In this case the second clause is entered

                                        // Imagine the following scenario
                                        // -- Genotype:
                                        // ID       REF ALT Dosage_1    Dosage_2    Dosage_3
                                        // rs123    A   C   2           1           0

                                        // With the following risk entry
                                        // -- Risk entry:
                                        // ID       ALT ES
                                        // rs123    C   0.05

                                        // Resulting scores
                                        // Sample_1: |2 - 2| = 0 -> 0 * 0.05 = 0.00
                                        // Sample_2: |2 - 1| = 1 -> 1 * 0.05 = 0.05
                                        // Sample_3: |2 - 0| = 2 -> 2 * 0.05 = 0.10

                                        // OR
                                        // -- Risk entry:
                                        // ID       ALT ES
                                        // rs123    A   -0.05

                                        // -- Resulting scores
                                        // Sample_1: 2 -> 2 * -0.05 = -0.10
                                        // Sample_2: 1 -> 1 * -0.05 = -0.05
                                        // Sample_3: 0 -> 0 * -0.05 = 0.00
//                                            Genos.append(referenceVariantOne.getSampleCalledDosages()[sample]).append(",");
//                                            Genos1.append(referenceVariantOne.getSampleVariants().get(sample).toString());Genos1.append(",");
//                                            Genos3.append(scores.getMatrix().getQuick(rowNr, sample)).append(",");
                                    }
                                }

                                //I have now removed the conversion from 0 to 2 and vice versa from the calculations. Why is this needed?
//                                    System.out.println("");
//                                    System.out.println("SNP: "+riskE.getRsName());
//                                    System.out.println("Allele in data: "+referenceVariantOne.getRefAllele().toString());
//                                    System.out.println("Allele: "+riskE.getAllele());
//                                    System.out.println("Risk coded as two: "+riskCodedAsTwo);
//                                    System.out.println("Or: "+or);
//                                    System.out.println("Dosages: "+Genos.toString());
//                                    System.out.println("Genotypes: "+Genos1.toString());
//                                    System.out.println("Scores: "+Genos2.toString());
//                                    System.out.println("ScoresT: "+Genos3.toString());
//                                    System.out.println("");

                                nrSNPs++;

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getStartPos() - riskE.getStartPos()) <= windowSize[1]) {
                                            GeneticVariant var2 = referenceGenotypeData.getSnpVariantByPos(riskE2.getSequenceName(), riskE2.getStartPos());
                                            try {
                                                if (calculateRsquare(referenceVariantOne, var2) >= rSquared) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getPrimaryVariantId());
                                                }
                                            } catch (LdCalculatorException ex) {
                                                LOGGER.error(ex);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (nrSNPs > 0) {
                        LOGGER.info(String.format("Risk entries used for chromosome %s, P-value %f: %d",
                                chrOrder[counter], pVal, nrSNPs));
                    }
                }
                p.iterate();
            }
        }
        p.close();
        return scores;
    }

    public DoubleMatrixDataset<String, String> initializePolygenicScoreMatrix(
            RandomAccessGenotypeData genotypeData,
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks) {

        // Initialize an array list with concatenated strings with a combination of P-values and gwas id.
        ArrayList<String> keys = new ArrayList<String>();
        for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, THashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno.getKey() + riskScorePheno2.getKey());
            }
        }

        // Initialize
        return new DoubleMatrixDataset<>(keys, Arrays.asList(genotypeData.getSampleNames()));
    }

    public DoubleMatrixDataset<String, String> initializePolygenicScoreMatrix(String phenotype) {
        ArrayList<String> keys = new ArrayList<>();
        for (Double pValue : pValueThresholds) {
            keys.add(phenotype + pValue);
        }

        return new DoubleMatrixDataset<>(keys, Arrays.asList(genotypeData.getSampleNames()));
    }

    public String[] getGenomicRangesToExclude() {
        return genomicRangesToExclude;
    }

    public List<Double> getpValueThresholds() {
        return pValueThresholds;
    }

    public THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskEntries(
            GwasSummaryStatistics summaryStatistics,
            RandomAccessGenotypeData genotypeData,
            double[] pValueThreshold, String[] genomicRangesToExclude, boolean unweighted) {
        THashMap<String, ArrayList<Pair<Integer, Integer>>> exclussionRanges = new THashMap<>();

        if (genomicRangesToExclude != null) {
            System.out.println("Trying to exclude genomic ranges.");
            int ranges = 0;
            for (String s : genomicRangesToExclude) {
                String[] parts = s.split(":");
                String key = parts[0];
                parts = parts[1].split("-");

                if (!exclussionRanges.contains(key)) {
                    exclussionRanges.put(key, new ArrayList<Pair<Integer, Integer>>());
                }
                exclussionRanges.get(key).add(new Pair(Integer.parseInt(parts[0]), Integer.parseInt(parts[1])));
                ranges++;
            }
            if (LOGGER.isDebugEnabled()) {
                System.out.println("Number of ranges excluded: " + ranges + " on: " + exclussionRanges.size() + " chromosomes");
            }
        }
        // A risk entry (value?) per variant? per sequence? per pval threshold? per file?
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = new THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>>();

        THashSet<String> chromosomesExcluded = new THashSet<>();
        int snpsExcluded = 0;
        try {
            
            String name = summaryStatistics.getGwasId();

            THashMap<String, THashMap<String, ArrayList<RiskEntry>>> filehash = new THashMap<>();

            for (double p : pValueThreshold) {
                String name2 = "_P" + p;
                if (!filehash.containsKey(name2)) {
                    filehash.put(name2, new THashMap<>());
                }
            }

            int numberOfVariantsWithoutSameAlleles = 0;
            int numberOfVariantsWithoutSameAllelesComplement = 0;

            for (EffectAllele effectAllele : summaryStatistics) {
                // Check if there is a effectAllele in the genotype data corresponding to this effectAllele
                GeneticVariant snpVariantByPos = genotypeData.getSnpVariantByPos(
                        effectAllele.getSequenceName(),
                        effectAllele.getStartPos());

                if (!effectAllele.matchesVariant(snpVariantByPos)) {
                    numberOfVariantsWithoutSameAlleles++;
                    continue;
                }

//                    System.out.println(s);
//                        System.out.print(snpObject.getSequenceName() + "\t" + snpObject.getStartPos() + "\n");
                double currentP = effectAllele.getPValue();
                boolean addEntry = true;
                float partsTwo = (float) effectAllele.getEffectSize();

                if (unweighted) {
                    if (partsTwo < 0) {
                        partsTwo = -1;
                    } else {
                        partsTwo = 1;
                    }
                }

                if (exclussionRanges.contains(effectAllele.getSequenceName())) {
                    chromosomesExcluded.add(effectAllele.getSequenceName());
                    for (Pair<Integer, Integer> p : exclussionRanges.get(effectAllele.getSequenceName())) {
                        if (p.getLeft() <= effectAllele.getStartPos() && p.getRight() >= effectAllele.getStartPos()) {
                            addEntry = false;
                            snpsExcluded++;
                        }
                    }
                }

                if (addEntry) {
                    for (double p : pValueThreshold) {
                        if (currentP < p) {
                            String name2 = "_P" + p;

                            if (!filehash.get(name2).containsKey(effectAllele.getSequenceName())) {
                                filehash.get(name2).put(effectAllele.getSequenceName(), new ArrayList<>());
                            }
                            filehash.get(name2).get(effectAllele.getSequenceName()).add(new RiskEntry(effectAllele.getPrimaryVariantId(),
                                    effectAllele.getSequenceName(), effectAllele.getStartPos(),
                                    effectAllele.getAllele().getAlleleAsSnp(), partsTwo, currentP));
                        }
                    }
                }
            }
            synchronized (risks) {
                risks.put(name, filehash);
            }

            if (numberOfVariantsWithoutSameAlleles > 0) {
                String message = String.format(
                        "%d variants have different alleles than the alleles of " +
                                "the matched variant from the input genotype data.%n" +
                                "The complement of these alleles is also different for %d variants, " +
                                "which will be removed.",
                        numberOfVariantsWithoutSameAlleles, numberOfVariantsWithoutSameAllelesComplement);
                System.err.println(message);
                LOGGER.warn(message);
            }

            if (LOGGER.isDebugEnabled()) {
                System.out.println("Chromosomes where regions are excluded: " + chromosomesExcluded);
                System.out.println("Number of SNPs excluded: " + snpsExcluded);
            }
        } catch (GwasSummaryStatisticsException ex) {
            LOGGER.error(ex);
        }

        for (Map.Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> e : risks.entrySet()) {

            for (Map.Entry<String, THashMap<String, ArrayList<RiskEntry>>> e2 : e.getValue().entrySet()) {
                int entries = 0;
                for (Map.Entry<String, ArrayList<RiskEntry>> e3 : e2.getValue().entrySet()) {
                    Collections.sort(e3.getValue());
                    entries += e3.getValue().size();
                }
//                System.out.println(e.getKey()+e2.getKey()+" has: "+entries+" entries");
            }

        }
        return risks;
    }
}
