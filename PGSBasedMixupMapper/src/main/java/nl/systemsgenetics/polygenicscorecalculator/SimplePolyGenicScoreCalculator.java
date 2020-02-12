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
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
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
public class SimplePolyGenicScoreCalculator {

    private static final Logger LOGGER = Logger.getLogger(SimplePolyGenicScoreCalculator.class);
    private static final String[] chrOrder = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

    public static DoubleMatrixDataset<String, String> calculate(
            RandomAccessGenotypeData genotypeData,
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks,
            double rSquare, double windowSize, double[] pValueThreshold, boolean sumRisk) {

        // Select the keys that correspond to the different combinations to test for
        ArrayList<String> keys = new ArrayList<>();
        for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, THashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno2.getKey());
            }
        }

        // Initialize the scores matrix
        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(keys, Arrays.asList(genotypeData.getSampleNames()));

        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        // Loop through the chromosomes
        for (int counter = 0; counter < chrOrder.length; counter++) {
            // Loop through the risk entries for every file
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                // Loop through the p-values
                HashSet<String> excludeList = new HashSet<String>();
                for (double pVal : pValueThreshold) {
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
                                if (excludeList.contains(valueE2.get(snp).getRsName())) {
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
                                GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                                //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                if (var1.getCallRate() < 0.75) {
                                    excludeSNPs[snp] = true;
                                    excludeList.add(riskE.getRsName());
                                    continue;
                                }

                                if (LOGGER.isDebugEnabled()) {
                                    LOGGER.debug(riskE.InfoToString());
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

                                double or = riskE.getOr();
//                                System.out.println("or = " + or);
                                boolean riskCodedAsTwo;
                                if (sumRisk && or < 0) {
                                    or = or * -1;
                                    riskCodedAsTwo = false;
                                    if (!(riskE.getAllele() == (var1RefAllele.getAlleleAsSnp()) || riskE.getAllele() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = true;
                                    }
                                } else {
                                    riskCodedAsTwo = true;
                                    if (!(riskE.getAllele() == (var1RefAllele.getAlleleAsSnp()) || riskE.getAllele() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
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
                                        if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize) {
                                            GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                            if (var2.getCallRate() < 0.75) {
                                                excludeSNPs[t] = true;
                                                excludeList.add(riskE2.getRsName());
                                                continue;
                                            }
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquare) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getRsName());
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

    public static DoubleMatrixDataset<String, String> calculateTwoStages(RandomAccessGenotypeData genotypeData, THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks, double rSquare, int[] windowSize, double[] pValueThreshold, boolean sumRisk) {
        ArrayList<String> keys = new ArrayList<String>();
        for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
            for (Entry<String, THashMap<String, ArrayList<RiskEntry>>> riskScorePheno2 : riskScorePheno.getValue().entrySet()) {
                keys.add(riskScorePheno.getKey() + riskScorePheno2.getKey());
            }
        }

        DoubleMatrixDataset<String, String> scores = new DoubleMatrixDataset<String, String>(keys, Arrays.asList(genotypeData.getSampleNames()));
        ProgressBar p = new ProgressBar(risks.size() * chrOrder.length);

        for (int counter = 0; counter < chrOrder.length; counter++) {
            for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskScorePheno : risks.entrySet()) {
                HashSet<String> excludeList = new HashSet<String>();

                for (double pVal : pValueThreshold) {
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
                                if (excludeList.contains(valueE2.get(snp).getRsName())) {
                                    excludeSNPs[snp] = true;
                                }
                            }
                        }
                        //Loop 1, pre-filtering.
                        for (int snp = 0; snp < nrSNPsThisChr; snp++) {
                            if (!excludeSNPs[snp]) {
                                RiskEntry riskE = valueE2.get(snp);
                                //System.out.println(snpID + "\t" + c + "\t" + chrPos + "\t" + object.doubleValue);
                                GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                                //Check if at least 75% of the sampels have information for the SNP otherwise it is removed by default.
                                if (var1 == null || var1.getCallRate() < 0.75) {
                                    excludeSNPs[snp] = true;
                                    excludeList.add(riskE.getRsName());
                                    continue;
                                }

                                for (int t = snp + 1; t < nrSNPsThisChr; t++) {
                                    if (!excludeSNPs[t]) {
                                        RiskEntry riskE2 = valueE2.get(t);
                                        if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize[0]) {
                                            GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                            if (var2 == null || var2.getCallRate() < 0.75) {
                                                excludeSNPs[t] = true;
                                                excludeList.add(riskE2.getRsName());
                                                continue;
                                            }
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquare) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getRsName());
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
                                GeneticVariant var1 = genotypeData.getSnpVariantByPos(riskE.getChr(), riskE.getPos());
                                if (LOGGER.isDebugEnabled()) {
                                    LOGGER.debug(riskE.InfoToString());
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

                                double or = riskE.getOr();
                                boolean riskCodedAsTwo;
                                if (sumRisk && or < 0) {
                                    or = or * -1; // NOTE: please make sure we're using betas here, and not ORS
                                    riskCodedAsTwo = false;
                                    if (!(riskE.getAllele() == (var1RefAllele.getAlleleAsSnp()) || riskE.getAllele() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
                                        riskCodedAsTwo = true;
                                    }
                                } else {
                                    riskCodedAsTwo = true;
                                    if (!(riskE.getAllele() == (var1RefAllele.getAlleleAsSnp()) || riskE.getAllele() == (var1RefAllele.getComplement().getAlleleAsSnp()))) {
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
                                        if (Math.abs(riskE2.getPos() - riskE.getPos()) <= windowSize[1]) {
                                            GeneticVariant var2 = genotypeData.getSnpVariantByPos(riskE2.getChr(), riskE2.getPos());
                                            try {
                                                if (calculateRsquare(var1, var2) >= rSquare) {
                                                    excludeSNPs[t] = true;
                                                    excludeList.add(riskE2.getRsName());
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
}
