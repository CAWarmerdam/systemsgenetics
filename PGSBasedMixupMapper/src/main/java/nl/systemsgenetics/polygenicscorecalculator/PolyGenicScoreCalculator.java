package nl.systemsgenetics.polygenicscorecalculator;

import nl.systemsgenetics.gwassummarystatistics.EffectAllele;
import nl.systemsgenetics.gwassummarystatistics.GwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.MultiStudyGwasSummaryStatistics;
import nl.systemsgenetics.gwassummarystatistics.VcfGwasSummaryStatistics;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class PolyGenicScoreCalculator {

    private static final Logger LOGGER = Logger.getLogger(PolyGenicScoreCalculator.class);
    private LDHandler ldHandler;
    private String[] genomicRangesToExclude;
    private List<Double> pValueThresholds;

    public PolyGenicScoreCalculator(LDHandler ldHandler,
                                    String[] genomicRangesToExclude) {
        this.ldHandler = ldHandler;
        this.genomicRangesToExclude = genomicRangesToExclude;
    }

    public DoubleMatrixDataset<String, Double> calculate(RandomAccessGenotypeData genotypeData,
                                                         GwasSummaryStatistics summaryStatistics) throws PolyGenicScoreCalculatorException {
        // Iterate over methods from the shrinkage strategy

        // Dependent on the LD handler, either get a list of unlinked variants, or
        // use the LD matrix and use it to correct for linkage

        Iterator<EffectAllele> effectAlleles = summaryStatistics.effectAlleles();

        return calculatePolyGenicScores(genotypeData, effectAlleles);
    }

    private DoubleMatrixDataset<String, Double> calculatePolyGenicScores(
            RandomAccessGenotypeData genotypeData,
            Iterator<EffectAllele> effectAlleles)
            throws PolyGenicScoreCalculatorException {

        // Initialize a scores array
        DoubleMatrixDataset<String, Double> scores =
                new DoubleMatrixDataset<>(Arrays.asList(genotypeData.getSampleNames()), pValueThresholds);

        // Make sure that the effect alleles are sorted by p-values, per sequence
        Iterator<EffectAllele> sortedEffectAlleleIterator = EffectAllele.sortEffectAllelesPerSequence(effectAlleles);

        // Get the iterator from the LD handler, making sure that only unlinked effect alleles are used.
        Iterator<EffectAllele> unlinkedEffectAlleleIterator = ldHandler.effectAlleleIterator(sortedEffectAlleleIterator);

        // Loop through the risk entries that are not linked
        for (; unlinkedEffectAlleleIterator.hasNext(); ) {
            EffectAllele effectAllele = unlinkedEffectAlleleIterator.next();


            // Get the genetic variant corresponding to the genotype data.
            GeneticVariant sampleGeneticVariant = extractSampleGeneticVariant(genotypeData,
                    effectAllele.getVariant());
            if (sampleGeneticVariant == null) continue;

            // Check if the risk entry / genetic variant is valid for use
            if (isValidEffectAllele(sampleGeneticVariant)) {
                // The effect allele is valid, calculate the PGSs for every p-value threshold.
                performThresholdedPgsCalculation(scores, effectAllele, sampleGeneticVariant);
            }
        }

        // if the effect sizes used as weights are reported as log odds ratios (log(OR)),
        // The PRS will take the same units as the logarithmic scale so this has to be transformed back to an OR scale.

        // The PRS can be standardised by dividing by the number of SNPs,
        // or standardised to a standard normal distribution
        return scores;
    }

    private void performThresholdedPgsCalculation(DoubleMatrixDataset<String, Double> scores,
                                                  EffectAllele effectAllele, GeneticVariant sampleGeneticVariant) {
        for (int thresholdIndex = 0; thresholdIndex < pValueThresholds.size(); thresholdIndex++) {
            double threshold = -Math.log10(pValueThresholds.get(thresholdIndex));
            if (effectAllele.getLogTransformedPValue() > threshold) {

                System.out.println(effectAllele.toString());
                // Get the dosages per sample for this variant
                byte[] sampleCalledDosages = sampleGeneticVariant.getSampleCalledDosages();

                // Get the effect size / odds ratio for this risk entry
                double effectSize = effectAllele.getEffectSize();

                boolean refAlleleMatchesEffectAllele = isRefAlleleMatchingEffectAllele(sampleGeneticVariant, effectAllele);

                // Loop through the samples to calculate the PGSs
                for (int sampleIndex = 0; sampleIndex < sampleCalledDosages.length; sampleIndex++) {
                    byte dosage = sampleCalledDosages[sampleIndex];
                    if (dosage != -1) {

                        // Calculate the PGS as the product of the dosage and the effect size summed for every risk entry
                        double variantEffect;

                        if (refAlleleMatchesEffectAllele) {
                            variantEffect = effectSize * dosage;
                        } else {
                            // Assume the other allele matches the effect allele
                            variantEffect = effectSize * Math.abs(dosage - 2);
                        }
                        scores.setElementQuick(sampleIndex, thresholdIndex,
                                scores.getElementQuick(sampleIndex, thresholdIndex) + variantEffect);
                    }
                }
            }
        }
    }

    private boolean isRefAlleleMatchingEffectAllele(GeneticVariant sampleGeneticVariant, EffectAllele effectAllele) {
        Allele refAllele = sampleGeneticVariant.getVariantAlleles().get(0);
        return refAllele.equals(effectAllele.getAllele());
    }

    private boolean isValidEffectAllele(GeneticVariant sampleGeneticVariant) {
        return !(sampleGeneticVariant.getCallRate() < 0.75);
    }

    private GeneticVariant extractSampleGeneticVariant(RandomAccessGenotypeData genotypeData,
                                                       GeneticVariant geneticVariant) throws PolyGenicScoreCalculatorException {
        for (GeneticVariant candidateVariant : genotypeData.getVariantsByPos(
                geneticVariant.getSequenceName(), geneticVariant.getStartPos())) {
            if (candidateVariant.getVariantAlleles().sameAlleles(geneticVariant.getVariantAlleles())) {
                return candidateVariant;
            }
        }
        return null;
//        throw new PolyGenicScoreCalculatorException(
//                "Not able to find variant in genotype data corresponding to risk allele");
    }

    public void setPValueThresholds(List<Double> pValueThresholds) {
        this.pValueThresholds = pValueThresholds;
    }

    public List<Double> getpValueThresholds() {
        return pValueThresholds;
    }

    public LDHandler getLdHandler() {
        return ldHandler;
    }

    public String[] getGenomicRangesToExclude() {
        return genomicRangesToExclude;
    }

    public void setGenomicRangesToExclude(String[] genomicRangesToExclude) {
        this.genomicRangesToExclude = genomicRangesToExclude;
    }
}
