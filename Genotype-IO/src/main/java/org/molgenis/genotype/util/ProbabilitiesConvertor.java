/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import java.util.*;
import java.util.stream.IntStream;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.bgen.BgenGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public class ProbabilitiesConvertor {

	/**
	 * Convert alleles to probability values. If more than two alleles will
	 * return all missing. If sample alleles found not specified by alleles then
	 * sample will be set to missing
	 *
	 *
	 * @param sampleAlleles
	 * @param alleles
	 * @return
	 */
	public static float[][] convertCalledAllelesToProbability(List<Alleles> sampleAlleles, Alleles alleles) {

		float[][] probs = new float[sampleAlleles.size()][3];

		if (alleles.getAlleleCount() == 0) {
			throw new GenotypeDataException("Error converting alleles to probabilities. No alleles detected");
		} else if (alleles.getAlleleCount() > 2) {

			Arrays.fill(probs, new float[]{0,0,0});
			return probs;

		} else {

			Alleles AA = Alleles.createAlleles(alleles.get(0), alleles.get(0));
			Alleles AB;
			Alleles BA;
			Alleles BB;

			if (alleles.getAlleleCount() == 2) {
				AB = Alleles.createAlleles(alleles.get(0), alleles.get(1));
				BA = Alleles.createAlleles(alleles.get(1), alleles.get(0));
				BB = Alleles.createAlleles(alleles.get(1), alleles.get(1));
			} else {
				AB = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
				BA = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
				BB = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
			}

			for (int i = 0; i < probs.length; ++i) {
				Alleles sampleVariant = sampleAlleles.get(i);

				if (sampleVariant.contains(Allele.ZERO)) {
					probs[i] = new float[]{0,0,0};
				} else if (sampleVariant == AA) {
					probs[i] = new float[]{1,0,0};
				} else if (sampleVariant == AB) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = new float[]{0,1,0};
				} else if (sampleVariant == BA) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = new float[]{0,1,0};
				} else if (sampleVariant == BB) { //We test above for missing so is never true for monomorphic SNPs
					probs[i] = new float[]{0,0,1};
				} else {
					probs[i] = new float[]{0,0,0};
				}

			}

			return probs;

		}

	}

	/**
	 * Uses inexact heuristic as defined by plinkseq
	 * (thtp://atgu.mgh.harvard.edu/plinkseq/dosage.shtml) to convert dosage to
	 * probabilities.
	 *
	 * @param sampleDosages assuming count of ref allele.
	 * @return
	 */
	public static float[][] convertDosageToProbabilityHeuristic(float[] sampleDosages) {

		float[][] probs = new float[sampleDosages.length][3];

		for (int i = 0; i < sampleDosages.length; ++i) {

			float sampleDosage = sampleDosages[i];

			if (sampleDosage > 2 || sampleDosage < 0) {
				probs[i] = new float[]{0f, 0f, 0f};
			} else if (sampleDosage < 1) {
				probs[i] = new float[]{0, sampleDosage, 1 - sampleDosage};
			} else {
				//sampleDosage >= 1 && sampleDosage <= 2
				probs[i] = new float[]{sampleDosage - 1, 2 - sampleDosage, 0};
			}

		}

		return probs;

	}

	/**
	 *
	 * @param probs
	 * @param variantAlleles the two alleles for this variant
	 * @param minProbability to call a genotype
	 * @return
	 */
	public static List<Alleles> convertProbabilitiesToAlleles(float[][] probs, Alleles variantAlleles, double minProbability) {

		ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>(probs.length);

		final int alleleCount = variantAlleles.getAlleleCount();

		if (alleleCount > 2 || alleleCount == 0) {
			throw new GenotypeDataException("Error converting posterior probabilities to called alleles. Found non biallelic SNP");
		}

		Alleles aa = Alleles.createAlleles(variantAlleles.get(0), variantAlleles.get(0));
		Alleles bb;

		if (alleleCount == 2) {
			bb = Alleles.createAlleles(variantAlleles.get(1), variantAlleles.get(1));
		} else {
			bb = null;
		}

		Alleles missing = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);

		for (float[] sampleProbs : probs) {

			int maxProbIndex = -1;
			float maxProb = 0;

			int i = 0;
			for (float prob : sampleProbs) {
				if (prob > 0 && prob >= minProbability && prob > maxProb) {
					maxProbIndex = i;
					maxProb = prob;
				}
				++i;
			}

			if (alleleCount == 1 && maxProbIndex >= 1) {
				throw new GenotypeDataException("Error converting posterior probabilities to called alleles. Illigale probability.");
			}

			switch (maxProbIndex) {
				case -1:
					sampleAlleles.add(missing);
					break;
				case 0:
					sampleAlleles.add(aa);
					break;
				case 1:
					sampleAlleles.add(variantAlleles);
					break;
				case 2:
					sampleAlleles.add(bb);
					break;
				default:
					throw new GenotypeDataException("Error converting posterior probabilities to called alleles. This should not happen, please report this bug.");
			}

		}

		return sampleAlleles;

	}

	public static float[] convertProbabilitiesToDosage(float[][] probs, double minProbability) {

		float[] dosages = new float[probs.length];

		for (int i = 0; i < probs.length; ++i) {

			boolean containsMinProbability = false;

			for (float prob : probs[i]) {
				if (prob >= minProbability) {
					containsMinProbability = true;
					break;
				}
			}

			if (containsMinProbability) {

				dosages[i] = (probs[i][0] * 2) + probs[i][1];
				if (dosages[i] > 2) {
					dosages[i] = 2;
				}

			} else {
				dosages[i] = -1;
			}
		}

		return dosages;


	}

	/**
	 * Method for converting complex probabilities to regular posterior probabilities.
	 * unlike the regular posterior probabilities, probabilities in BGEN can be stored
	 * in an array of arbitrary length corresponding to the number of ordered combinations of
	 * alleles for a given ploidy. Any number of alleles can be represented.
	 *
	 * If the complex probabilities for a sample represent 3 possible genotypes
	 * (diploid samples for a biallelic variant) This is returned.
	 *
	 * If the complex probabilities represent less than 3 or more than 3 possible genotypes
	 * an empty array of size 3 is returned for the sample.
	 *
	 * If the complex probabilities represent 2 possible genotypes (biallelic for haploid samples)
	 * this is <b>not</b> expanded as if diploid with a zero probability for heterozygosity.
	 * Instead missingness is returned.
	 *
	 * @param complexProbabilities The probabilities returned by a getSampleProbabilitiesComplex method, only for
	 *                             biallelic variants.
	 * @return An array of arrays of size 3 with posterior probabilities.
	 */
	public static float[][] convertBiallelicComplexProbabilitiesToProbabilities(double[][] complexProbabilities) {
		// Define an array consisting of an array of posterior complexProbabilities for each genotype
		float[][] probabilities = new float[complexProbabilities.length][3];

		for (int sampleIndex = 0; sampleIndex < complexProbabilities.length; sampleIndex++) {
			// Get the array of doubles
			double[] sampleProbabilitiesBgen = complexProbabilities[sampleIndex];
			// Get the length of the probabilities array
			int probabilitiesArrayLength = sampleProbabilitiesBgen.length;
			// Initialize empty array for float values.
			float[] sampleProbabilities = new float[3];
			// Check for the length of the sample probabilities.
			if (probabilitiesArrayLength == 3) {
				// Convert the array of doubles to the array of floats
				IntStream.range(0, probabilitiesArrayLength)
						.forEach(index -> sampleProbabilities[index] = (float) sampleProbabilitiesBgen[index]);
			}
			// Currently returning missing probabilities when the number of probabilities is equal to 2, as .
			// Uncomment to recode the 2 probabilities as if it represents a diploid sample.
//			else if (probabilitiesArrayLength == 2) {
//				sampleProbabilities[0] = (float) sampleProbabilitiesBgen[0];
//				sampleProbabilities[2] = (float) sampleProbabilitiesBgen[1];
//			}
			// If probabilities array length is less than 2 or greater than 3,
			// just return an array of zeros. [0.0, 0.0, 0.0]

			// Insert the probabilities for this sample
			probabilities[sampleIndex] = sampleProbabilities;
		}

		return probabilities;
	}

	/**
	 * Method that converts between posterior probabilities in its legacy format
	 * (3 probabilities, float values)
	 * and complex posterior probabilities that can contain an arbitrary number of probabilities.
	 * This method effectively only converts between float and double values.
	 *
	 * @param probabilities The probabilities to convert.
	 * @return An array of probabilities as doubles for every individual in an array.
	 */
	public static double[][] convertProbabilitiesToComplexProbabilities(float[][] probabilities) {
		double[][] bgenProbabilities = new double[probabilities.length][probabilities[0].length];
		for (int i = 0; i < probabilities.length; i++) {
			for (int j = 0; j < probabilities[0].length; j++) {
				bgenProbabilities[i][j] = probabilities[i][j];
			}
		}
		return bgenProbabilities;
	}

	/**
	 * Converts phased probabilities to complex / bgen probabilities.
	 *
	 * @param haplotypeProbabilities The haplotype probabilities
	 * @return The probabilities per genotype per sample.
	 */
	public static double[][] convertPhasedProbabilitiesToComplexProbabilities(double[][][] haplotypeProbabilities) {
		double[][] probabilities;
		probabilities = new double[haplotypeProbabilities.length][];
		// Calculate the probability for homozygous genotype 'AA',
		// the probability for 'AB', and the probability for 'BB'
		for (int sampleIndex = 0; sampleIndex < haplotypeProbabilities.length; sampleIndex++) {
			// Convert the probabilities for this particular sample
			probabilities[sampleIndex] = phasedSampleProbabilitiesToGenotypeProbabilities(
					haplotypeProbabilities[sampleIndex], haplotypeProbabilities[sampleIndex][0].length);
		}
		return probabilities;
	}

	/**
	 * Converts an array of arrays with probabilities per allele for particular haplotypes, to
	 * regular genotype probability values per genotype.
	 *
	 * @param phasedProbabilities An array, with its size equal to the number haplotypes, with arrays
	 *                            with a probability for every possible allele for the haplotype.
	 * @param numberOfAlleles     The number of alleles.
	 * @return The number of probabilities for every possible genotype.
	 */
	private static double[] phasedSampleProbabilitiesToGenotypeProbabilities(double[][] phasedProbabilities,
																			 int numberOfAlleles) {
		// Get all possible combinations of alleles for the haplotypes (all permutations)
		List<List<Integer>> haplotypeCombinations = BgenGenotypeData.getHaplotypeCombinations(
				numberOfAlleles, phasedProbabilities.length);

		// Initialize an array of probabilities.
		double[] genotypeProbabilities = new double[BgenGenotypeData.numberOfProbabilitiesForPloidyAlleleCountCombination(
				phasedProbabilities.length,
				numberOfAlleles - 1)]; // number of alleles is equal to n, n-1 = numberOfAlleles - 1

		// Loop through the combinations of haplotypes
		for (List<Integer> haplotypeCombination : haplotypeCombinations) {
			double probability = 1;
			// Multiply the probabilities that correspond to the indices within the haplotype combination
			for (int haplotypeIndex = 0; haplotypeIndex < haplotypeCombination.size(); haplotypeIndex++) {
				probability *= phasedProbabilities[haplotypeIndex][haplotypeCombination.get(haplotypeIndex)];
			}
			// Add the multiplied probability of this combination to the other combinations with the same genotype
			Collections.sort(haplotypeCombination);
			genotypeProbabilities[BgenGenotypeData.getIndexOfOrderedCombination(
					Collections.unmodifiableList(haplotypeCombination))] += probability;
		}
		return genotypeProbabilities;
	}

	/**
	 * Converts between allele dosages per haplotype per sample, and probabilities per haplotype per sample.
	 * The dosage value is assumed to correspond to the alternative allele.
	 *
	 * @param haplotypeDosages An array of dosages for the alternative allele per haplotype, per sample.
	 * @return An array of probabilities per haplotype, per sample
	 */
	public static double[][][] haplotypeDosagesToHaplotypeProbabilities(double[][] haplotypeDosages) {
		double[][][] probs = new double[haplotypeDosages.length][][];

		// Loop through the samples
		for (int sampleIndex = 0; sampleIndex < haplotypeDosages.length; sampleIndex++) {

			// Get the dosages per haplotype for this sample
			double[] haplotypeDosage = haplotypeDosages[sampleIndex];
			// Initialize a nested array for probabilities corresponding to this sample
			double[][] sampleProbabilities = new double[haplotypeDosage.length][2];
			for (int haplotypeIndex = 0; haplotypeIndex < haplotypeDosage.length; haplotypeIndex++) {
				// Get the allele dosage for this sample and haplotype
				double alternativeAlleleDosage = haplotypeDosage[haplotypeIndex];
				// assign probabilities for this sample and haplotype
				sampleProbabilities[haplotypeIndex][0] = 1 - alternativeAlleleDosage;
				sampleProbabilities[haplotypeIndex][1] = alternativeAlleleDosage;
			}
			// assign the probabilities for this sample.
			probs[sampleIndex] = sampleProbabilities;
		}
		return probs;
	}

	/**
	 * Method that converts between called alleles and phased probabilities.
	 *
	 * @param sampleAlleles A list with called alleles for every sample
	 * @param alleles The variant's alleles.
	 * @return An array with for every sample an array, with for every chromosome another array,
	 * with for every possible allele in the variant the probability of that allele being represented on the
	 * the chromosome of the sample.
	 */
	public static double[][][] convertCalledAllelesToPhasedProbabilities(List<Alleles> sampleAlleles, Alleles alleles) {
		double[][][] probs = new double[sampleAlleles.size()][][];

		if (alleles.getAlleleCount() == 0) {
			throw new GenotypeDataException("Error converting alleles to probabilities. No alleles detected");
		}

		// Loop through samples
		for (int sampleIndex = 0; sampleIndex < sampleAlleles.size(); sampleIndex++) {
			// Get the alleles that are called for this sample.
			Alleles allelesForSample = sampleAlleles.get(sampleIndex);
			// The number of alleles called is assumed to be the ploidy of the sample for this variant.
			int ploidy = allelesForSample.getAlleleCount();
			// Initialize a nested array of probabilities for the current sample.
			double[][] sampleProbabilities = new double[ploidy][alleles.getAlleleCount()];
			for (int i = 0; i < ploidy; i++) {
				// Set, at the position of the called allele in the list of alleles from the variant, the probability of 1.
				Allele allele = allelesForSample.get(i);
				int indexOfSampleAllele = alleles.getAlleles().indexOf(allele);
				if (indexOfSampleAllele != -1) {
					sampleProbabilities[i][indexOfSampleAllele] = 1;
				}
			}
			probs[sampleIndex] = sampleProbabilities;
		}

		return probs;
	}

	/**
	 * Method that converts between called alleles and complex probabilities.
	 *
	 * @param allelesPerSample The called alleles per sample.
	 * @param variantAlleles The all alleles.
	 * @return The called alleles converted to probabilities.
	 */
    public static double[][] convertCalledAllelesToComplexProbabilities(
    		List<Alleles> allelesPerSample, Alleles variantAlleles) {

		double[][] probs = new double[allelesPerSample.size()][];

		int numberOfAlleles = variantAlleles.getAlleleCount();
		if (numberOfAlleles == 0) {
			throw new GenotypeDataException("Error converting alleles to probabilities. No alleles detected");
		}

		// Loop through samples
		for (int sampleIndex = 0; sampleIndex < allelesPerSample.size(); sampleIndex++) {
			// Get the alleles that are called for this sample.
			Alleles allelesForSample = allelesPerSample.get(sampleIndex);
			// The number of alleles called is assumed to be the ploidy of the sample for this variant.
			int ploidy = allelesForSample.getAlleleCount();
			// Initialize a nested array of probabilities for the current sample.
			double[][] sampleProbabilities = new double[ploidy][numberOfAlleles];
			for (int i = 0; i < ploidy; i++) {
				// Set, at the position of the called allele in the list of alleles from the variant, the probability of 1.
				Allele allele = allelesForSample.get(i);
				int indexOfSampleAllele = variantAlleles.getAlleles().indexOf(allele);
				if (indexOfSampleAllele != -1) {
					sampleProbabilities[i][indexOfSampleAllele] = 1;
				}
			}
			probs[sampleIndex] = phasedSampleProbabilitiesToGenotypeProbabilities(
					sampleProbabilities, numberOfAlleles);
		}

		return probs;
	}
}
