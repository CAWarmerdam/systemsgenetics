/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.ArrayList;

/**
 *
 * @author patri
 */
public class PearsonRToPValueBinned {

	private final double[] pValueLookupTable;
	private final int totalNumberOfBins;
	private final int binsPerSide;
	private final double halfStep;
	private final int maxBin;

	/**
	 * It is recommended to name instances of this class r2Pvalue ;)
	 *
	 * @param numberOfBins determines precision.
	 * @param samplesUsedForCor
	 */
	public PearsonRToPValueBinned(final int numberOfBins, final int samplesUsedForCor) {

		final int df = samplesUsedForCor - 2;
		DoubleRandomEngine randomEngine = new DRand();

		this.binsPerSide = numberOfBins;
		this.totalNumberOfBins = (numberOfBins * 2) + 1;
		this.maxBin = totalNumberOfBins - 1;
		pValueLookupTable = new double[totalNumberOfBins];

		final double stepSize = 2d / totalNumberOfBins;
		halfStep = stepSize / 2;

		for (int bin = 0; bin < totalNumberOfBins; ++bin) {

			if (bin == this.binsPerSide) {

//				//Enforce r=0 results in z=0
//				zscoreLookupTable[bin] = 0;
				
				//Enforce r=0 results in p=1
				pValueLookupTable[bin] = 1;
				
			} else {

				final double corBinCenter = -1 + (stepSize * bin) + halfStep;

				StudentT tDistColt = new StudentT(df, randomEngine);
				double t = corBinCenter / (Math.sqrt((1 - corBinCenter * corBinCenter) / (double) (df)));
				double pValue;
				if (t < 0) {
					pValue = tDistColt.cdf(t);
					if (pValue < 2.0E-323) {
						pValue = 2.0E-323;
					}
				} else {
					pValue = tDistColt.cdf(-t);
					if (pValue < 2.0E-323) {
						pValue = 2.0E-323;
					}
				}

				pValueLookupTable[bin] = pValue;
				
				//System.out.println("Bin: " + bin + "\tcenter r: " + corBinCenter + "\tZscore: " + zScore);
				
			}

			

		}

	}

	/**
	 * 
	 * 
	 * @param r pearson r. An r > 1 or r small -1 is accepted for imprecise r calculations
	 * @return  P-value
	 */
	public double lookupPValueForR(double r) {

		final long bin;

		
		if (r >= 1) {
			//this is needed because due to rounding a r of 1 will not fall into any bin. 
			bin = maxBin;
		} else if (r < -1) {
			//-1 exactly will always round properly, no need to test this.
			bin = 0;
		} else {
			bin = Math.round((r + 1) * binsPerSide - halfStep);
		}
//
//		System.out.println("r: " + r);
//		System.out.println("bin: " + bin);
//		System.out.println("p: " + pValueLookupTable[(int) bin]);
//		System.out.println("---");
		return pValueLookupTable[(int) bin];

	}

	/**
	 * Inplace replacement of pearson r to P-values
	 *
	 * @param dataset
	 */
	public void inplaceRToPValue(DoubleMatrixDataset dataset) {

		DoubleMatrix2D matrix = dataset.getMatrix();

		final int rows = matrix.rows();
		final int cols = matrix.columns();

		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < cols; ++c) {
				matrix.setQuick(r, c, lookupPValueForR(matrix.getQuick(r, c)));
			}
		}

	}

}
