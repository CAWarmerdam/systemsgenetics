package nl.systemsgenetics.pgsbasedmixupmapper;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.*;

import static org.testng.Assert.*;

public class PhenotypeDataTest {

    private File examplePhenotypeFile = new File(this.getClass()
            .getResource("/phenotypedata/phenotypeswithmissingness.tsv").toURI());

    private Set<String> samples = new HashSet<>(
            Arrays.asList("S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11"));
    private List<String> expectedSamples = new ArrayList<>(
            Arrays.asList("S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09"));

    private Set<String> numericTraits = new HashSet<>(
            Arrays.asList("hdc", "ldc"));
    private List<String> expectedNumericTraits = new ArrayList<>(
            Arrays.asList("hdc", "ldc"));

    private Set<String> nonNumericData = new HashSet<>(
            Arrays.asList("sex"));

    private double[] expectedLdcValues = new double[]
            {-0.70595576, -1.47496224, 0.24192802, 0.67023177, 0.03308281, 1.40277436,  0.56280245, -0.72990141};
    private double[] expectedHdcValues = new double[]
            {-1.2425349, 1.0499510, -0.2030364, 0.4419174, 0.9101010, 0.7012716, -1.2888321, -0.3688377};

    private static final char CSV_DELIMITER = '\t';


    public PhenotypeDataTest() throws URISyntaxException {
    }

    @org.testng.annotations.Test
    public void testFromFile() throws IOException {

        PhenotypeData phenotypeData = PhenotypeData.fromFile(
                examplePhenotypeFile, samples, numericTraits, nonNumericData, CSV_DELIMITER);

        assertEquals(phenotypeData.getNumberOfSamples(), 8);

        assertEquals(phenotypeData.getSamples(), expectedSamples);

        assertEquals(phenotypeData.getNumericPhenotypes(), expectedNumericTraits);

        DoubleMatrixDataset<String, String> sexNormalizedPhenotypeData =
                phenotypeData.getNormalizedPhenotypeValuesOfCompleteSamples("sex");
        assertEquals(sexNormalizedPhenotypeData.getCol("hdc").toArray(), expectedHdcValues, 1e-7);
        assertEquals(sexNormalizedPhenotypeData.getCol("ldc").toArray(), expectedLdcValues, 1e-7);
    }
}