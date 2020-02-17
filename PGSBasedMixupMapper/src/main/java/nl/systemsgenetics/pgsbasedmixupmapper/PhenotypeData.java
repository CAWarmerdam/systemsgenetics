package nl.systemsgenetics.pgsbasedmixupmapper;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

class PhenotypeData {
    private final DoubleMatrixDataset<String, String> values;
    private final DoubleMatrixDataset<String, String> presenceMatrix;

    private PhenotypeData(DoubleMatrixDataset<String, String> values, DoubleMatrixDataset<String, String> presenceMatrix) {
        this.values = values;
        this.presenceMatrix = presenceMatrix;
    }

    /**
     * Load a CSV of phenotype values as a list of samples.
     *
     * @param inputPhenotypePath The path that points towards the file with the phenotype data.
     * @param phenotypeSampleIdentifiersToInclude A set of strings that contains all the sample identifiers to
     *                                            load phenotypes for.
     * @param delimiter The delimiter with which columns are separated.
     * @return A List of samples with the phenotypes annotated within these samples.
     * @throws IOException If an I/O error occurs while reading the phenotype data.
     */
    static PhenotypeData fromFile(File inputPhenotypePath,
                                  Set<String> phenotypeSampleIdentifiersToInclude,
                                  Set<String> traitsToInclude,
                                  char delimiter) throws IOException {

        // Create the CSV reader object with which the phenotypes can be read through
        final CSVParser parser = new CSVParserBuilder()
                .withSeparator(delimiter)
                .withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(
                new FileReader(inputPhenotypePath)))
                .withCSVParser(parser).build();

        // Create a phenotype matrix with the phenotype sample identifiers as
        DoubleMatrixDataset<String, String> phenotypeMatrix = new DoubleMatrixDataset<>(
                phenotypeSampleIdentifiersToInclude,
                traitsToInclude);

        // Create a matrix that represents the presence of phenotype values
        DoubleMatrixDataset<String, String> presenceMatrix = new DoubleMatrixDataset<>(
                phenotypeSampleIdentifiersToInclude,
                traitsToInclude);

        // Get the header containing column names.
        String[] header = reader.readNext();
        int numberOfColumns = header.length;

        // Create a list from the set of traits to include to be able to iterate through these traits
        // by using an indices
        List<String> traitsToIncludeList = new ArrayList<>(traitsToInclude);
        // Create a list from the first row (header) in the file to look up which column id corresponds to
        // which trait using .indexOf()
        List<String> headerAsList = Arrays.asList(header);

        // Get the index for every trait to include
        int[] traitIndices = new int[traitsToInclude.size()];
        for (int i = 0; i < traitsToInclude.size(); i++) {
            traitIndices[i] = headerAsList
                    .indexOf(traitsToIncludeList.get(i));
        }

        // Initialize a linked hash set to store the identifiers in to include in the analysis
        LinkedHashSet<String> ids = new LinkedHashSet<>();
        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {

            // The first columns should contain the individual id
            String individual_id = nextLine[0];

            // If the individual id is in the set indicating which samples to include,
            // include it.
            if (!phenotypeSampleIdentifiersToInclude.contains(individual_id)) {
                // Otherwise, continue on the next iteration with a new next line.
                continue;
            }

            // Check if the number of columns remains consistent
            if (numberOfColumns != nextLine.length) {
                throw new IllegalArgumentException("Different number of ids");
            }

            // Fill the phenotype matrix
            for (int traitIndex : traitIndices) {
                try {
                    String value = nextLine[traitIndex];
                    phenotypeMatrix.setElement(individual_id, header[traitIndex], Double.parseDouble(value));

                    // Set the presence for this individual id and phenotype
                    presenceMatrix.setElement(individual_id, header[traitIndex], 1);
                } catch (NumberFormatException ignored) {
                }
            }
            // The individual id should not have been added already
            if (!ids.add(individual_id)) {
                throw new IllegalArgumentException("Duplicate individual id name: " + individual_id);
            }

        }
        return new PhenotypeData(phenotypeMatrix.viewRowSelection(ids), presenceMatrix.viewRowSelection(ids));
    }

    int getSampleCount() {
        return values.rows();
    }

    List<String> getSamples() {
        return values.getRowObjects();
    }

    DoubleMatrixDataset<String, String> getValueMatrix() {
        return values;
    }

    void orderSamples(Collection<String> sampleOrder) {
        values.viewRowSelection(sampleOrder);
        presenceMatrix.viewRowSelection(sampleOrder);
    }

    DoubleMatrixDataset<String, String> getColumnNormalizedOfCompleteSamples() {
        DoubleMatrixDataset<String, String> copiedValueMatrix = values.duplicate();

        int numberOfPhenotypes = values.columns();
        Set<String> sampleIdentifiersOfPresentSamples = new LinkedHashSet<>();

        for (int i = 0; i < copiedValueMatrix.rows(); i++) {
            if (copiedValueMatrix.viewRow(i).zSum() == numberOfPhenotypes) {
                sampleIdentifiersOfPresentSamples.add(copiedValueMatrix.getRowObjects().get(i));
            }
        }

        copiedValueMatrix.viewRowSelection(sampleIdentifiersOfPresentSamples);
        copiedValueMatrix.normalizeColumns();
        return copiedValueMatrix;
    }

    private LinkedList<float[]> calculateRankedPhenotypes() {
        // Duplicate phenotype matrix to prevent ranking being reflected in the original matrix.
        DoubleMatrixDataset<String, String> processedPhenotypeMatrix = this.values.duplicate();
        for (int columnIndex = 0; columnIndex < processedPhenotypeMatrix.columns(); columnIndex++) {
            DoubleMatrix1D col = processedPhenotypeMatrix.getCol(columnIndex);
            double[] rank = rankArray(col.toArray());
            processedPhenotypeMatrix.viewCol(columnIndex).assign(rank);
        }
        return processedPhenotypeMatrix;
    }

    public int numberOfPhenotypes() {
        return values.columns();
    }

    public List<String> getPhenotypes() {
        return values.getColObjects();
    }

    public int getNumberOfSamples() {
        return values.rows();
    }
}
