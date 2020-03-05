package nl.systemsgenetics.pgsbasedmixupmapper;

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
import java.util.stream.Collectors;

/**
 * Represents phenotype data.
 * This holds a matrix of doubles and a matrix of non-numeric data.
 * @author Robert Warmerdam
 */
class PhenotypeData {
    private final DoubleMatrixDataset<String, String> numericValues;
    private final DoubleMatrixDataset<String, String> presenceMatrix;
    private LinkedHashMap<String, String[]> nonNumericData;
    private final LinkedHashMap<String, Integer> nonNumericDataFieldHashMap;

    private PhenotypeData(
            DoubleMatrixDataset<String, String> numericValues,
            DoubleMatrixDataset<String, String> presenceMatrix,
            LinkedHashMap<String, String[]> nonNumericData,
            LinkedHashMap<String, Integer> nonNumericDataFieldHashMap) {
        this.numericValues = numericValues;
        this.presenceMatrix = presenceMatrix;
        this.nonNumericData = nonNumericData;
        this.nonNumericDataFieldHashMap = nonNumericDataFieldHashMap;
    }

    /**
     * Load a CSV of phenotype values as a list of samples.
     *
     * @param inputPhenotypePath The path that points towards the file with the phenotype data.
     * @param phenotypeSampleIdentifiersToInclude A set of strings that contains all the sample identifiers to
     *                                            load phenotypes for.
     * @param numericTraitsToInclude A set of traits to include from the given file.
     *                               These traits should match the column names in the file.
     * @param nonNumericTraitsToInclude Other columns to include from the given file.
     *                                  These columns should match the column names in the file.
     * @param delimiter The delimiter with which columns are separated.
     * @return A List of samples with the phenotypes annotated within these samples.
     * @throws IOException If an I/O error occurs while reading the phenotype data.
     */
    static PhenotypeData fromFile(File inputPhenotypePath,
                                  Set<String> phenotypeSampleIdentifiersToInclude,
                                  Set<String> numericTraitsToInclude,
                                  Set<String> nonNumericTraitsToInclude,
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
                numericTraitsToInclude);

        // Create a phenotype matrix
        LinkedHashMap<String, String[]> nonNumericDataMap =
                new LinkedHashMap<>();

        LinkedHashMap<String, Integer> nonNumericDataHashMap = generateLinkedHashMap(nonNumericTraitsToInclude);

        // Create a matrix that represents the presence of phenotype values
        DoubleMatrixDataset<String, String> presenceMatrix = new DoubleMatrixDataset<>(
                phenotypeSampleIdentifiersToInclude,
                numericTraitsToInclude);

        // Get the header containing column names.
        String[] header = reader.readNext();
        int numberOfColumns = header.length;
        int[] traitIndices = determineTraitIndices(numericTraitsToInclude, header);
        int[] nonNumericTraitIndices = determineTraitIndices(nonNumericTraitsToInclude, header);

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
            try {
                for (int traitIndex : traitIndices) {

                    String value = nextLine[traitIndex];
                    phenotypeMatrix.setElement(individual_id, header[traitIndex], Double.parseDouble(value));

                    // Set the presence for this individual id and phenotype
                    presenceMatrix.setElement(individual_id, header[traitIndex], 1);
                }
            } catch (NumberFormatException e) {
                System.err.println("Encountered an error: " + e.getMessage());
                System.err.println("Skipping: " + individual_id);
                continue;
            }

            try {
                // Get the non-numeric data from this line.
                String[] nonNumericData = extractNonNumericData(nonNumericDataHashMap, nonNumericTraitIndices,
                        header, nextLine);

                // Insert the non-numeric data.
                nonNumericDataMap.put(individual_id, nonNumericData);

            } catch (PhenotypeDataException e) {
                System.err.println("Encountered an error: " + e.getMessage());
                System.err.println("Skipping: " + individual_id);
                continue;
            }

            // The individual id should not have been added already
            if (!ids.add(individual_id)) {
                throw new IllegalArgumentException("Duplicate individual id name: " + individual_id);
            }

        }

        LinkedHashMap<String, String[]> filteredNonNumericData = ids.stream().collect(
                LinkedHashMap::new, (map, item) -> map.put(item, nonNumericDataMap.get(item)),
                Map::putAll);

        reportNonNumericValues(nonNumericDataHashMap, filteredNonNumericData);

        return new PhenotypeData(
                phenotypeMatrix.viewRowSelection(ids),
                presenceMatrix.viewRowSelection(ids),
                filteredNonNumericData, nonNumericDataHashMap);
    }

    private static void reportNonNumericValues(LinkedHashMap<String, Integer> nonNumericDataHashMap,
                                               LinkedHashMap<String, String[]> filteredNonNumericData) {

        System.out.printf("Read %d non-numeric fields:%n", nonNumericDataHashMap.size());

        for (Map.Entry<String, Integer> entry : nonNumericDataHashMap.entrySet()) {
            Set<String> setOfNonNumericValues = filteredNonNumericData.values().stream()
                    .map(array -> array[entry.getValue()])
                    .collect(Collectors.toSet());

            System.out.printf("\t%s: %d unique values (%s)%n",
                    entry.getKey(),
                    setOfNonNumericValues.size(),
                    String.join(", ", setOfNonNumericValues));
        }
    }

    private static String[] extractNonNumericData(LinkedHashMap<String, Integer> nonNumericDataHashMap,
                                                  int[] nonNumericTraitIndices,
                                                  String[] header, String[] nextLine) throws PhenotypeDataException {
        // Initialize an array for non numeric data for this individual.
        String[] nonNumericData = new String[nonNumericDataHashMap.size()];

        // Fill non numeric matrix
        for (int traitIndex : nonNumericTraitIndices) {
            String value = nextLine[traitIndex];
            if (value.equals("NA")) {
                throw new PhenotypeDataException(String.format(
                        "The value for '%s' equals 'NA'", header[traitIndex]));
            }
            nonNumericData[nonNumericDataHashMap.get(header[traitIndex])] = value;
        }
        return nonNumericData;
    }

    private static LinkedHashMap<String, Integer> generateLinkedHashMap(Set<String> nonNumericTraitsToInclude) {
        LinkedHashMap<String, Integer> nonNumericTraitHashMap = new LinkedHashMap<>(nonNumericTraitsToInclude.size());
        int i = 0;
        for (String trait : nonNumericTraitsToInclude) {
            nonNumericTraitHashMap.put(trait, i);
            ++i;
        }
        return nonNumericTraitHashMap;
    }

    private static int[] determineTraitIndices(Set<String> traitsToInclude, String[] header) {
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
        return traitIndices;
    }

    int getSampleCount() {
        return numericValues.rows();
    }

    List<String> getSamples() {
        return numericValues.getRowObjects();
    }

    DoubleMatrixDataset<String, String> getValueMatrix() {
        return numericValues;
    }

    void orderSamples(Collection<String> sampleOrder) {
        numericValues.viewRowSelection(sampleOrder);
        presenceMatrix.viewRowSelection(sampleOrder);
        nonNumericData = sampleOrder.stream().collect(
                LinkedHashMap::new, (map, item) -> map.put(item, nonNumericData.get(item)),
                Map::putAll);
    }

    DoubleMatrixDataset<String, String> getNormalizedPhenotypeValuesOfCompleteSamples(String groupField) {
        // First get a copy of the original data and isolate the samples that have complete records.
        LinkedHashSet<String> sampleIdentifiersOfCompleteSamples = getCompleteSamples(presenceMatrix);
        DoubleMatrixDataset<String, String> copiedValueMatrix = numericValues.duplicate()
                .viewRowSelection(sampleIdentifiersOfCompleteSamples);

        // For every group in the field to normalize separately for,
        // get sample identifiers that belong to the group
        Map<String, LinkedHashSet<String>> groupedSampleIdentifiers = new HashMap<>();

        // Loop through all samples.
        for (Map.Entry<String, String[]> entry : nonNumericData.entrySet()) {
            // Get the sample identifier
            String sampleIdentifier = entry.getKey();

            // Check if this sample is complete.
            if (sampleIdentifiersOfCompleteSamples.contains(sampleIdentifier)) {
                // Get the value this individual has for the given field.
                String groupValue = entry.getValue()[nonNumericDataFieldHashMap.get(groupField)];
                // Put a linked hash set in the grouped sample identifiers map if it is not present yet.
                groupedSampleIdentifiers.putIfAbsent(groupValue, new LinkedHashSet<>());
                // Add the sample identifier if it is not present yet.
                groupedSampleIdentifiers.get(groupValue).add(sampleIdentifier);
            }
        }

        for (String groupValue : groupedSampleIdentifiers.keySet()) {
            copiedValueMatrix.viewRowSelection(groupedSampleIdentifiers.get(groupValue)).normalizeColumns();
        }

        return copiedValueMatrix;
    }

    DoubleMatrixDataset<String, String> getNormalizedPhenotypeValuesOfCompleteSamples() {
        LinkedHashSet<String> sampleIdentifiersOfCompleteSamples = getCompleteSamples(presenceMatrix);
        return getNormalizedPhenotypeValues(sampleIdentifiersOfCompleteSamples);
    }

    private DoubleMatrixDataset<String, String> getNormalizedPhenotypeValues(LinkedHashSet<String> samples) {
        DoubleMatrixDataset<String, String> selection = numericValues.duplicate().viewRowSelection(samples);
        selection.normalizeColumns();
        return selection;
    }

    private LinkedHashSet<String> getCompleteSamples(DoubleMatrixDataset<String, String> presenceMatrix) {
        LinkedHashSet<String> sampleIdentifiersOfPresentSamples = new LinkedHashSet<>();

        for (int i = 0; i < presenceMatrix.rows(); i++) {
            if (presenceMatrix.viewRow(i).zSum() == numberOfNumericPhenotypes()) {
                sampleIdentifiersOfPresentSamples.add(presenceMatrix.getRowObjects().get(i));
            }
        }
        return sampleIdentifiersOfPresentSamples;
    }

    public int numberOfNumericPhenotypes() {
        return numericValues.columns();
    }

    public List<String> getNumericPhenotypes() {
        return numericValues.getColObjects();
    }

    public int getNumberOfSamples() {
        return numericValues.rows();
    }
}
