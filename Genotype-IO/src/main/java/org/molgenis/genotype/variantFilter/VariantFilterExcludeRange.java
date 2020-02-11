package org.molgenis.genotype.variantFilter;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VariantFilterExcludeRange implements VariantFilter {

    private Map<String, TreeSet<Pair<Integer, Integer>>> genomicRangesToExclude;
    private static Pattern GENOMIC_RANGE_PATTERN = Pattern.compile("^(\\d+):(\\d+)-(\\d+)$");

    public static VariantFilterExcludeRange fromStrings(Set<String> genomicRangesToExclude) {
        // Initialize hash set
        Map<String, TreeSet<Pair<Integer, Integer>>> rangeMap = new HashMap<>();

        for (String genomicRange : genomicRangesToExclude) {
            Matcher matcher = GENOMIC_RANGE_PATTERN.matcher(genomicRange);
            // Check if the pattern matches the input.
            // If this is the case, process the genomic range
            // If this is not the case, throw an exception.
            if (matcher.matches()) {
                // Initialize the start and end position
                int startBasePosition;
                int endBasePosition;
                try {
                    startBasePosition = Integer.parseInt(matcher.group(2));
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException(String.format(
                            "Start position %s cannot be parsed to an integer", matcher.group(2)));
                }

                try {
                    endBasePosition = Integer.parseInt(matcher.group(3));
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException(String.format(
                            "End position %s cannot be parsed to an integer", matcher.group(3)));
                }

                if (startBasePosition > endBasePosition) {
                    throw new IllegalArgumentException(String.format(
                            "Start position is greater than the end position (%d > %d)",
                            startBasePosition,
                            endBasePosition));
                }
                Pair<Integer, Integer> immutablePair = new ImmutablePair<>(
                        startBasePosition, endBasePosition);

                String sequenceName = matcher.group(1);
                if (rangeMap.containsKey(sequenceName)) {
                    rangeMap.get(sequenceName).add(immutablePair);
                } else {
                    TreeSet<Pair<Integer, Integer>> positions = new TreeSet<>();
                    positions.add(immutablePair);
                    rangeMap.put(sequenceName, positions);
                }
            } else {
                throw new IllegalArgumentException(String.format(
                        "The given string '%s' does not comply with the format for a genomic range '%s'",
                        genomicRange, GENOMIC_RANGE_PATTERN));
            }
        }

        return new VariantFilterExcludeRange(rangeMap);
    }

    private VariantFilterExcludeRange(Map<String, TreeSet<Pair<Integer, Integer>>> genomicRangesToExclude) {
        this.genomicRangesToExclude = genomicRangesToExclude;
    }

    @Override
    public boolean doesVariantPassFilter(GeneticVariant variant) {

        // Check if an excluded genomic range is on the chromosome of the given
        // variant. If this is not(!) the case the variant is not in an excluded
        // region.
        if (!genomicRangesToExclude.containsKey(variant.getSequenceName())) {
            return true;
        }

        // If this is the case, search through the ranges for the specific chromosome / sequence
        // for a range that in which the variant is located
        TreeSet<Pair<Integer, Integer>> ranges = genomicRangesToExclude.get(variant.getSequenceName());
        for (Pair<Integer, Integer> range : ranges) {
            // The variant does not pass if its start position
            // is between the left and right integer (inclusive)
            if (range.getLeft() <= variant.getStartPos()
                    && range.getRight() >= variant.getStartPos()) {
                return false;
            }
        }

        // The variant is not in a genomic range to exclude if the for loop is completed.
        return true;
    }

    @Override
    public boolean doesIdPassFilter(String id) {
        return true;
    }
}
