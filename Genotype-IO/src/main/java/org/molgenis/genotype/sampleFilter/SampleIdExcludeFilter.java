/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.sampleFilter;

import org.molgenis.genotype.Sample;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

/**
 *
 * @author Patrick Deelen
 */
public class SampleIdExcludeFilter implements SampleFilter {

	private final HashSet<String> excludedSampleIds;

	public SampleIdExcludeFilter(Collection<String> excludedSampleIds) {
		this.excludedSampleIds = new HashSet<String>(excludedSampleIds);
	}

    public SampleIdExcludeFilter(String... ids){
		this.excludedSampleIds = new HashSet<String>();
		excludedSampleIds.addAll(Arrays.asList(ids));
	}

	public SampleIdExcludeFilter(String sample){
		this.excludedSampleIds = new HashSet<String>(1);
		excludedSampleIds.add(sample);
	}

	public SampleIdExcludeFilter(HashSet<String> excludedSampleIds) {
		this.excludedSampleIds = excludedSampleIds;
	}

	@Override
	public boolean doesSamplePassFilter(Sample sample) {
		return !excludedSampleIds.contains(sample.getId());
	}
}
