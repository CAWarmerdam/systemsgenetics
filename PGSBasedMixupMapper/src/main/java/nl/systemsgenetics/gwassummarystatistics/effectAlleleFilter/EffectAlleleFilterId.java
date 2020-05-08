/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gwassummarystatistics.effectAlleleFilter;

import nl.systemsgenetics.gwassummarystatistics.effectAllele.EffectAllele;
import org.molgenis.genotype.variant.GeneticVariant;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class EffectAlleleFilterId implements EffectAlleleFilter {

	private final Set<String> include;

	public EffectAlleleFilterId(Set<String> include) {
		if(include == null ){
			throw new IllegalArgumentException();
		}
		this.include = include;
	}

	public EffectAlleleFilterId(String... ids){
		this.include = new HashSet<>();
		include.addAll(Arrays.asList(ids));
	}

	@Override
	public boolean doesEffectAllelePassFilter(EffectAllele variant) {
		return include.contains(variant.getPrimaryVariantId());
	}
	
	public void addIdToInclude(String id){
		include.add(id);
	}

	@Override
	public boolean doesVariantIdPassFilter(String id) {
		return include.contains(id);
	}
	
}
