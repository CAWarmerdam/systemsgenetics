/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.polygenicscorecalculator;

/**
 *
 * @author MarcJan
 */
public class RiskEntry implements Comparable<RiskEntry>  {
    private final String variantId;
    private double effectSize;
    private final int startPos;
    private final int sequenceName;
    private final char allele;
    private final double pValue;

    RiskEntry(String variantId, int sequenceName, int startPos, char allele, String effectSize, double pValue) {
        this.variantId = variantId;
        this.sequenceName = sequenceName;
        this.startPos = startPos;
        this.allele = allele;
        this.pValue = pValue;
        this.effectSize = Double.parseDouble(effectSize);
    }
    
    RiskEntry(String variantId, String sequenceName, int startPos, char allele, String effectSize, double pValue) {
        this(variantId, sequenceName, startPos, allele, Double.parseDouble(effectSize), pValue);
    }

    public RiskEntry(String variantId, String sequenceName, int startPos, char allele, double effectSize, double pValue) {
        this.variantId = variantId;
        this.startPos = startPos;
        this.allele = allele;
        this.pValue = pValue;
        this.effectSize = effectSize;
        if(sequenceName.equals("X")){
            this.sequenceName = 23;
        } else if(sequenceName.equals("Y")){
            this.sequenceName = 24;
        } else {
            this.sequenceName =Integer.parseInt(sequenceName);
        }
    }

    public String getVariantId() {
        return variantId;
    }

    public double getEffectSize() {
        return effectSize;
    }

    public int getStartPos() {
        return startPos;
    }

    public String getSequenceName() {
        return String.valueOf(sequenceName);
    }

    public char getAllele() {
        return allele;
    }

    public double getpValue() {
        return pValue;
    }

    @Override
    public int compareTo(RiskEntry o) {
        if(this.pValue<o.getpValue()){
            return -1;
        } else if (this.pValue == o.getpValue()) {
            if(Math.abs(this.effectSize)>Math.abs(o.getEffectSize())){
                return -1;
            } else if (Math.abs(this.effectSize) == Math.abs(o.getEffectSize())) {
                return 0;
            } else {
                return 1;
            }
        } else {
            return 1;
        }
    }

    String InfoToString() {
        StringBuilder s = new StringBuilder();
        s.append(this.variantId).append("\t");
        s.append(this.sequenceName).append("\t");
        s.append(this.startPos).append("\t");
        s.append(this.allele).append("\t");
        s.append(this.effectSize).append("\t");
        s.append(this.pValue);
        return(s.toString());
    }

    public void setEffectSize(int i) {
        effectSize = i;
    }
}
